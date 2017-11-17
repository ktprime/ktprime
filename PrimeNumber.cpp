/*Fast Single-thread segmented sieving of prime number n < 2^64 ***/
const char* Benchmark =
"Mingw: g++ 5.1.0 bailuzhou@163\n"
"Windows 10 x64               i3-350M, i5-3470, i7-7500u i7-6700\n"
"pi(0,    1e10)  =  455052511     3.20    1.84    1.55   1.40\n"
"pi(1e11, 1e10)  =  394050419     4.45    2.56    2.05   1.90\n"
"pi(1e12, 1e10)  =  361840208     5.50    3.00    2.40   2.23\n"
"pi(1e13, 1e10)  =  334067230     6.70    3.50    2.90   2.66\n"
"pi(1e14, 1e10)  =  310208140     8.03    4.20    3.56   3.14\n"
"pi(1e15, 1e10)  =  289531946     10.3    5.10    4.40   3.78\n"
"pi(1e16, 1e10)  =  271425366     12.4    6.10    5.11   4.52\n"
"pi(1e17, 1e10)  =  255481287     15.0    7.09    5.95   5.31\n"
"pi(1e18, 1e10)  =  241272176     18.1    8.58    7.17   6.41\n"
"pi(1e19, 1e10)  =  228568014     25.2    11.6    9.86   8.72\n"
"pi(0-1e9,1e9)   =  22537866      8.75    4.28    3.92   3.62\n"
"pi(1e18, 1e6)   =  24280         0.75    0.46    0.34   0.50\n"
"pi(1e18, 1e8)   =  2414886       1.50    0.81    0.70   0.64\n"
"pi(1e18, 1e9)   =  24217085      3.70    1.80    1.48   1.40\n"
"pi(0,    1e12)  =  37607912018   520     280     234    210\n"
"pi(1e14, 1e12)  =  31016203073   810     430     354    314\n"
"pi(1e16, 1e12)  =  27143405794   1220    600     512    445\n"
"pi(1e18, 1e12)  =  24127637783   1570    760     622    590\n"
"pi(1e19, 1e12)  =  22857444126   1730    840     702    625\n"
"pi(1e18, 1e13)  =  241274558866  15700   7600    6200   5900\n"
"pi(1e19, 1e13)  =  228575545410  17300   8400    7010   6230\n";

#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <string.h>

//const
enum ECONST
{
	WHEEL30       = 30,
	WHEEL210      = 210,
	PRIME_PRODUCT = 210 * 11 * 13 * 17 * 19,
	FIRST_INDEX   = PRIME_PRODUCT / 9699690 + 7,
	WHEEL_SKIP    = 0x5A28A6,
	PATTERN_GAP   = 10,//max prime wheel 210 gap
	PRIME_SIZE    = 6542 + 10,//pi(2^16)
	MAX_SEGMENT   = 4 << 10,//4M max sieve_size
};

//performance marco
enum ECONFIG
{
	WP_BIT         = 6, //48 < 2^WP_BIT < WHEEL210
	SI_BIT         = 8, //8 -- 10
#if (L2_DCACHE_SIZE < 256 || L2_DCACHE_SIZE > 4096)
	L2_DCACHE_SIZE =  256,
#endif
#ifndef SIEVE_SIZE
	SIEVE_SIZE     =  2048,
#endif
#ifndef UINT_MAX
	UINT_MAX       = 0-1u,
#endif

	ERAT_SMALL     = 2, //4 - 16
	ERAT_MEDIUM    = 2, //2 - 6
};

enum EBUCKET
{
	UINT_PIMAX = 203280221, //= pi(2^32)
	MAX_BUCKET = 0xffffffff / (256 * (1 << 10) * 3) + 4, //5465 (sqrtp * (8 + 2) / sieve_size + 2);
	WHEEL_SIZE = 1 << 11, //=4096 11: 16k, [10 - 13]
	MEM_WHEEL  = WHEEL_SIZE * 8, //=32768
	MEM_BLOCK  = (1 << 19) / WHEEL_SIZE, //=256 1 << 20:8 MB
	MAX_STOCK  = UINT_PIMAX / WHEEL_SIZE + MAX_BUCKET + MEM_BLOCK, //=55221
	MAX_POOL   = UINT_PIMAX / (MEM_BLOCK * WHEEL_SIZE) + 100, //=487
};

enum EFLAG
{
	PRINT_TIME = 1 << ('T' - 'A'),
	SLOW_TEST  = 1 << ('A' - 'A'),
	SAVE_DATA  = 1 << ('F' - 'A'),
	PRINT_RET  = 1 << ('R' - 'A'),
};

enum ECMD
{
	COUNT_PRIME = 0,
	COPY_BITS,
	SAVE_BYTE,
	SAVE_PRIME,
	PCALL_BACK,
	SAVE_BYTEGAP,
	FIND_MAXGAP,
};

enum EBITMASK
{
	BIT0 = 1 << 0,
	BIT1 = 1 << 1,
	BIT2 = 1 << 2,
	BIT3 = 1 << 3,
	BIT4 = 1 << 4,
	BIT5 = 1 << 5,
	BIT6 = 1 << 6,
	BIT7 = 1 << 7,
};

#ifndef ERAT_BIG
# define ERAT_BIG         5 //2 - 6
#endif

#ifdef W30 //fast on some cpu i7
	# define WHEEL        WHEEL30
	# define WHEEL_MAP    Wheel30
	# define WHEEL_INIT   WheelInit30
	# define WHEEL_FIRST  WheelFirst30
	# define PATTERN      Pattern30
#else
	# define WHEEL        WHEEL210
	# define WHEEL_MAP    Wheel210
	# define WHEEL_INIT   WheelInit210
	# define WHEEL_FIRST  WheelFirst210
	# define PATTERN      Pattern210
#endif

#if __x86_64__ || _M_AMD64 || __amd64__
	# define X86_64       1
#endif

#if X86_64 && _MSC_VER
	# define ASM_X86      0
	# define BIT_SCANF    1
	#include<intrin.h>
#elif X86_64 || _M_IX86 || __i386__
	# define ASM_X86      1
	# define BIT_SCANF    1
#else
	# define ASM_X86      0
	# define BIT_SCANF    0
#endif

#if _WIN32 && _MSC_VER < 1500
	typedef unsigned __int64 uint64;
	typedef __int64 int64;
#else
	typedef unsigned long long uint64;
	typedef long long int64;
#endif

typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;

#define MIN(a, b)         (a < b ? a : b)
#define MAX(a, b)         (a > b ? a : b)

#if BIT_SCANF == 0
	#define PRIME_OFFSET(mask) Lsb[mask]
	typedef ushort stype;
#else
	#define PRIME_OFFSET(mask) Pattern30[bitScanForward(mask)]
#if X86_64
	typedef uint64 stype;
#else
	typedef uint stype;
#endif
#endif

static const char* Help = "\
	[B: Benchmark (0 - 12, 0 - 40)]\n\
	[D: D[T, R] dump time and result]\n\
	[M: Progress of calculating (0 - 20)]\n\
	[C: Cpu L1/L2 data cache size (L1:16-64, L2:256-4096)k]\n\
	[S: Set sieve segment size (32 - 4096)k]\n\
	[L: Set sieve cache segs L(2-12)1, L(2-8)2 L(2-6)3]\n\
	[I: Info of programming]\n\
Example:\n\
	1e16 10^10 s1024 c321 c2562";

static struct _Threshold
{
	uint L1Size; //cpu L1/L2 cache size
	uint L2Size;

	uint L1Maxp;
	uint L1Index;
	uint L1Segs;

	uint L2Maxp;
	uint L2Index;
	uint L2Segs;

	uint Medium;
	uint Msegs;

	uint64 BucketStart; //min bucket start
}
Threshold =
{
	32 << 10, L2_DCACHE_SIZE << 10,
	(32 << 10) / ERAT_SMALL, 0, ERAT_SMALL,
	(L2_DCACHE_SIZE << 10) / ERAT_MEDIUM, 0, ERAT_MEDIUM,
	SIEVE_SIZE * (WHEEL30 << 10) / ERAT_BIG, ERAT_BIG,
};

struct _Config
{
	uint SieveSize;
	uint Flag;
	uint Progress;
};

_Config Config =
{
	SIEVE_SIZE << 10,
	PRINT_RET | PRINT_TIME,
	(1 << 5) - 1,
};

struct WheelPrime
{
	//[0 - 6]: p % wheel, [7 - 31]: p / sieve_size
	uint Wp;
	//[0 - 7]: index % sieve_size, [8 - 31]: index / WHEEL30
	uint Si;
};

struct Stock
{
	WheelPrime* RWheel; //read only 4k size
	Stock* Next;
};

struct _Bucket
{
	WheelPrime* Wheel; //write only
	Stock* Head; //Head->Stock1->Stock2->....->Stockn
};

struct _BucketInfo
{
	uint CurStock;
	uint MaxBucket; //overflow for big range
	uint LoopSize;
	uint SieveSize;
	uint Log2Size;

	uint StockSize;
	uint PoolSize;
};

//thread ...
static _BucketInfo BucketInfo;
static _Bucket Bucket [MAX_BUCKET];
static Stock StockCache [MAX_STOCK];
static WheelPrime* WheelPool [MAX_POOL]; //2G vm

static Stock* StockHead;
static WheelPrime* MediumSieve;
static WheelPrime* SmallSieve;

//small prime
static uint Prime[PRIME_SIZE];

//presieved with prime <= 19.
static uchar PreSieved[PRIME_PRODUCT / WHEEL30];

#if BIT_SCANF == 0
static uchar Lsb[1 << 16];
#endif

//number of bits 1 binary representation table in Range[0-2^16)
static uchar WordNumBit1[1 << 16];

struct WheelElement
{
	char WheelIndex;
	uchar UnsetBit;
	uchar Correct;
	uchar NextMultiple;
};

struct WheelInit
{
	char WheelIndex;
	uchar UnsetBit;
	uchar PrimeIndex;
	uchar Reserved; //remove slow ?
};

typedef WheelElement WheelFirst;
static WheelInit WheelInit30[WHEEL30];
static WheelFirst WheelFirst30[WHEEL30][8];
static WheelElement Wheel30[8][8];

static WheelInit WheelInit210[WHEEL210];
static WheelFirst WheelFirst210[WHEEL210][48];
static WheelElement Wheel210[48][48];

static uchar Pattern30[64];
static uchar Pattern210[WHEEL210];

//api
struct Cmd
{
	int Oper;
	uchar* Data;
	uint64 Primes;
};

typedef void (*sieve_call)(uint64, uint64);
static void printInfo();
void initPrime(int sieve_size);
uint64 doSieve(const uint64 start, const uint64 end, Cmd* cmd);
int setSieveSize(uint sieve_size);
void setLevelSegs(uint level, uint size);
void setCpuCache(uint level, uint cpusize);

static int64 getTime()
{
#ifdef _WIN32
	return clock();
#else
	return clock() / 1000;
#endif
}

//n > 1
static int ilog(uint64 x, const uint n)
{
	int logn = 0;
	while (x / n) {
		logn ++;
		x /= n;
	}

	return logn;
}

//x^n < 2^64, n < 64
static uint64 ipow(const uint x, uint n)
{
	uint64 pown = 1;
	while (n --) {
		pown *= x;
	}

	return pown;
}

//x > 1
static uint isqrt(const uint64 x)
{
	const uint s = ilog(x - 1, 2);
	uint64 g0 = (uint64)1 << s;
	uint64 g1 = (g0 + (x >> s)) >> 1;

	while (g1 < g0) {
		g0 = g1;
		g1 = (g1 + x / g1) >> 1;
	}

	return (uint)g0;
}

#if BIT_SCANF
inline static uint bitScanForward(const stype n)
{
#if _MSC_VER
	unsigned long index;
	#if X86_64
	_BitScanForward64(&index, n);
	#else
	__asm
	{
		bsf eax, n
		mov index, eax
	}
	#endif
	return index;
#else
	#if X86_64
	return __builtin_ffsll(n) - 1;
	#else
	return __builtin_ffsl(n) - 1;
	#endif
#endif
}
#endif

//n % p < 2^32
inline static uint fastMod(const uint64 n, uint p)
{
#if ASM_X86 == 0
	p = (uint)(n % p);
#elif __GNUC__
	uint loww = n, higw = n >> 32;
	__asm
	(
		"divl %%ecx\n"
		: "=d" (p)
		: "d"(higw), "a"(loww), "c"(p)
	);
#else
	uint loww = (uint)n, higw = (uint)(n >> 32);
	__asm
	{
		mov eax, loww
		mov edx, higw
		div p
		mov p, edx
	}
#endif

	return p;
}

//valid is [number][e][+-*][number]
//ex: 1000, 2e9, 300-2e1, 2^32*2-10^3, 2^30-1E2, 2e9+2^20
uint64 atoint64(const char* str)
{
	uint64 ret = 0;
	while (isdigit(*str)) {
		ret = ret * 10 + *str ++ - '0';
	}

	if (*str && isdigit(str[1])) {
		if (str[0] == '^') {
			ret = ipow((uint)ret, atoi(str + 1));
		} else if (str[0] == 'e' || str[0] == 'E') {
			if (ret == 0) {
				ret = 1;
			}
			ret *= ipow(10, atoi(str + 1));
		}
	}

	const char* ps = str;
	if ((ps = strchr(str, '+')) != NULL) {
		ret += atoint64(ps + 1);
	} else if ((ps = strchr(str, '-')) != NULL) {
		ret -= atoint64(ps + 1);
	} else if ((ps = strchr(str, '*')) != NULL) {
		ret *= atoi(ps + 1);
	} else if ((ps = strchr(str, '/')) != NULL) {
		ret /= atoi(ps + 1);
	}

	return ret;
}

//count number of bit 0 in binary representation of array
static uint countBit0sArray(const uint64 bitarray[], const uint bitsize)
{
	int loops = bitsize / 64;
	uint bit0s = (1 + loops) * 64;

	while (loops -- >= 0) {
		const uint hig = (uint)(*bitarray >> 32);
		const uint low = (uint)*bitarray++;
		bit0s -= WordNumBit1[(ushort)hig] + WordNumBit1[hig >> 16];
		bit0s -= WordNumBit1[(ushort)low] + WordNumBit1[low >> 16];
	}

	return bit0s;
}

static void allocWheelBlock(const uint blocks)
{
	WheelPrime *pwheel = (WheelPrime*) malloc((MEM_BLOCK + 1) * MEM_WHEEL);
	WheelPool[BucketInfo.PoolSize ++] = pwheel;
	//assert (BucketInfo.PoolSize < sizeof(WheelPool) / sizeof(WheelPool[0]));
	//assert (BucketInfo.StockSize + blocks < sizeof(StockCache) / sizeof(StockCache[0]));

	//align by MEM_WHEEL
	pwheel = (WheelPrime*)((size_t)pwheel + MEM_WHEEL - (size_t)pwheel % MEM_WHEEL);
	Stock* pStock = StockCache + BucketInfo.StockSize;
	for (uint i = 0; i < MEM_BLOCK; i ++) {
		pStock->RWheel = pwheel + WHEEL_SIZE * i;
		pStock->Next = pStock + 1;
		pStock ++;
	}
	pStock[-1].Next = StockHead;
	StockHead = pStock - MEM_BLOCK;

	BucketInfo.StockSize += MEM_BLOCK;
	BucketInfo.CurStock += MEM_BLOCK;
}

#define PI(x, r) x / log((double)x) * (1 + r / log((double)x))
static int setBucketInfo(const uint sieve_size, const uint sqrtp, const uint64 range)
{
	//assert (range >> 32 < sieve_size);
	BucketInfo.MaxBucket = range / sieve_size + 1;
	BucketInfo.LoopSize = sqrtp / (sieve_size / PATTERN_GAP) + 2; //add the first, last bucket
	BucketInfo.Log2Size = ilog(sieve_size / WHEEL30, 2);
	BucketInfo.SieveSize = (1 << BucketInfo.Log2Size) - 1;
	//assert (BucketInfo.Log2Size <= (32 - SI_BIT) && SI_BIT >= WP_BIT);
	return 0;
}

static int setWheelSmall(const uint64 start, uint l1_size = Threshold.L1Size * WHEEL30)
{
	const uint align210 = start % WHEEL210;
	uint64 offset = start - align210;

	for (uint p = Prime[FIRST_INDEX], j = FIRST_INDEX; p < Threshold.L1Maxp; p = Prime[++j]) {
		const uint64 p2 = (uint64)p * p;
		if (p2 > offset && p2 > offset + l1_size) { //overflow
			offset += (p2 - offset) / l1_size * l1_size;
		}
		//assert(p2 < offset);
		uint sieve_index = p - (uint)(offset % p);
		if (p2 > offset) {
			sieve_index = (uint)(p2 - offset);
		}

		const WheelFirst& wf = WheelFirst30[sieve_index % WHEEL30][WheelInit30[p % WHEEL30].PrimeIndex];
		const uint multiples = (WHEEL_SKIP >> (wf.NextMultiple * 3)) | (WHEEL_SKIP << (24 - wf.NextMultiple * 3));
		sieve_index += wf.Correct * p;

		SmallSieve[j].Wp = multiples;
#ifndef SM0
		SmallSieve[j].Si = 0 - sieve_index;
#else
		SmallSieve[j].Si = sieve_index;
#endif
	}

	return align210;
}

static int segmentedSieve2(uchar bitarray[], uint start, uint sieve_size);
static void setWheelMedium(uchar* bitarray, const uint sieve_size, const uint medium, const uint64 start)
{
	uint j = Threshold.L1Index, l1_maxp = Threshold.L1Maxp - Threshold.L1Maxp % WHEEL210;
	Threshold.L2Index = 0;
	uint msize = MIN(sieve_size, Threshold.L2Size * WHEEL30);
	const int bytes = segmentedSieve2(bitarray, l1_maxp, medium - l1_maxp);

	stype mask = 0, *sbitarray = (stype*)bitarray;
	uint noffset = l1_maxp - sizeof(mask) * WHEEL30;
	const uint pn = 2 + bytes / sizeof(mask);

	for (uint i = 0; i < pn; ) {
		if (mask == 0) {
			mask = ~sbitarray[i ++];
			noffset += sizeof(mask) * WHEEL30;
			continue;
		}

		const uint p = noffset + PRIME_OFFSET(mask); mask &= mask - 1;
		if (p >= Threshold.L2Maxp && Threshold.L2Index == 0 /*&& p % WHEEL30 == 1*/) {
			Threshold.L2Maxp = p;
			Threshold.L2Index = j;
			msize = sieve_size;
		}

		uint64 offset = start;
		const uint64 p2 = (uint64)p * p;
		if (p2 - msize >= offset) { //overflow
			offset += (p2 - offset) / msize * msize;
		}

		//assert (p2 < offset + msize);
		uint sieve_index = p - (uint)(offset % p);
		if (p2 > offset) {
			sieve_index = (uint)(p2 - offset);
		}

		const int pi = WHEEL_INIT[p % WHEEL].PrimeIndex;
		const WheelFirst& wf = WHEEL_FIRST[(sieve_index + uint(offset % WHEEL)) % WHEEL][pi];
		sieve_index += wf.Correct * p;
		assert (sieve_index / WHEEL30 < (-1u >> SI_BIT)); //TODO
		MediumSieve[j].Wp = (p / WHEEL << SI_BIT) + pi;
		MediumSieve[j ++].Si = (sieve_index / WHEEL30 << SI_BIT) + wf.WheelIndex;
	}

	MediumSieve[j + 0].Wp = UINT_MAX;
	MediumSieve[j + 1].Wp = UINT_MAX;
}

static void pushBucket(const uint offset, const uint wp, const uchar wi)
{
	const uint next_bucket = offset >> BucketInfo.Log2Size;
#ifndef B_R
	if (next_bucket >= BucketInfo.MaxBucket) {
		return;
	}
#endif

	_Bucket* pbucket = Bucket + next_bucket;
	WheelPrime* wheel = pbucket->Wheel ++;
	if ((uint)(size_t)(wheel) % MEM_WHEEL == 0) {
		BucketInfo.CurStock --;
//		if (BucketInfo.CurStock -- == 0) allocWheelBlock(1);
		Stock* phead = StockHead; StockHead = phead->Next;
		phead->Next = pbucket->Head;
		pbucket->Head = phead;
		wheel = phead->RWheel;
		pbucket->Wheel = wheel + 1;
	}

	wheel->Wp = wp;
	wheel->Si = (offset & BucketInfo.SieveSize) << SI_BIT | wi;
}

static void setWheelBig(uchar* bitarray, uint medium, uint sqrtp, const uint64 start, const uint64 range)
{
	uint nextp = 0; uint64 remp = 0;
	if (sqrtp < UINT_MAX) sqrtp ++; //watch overflow if sqrtp = 2^32 - 1

	const uint irange = (uint)((range >> 32) > WHEEL30 ? UINT_MAX : range / WHEEL30);
	for (uint l2size = L2_DCACHE_SIZE * WHEEL30 << 10; medium < sqrtp; medium += l2size) {
		if (l2size > sqrtp - medium)
			l2size = sqrtp - medium;

		const int bytes = segmentedSieve2(bitarray, medium, l2size);
		stype mask = 0, *sbitarray = (stype*)bitarray; //little endian
		uint offset = medium - sizeof(mask) * WHEEL30;
		const uint pn = 2 + bytes / sizeof(mask);

		for (uint j = 0; j < pn; ) {
			if (mask == 0) {
				mask = ~sbitarray[j ++];
				offset += sizeof(mask) * WHEEL30;
				continue;
			}

			const uint p = offset + PRIME_OFFSET(mask); mask &= mask - 1;
#if ASM_X86
			if (p > nextp) { //ugly & difficult to understand but efficient
				remp = start / (nextp = p + (uint64)p * p / (uint)(start >> 32));
				if (p > nextp) //overflow
					remp = start >> 32, nextp = UINT_MAX;
			}
			uint sieve_index = p - fastMod(start - remp * p, p);
#else
			uint sieve_index = p - (uint)(start % p);
#endif

#if 0
			if (sieve_index > range) continue;
#endif
			const uint wp = (p / WHEEL210 << WP_BIT) + WheelInit210[p % WHEEL210].PrimeIndex;
			const uint module_wheel = sieve_index % WHEEL210;
			const WheelFirst& wf = WheelFirst210[module_wheel][wp % (1 << WP_BIT)];
#if X86_64
			sieve_index = (sieve_index + wf.Correct * (uint64)p) / WHEEL30;
#else
			sieve_index = (sieve_index / WHEEL210 + (wp >> WP_BIT) * wf.Correct) * (WHEEL210 / WHEEL30);
			sieve_index += (wf.Correct * (p % WHEEL210) + module_wheel) / WHEEL30;
#endif

#ifndef B_R
			if (sieve_index > irange) continue;
#endif
	//		pushBucket(sieve_index, wp, wf.WheelIndex);

			_Bucket* pbucket = Bucket + (sieve_index >> BucketInfo.Log2Size);
			WheelPrime* wheel = pbucket->Wheel ++;
			if ((uint)(size_t)(wheel) % MEM_WHEEL == 0) {
				if (BucketInfo.CurStock -- == 0) allocWheelBlock(1);
				Stock* phead = StockHead; StockHead = phead->Next;
				phead->Next = pbucket->Head;
				pbucket->Head = phead;
				wheel = phead->RWheel;
				pbucket->Wheel = wheel + 1;
			}
			wheel->Wp = wp;
			wheel->Si = (sieve_index & BucketInfo.SieveSize) << SI_BIT | wf.WheelIndex;
		}
	}
}

static int sieveSmall0(uchar bitarray[], const uchar* pend, const uint p, uint offset, uint multiples)
{
	uchar* ppbeg[8], wi;
	for (int i = 0; i < 8; i ++) {
		wi = WheelInit30[offset % WHEEL30].WheelIndex;
		ppbeg[wi] = bitarray + offset / WHEEL30;
		offset += (multiples % 8) * p; multiples /= 8;
	}

	#define OR_ADD(n) *ps##n |= BIT##n, ps##n += p
	uchar* ps0 = ppbeg[0], *ps1 = ppbeg[1], *ps2 = ppbeg[2], *ps3 = ppbeg[3];
	uchar* ps4 = ppbeg[4], *ps5 = ppbeg[5], *ps6 = ppbeg[6], *ps7 = ppbeg[7];

	uchar* pmax = ppbeg[wi] - p;
	while (pmax <= pend) {
		pmax += p;
		OR_ADD(0); OR_ADD(5); OR_ADD(1); OR_ADD(4);
		OR_ADD(3); OR_ADD(6); OR_ADD(2); OR_ADD(7);
	}
	return (int)(pmax - ppbeg[wi]) * WHEEL30 + offset;
}

static void sieveSmall1(uchar bitarray[], const uchar* pend, const uint p, uint offset, uint multiples)
{
	for (int i = 0; i < 4; i ++) {
		uchar* ps0 = bitarray + offset / WHEEL30;
		const uchar masks0 = WheelInit30[offset % WHEEL30].UnsetBit;
		offset += (multiples % 8) * p; multiples /= 8;

		uchar* ps1 = bitarray + offset / WHEEL30;
		const uchar masks1 = WheelInit30[offset % WHEEL30].UnsetBit;
		offset += (multiples % 8) * p; multiples /= 8;

		while (ps1 <= pend) {
			*ps1 |= masks1, ps1 += p;
			*ps0 |= masks0, ps0 += p;
		}
		if (ps0 <= pend)
			*ps0 |= masks0;
	}
}

#define poSet(s1,o1,b1,s2,o2,b2,s3,o3,b3,s4,o4,b4,s5,o5,b5,s6,o6,b6,s7,o7,b7)\
	while (ps <= pend) {\
		ps[o * 00 + 00] |= BIT##0,  ps[o * s1 + o1] |= BIT##b1;\
		ps[o * s2 + o2] |= BIT##b2, ps[o * s3 + o3] |= BIT##b3;\
		ps[o * s4 + o4] |= BIT##b4, ps[o * s5 + o5] |= BIT##b5;\
		ps[o * s6 + o6] |= BIT##b6, ps[o * s7 + o7] |= BIT##b7;\
		ps += p;\
	}

static int sieveSmall3(uchar* ps, const uchar* pend, const uint p)
{
	const uchar* pbeg = ps;
	const uint o = p / WHEEL30;
	switch (WheelInit30[p % WHEEL30].PrimeIndex)
	{
		case 3 : poSet(4,1,6, 6,2,5,  10,4,2,  12,5,1,  16,6,7,  22,9,4,  24,10,3) break;
		case 7 : poSet(2,1,7, 8,7,6,  12,11,5, 14,13,4, 18,17,3, 20,19,2, 24,23,1) break;
		case 2 : poSet(2,0,6, 6,2,1,  8,2,7,   12,4,3,  18,6,5,  20,7,2,  26,9, 4) break;
		case 1 : poSet(4,0,7, 6,1,3,  10,2,2,  16,3,6,  18,4,1,  24,5,5,  28,6, 4) break;
		case 6 : poSet(2,1,4, 6,4,5,  12,9,1,  14,10,6, 20,15,2, 24,18,3, 26,19,7) break;
		case 4 : poSet(6,3,3, 8,4,4,  14,7,7,  18,10,1, 20,11,2, 24,13,5, 26,14,6) break;
		case 0 : poSet(6,0,1, 10,0,2, 12,0,3,  16,0,4,  18,0,5,  22,0,6,  28,0, 7) break;
		case 5 : poSet(4,2,4, 10,6,2, 12,7,5,  18,11,3, 22,13,7, 24,15,1, 28,17,6) break;
	}

	return (int)(ps - pbeg) * WHEEL30;
}

static void preSieve(uchar bitarray[], const uint64 start, const uint sieve_size)
{
	const uint offset = (uint)(start % PRIME_PRODUCT) / WHEEL30;
	const uint bits = sieve_size / WHEEL30 * 8 + WheelInit30[sieve_size % WHEEL30].WheelIndex;
	uint bytes = (bits + 7) / 8, remains = sizeof(PreSieved) - offset;
	if (remains > bytes) {
		memcpy(bitarray, PreSieved + offset, bytes);
	} else {
		memcpy(bitarray, PreSieved + offset, remains);
		for (; remains + sizeof(PreSieved) < bytes; remains += sizeof(PreSieved))
			memcpy(bitarray + remains, PreSieved, sizeof(PreSieved));
		memcpy(bitarray + remains, PreSieved, bytes - remains);
	}

	//1 is not prime, pattern < WHEEL30 is prime
	if (start == 0) {
		bitarray[0] = BIT0;
	}

	//fill the last byte with bit 1
	if (bits % 8 != 0)
		bitarray[bits >> 3] |= ~((1 << (bits & 7)) - 1);
}

static void eratSieveL0(uchar bitarray[], const uint64 start, const uint sieve_size, uint maxp)
{
	const uchar* pend = bitarray + sieve_size / WHEEL30;
	if (start + sieve_size < ((uint64)maxp) * maxp) {
		maxp = isqrt(start + sieve_size) + 1;
	}

	for (uint p = Prime[FIRST_INDEX], j = FIRST_INDEX; p < maxp; p = Prime[++j]) {
		uint offset = p - (uint)(start % p);
		if (start <= p) {
			offset = p * p - (uint)start;
		}

		const WheelFirst& wf = WheelFirst30[offset % WHEEL30][WheelInit30[p % WHEEL30].PrimeIndex];
		const uint multiples = (WHEEL_SKIP >> (wf.NextMultiple * 3)) | (WHEEL_SKIP << (24 - wf.NextMultiple * 3));
		offset += wf.Correct * p;
		sieveSmall0(bitarray, pend, p, offset, multiples);
	}
}

static void eratSieveL1(uchar bitarray[], const uint64 start, const uint sieve_size, uint maxp)
{
	const uchar* pend = bitarray + sieve_size / WHEEL30;
	if (start + sieve_size < ((uint64)maxp) * maxp) {
		maxp = isqrt(start + sieve_size) + 1;
	}

	for (uint p = Prime[FIRST_INDEX], j = FIRST_INDEX; p < maxp; p = Prime[++j]) {
		uint& offset = SmallSieve[j].Si;
		const uint multiples = SmallSieve[j].Wp;
#ifndef SM0
		if ((int)offset <= 0) {
			offset = 0 - offset;
			for (int i = 0; i < 24; i += 3) {
				const uchar mask = WheelInit30[offset % WHEEL30].UnsetBit;
				if (mask == BIT0) {
					break;
				}
				bitarray[offset / WHEEL30] |= mask;
				offset += (multiples >> i) % 8 * p;
			}
		}
		offset += sieveSmall3(bitarray + offset / WHEEL30, pend, p) - sieve_size;
#else
		offset = sieveSmall0(bitarray, pend, p, offset, multiples) - sieve_size;
#endif
	}
}

static void eratSieveSmall(uchar bitarray[], const uint64 start, const uint sieve_size)
{
	for (uint offset = 0, l1size = Threshold.L1Size * WHEEL30; offset < sieve_size; offset += l1size) {
		if (l1size + offset > sieve_size)
			l1size = sieve_size - offset;
		eratSieveL1(bitarray + offset / WHEEL30, start + offset, l1size, Threshold.L1Maxp);
	}
}

#define SAFE_SET(n) \
	we##n = wdata##n[(int)we##n.WheelIndex]; \
	bitarray[(int)offset##n] |= we##n.UnsetBit; \
	offset##n += we##n.Correct + (we##n.NextMultiple) * wi##n

//sieve 1 medium prime from array
static void sieveMedium1(uchar bitarray[], const uint sieve_byte, WheelPrime* pwheel)
{
	uint& wheel = pwheel->Si;
	int offset = (wheel >> SI_BIT) - sieve_byte;
	if (offset >= 0) {
		wheel -= sieve_byte << SI_BIT;
		return;
	}

	const uint wi = pwheel->Wp >> SI_BIT;
	WheelElement* wdata = WHEEL_MAP[pwheel->Wp % (1 << SI_BIT)];
	WheelElement we; we.WheelIndex = wheel % (1 << SI_BIT);

	do {
		we = wdata[we.WheelIndex];
		bitarray[(int)offset] |= we.UnsetBit;
		offset += we.Correct + (we.NextMultiple) * wi;
	} while (offset < 0);

	wheel = (offset << SI_BIT) | we.WheelIndex;
}

//sieve 3 medium prime from array
inline static int sieveMedium3(uchar bitarray[], const uint sieve_byte, WheelPrime* pwheel)
{
	uint& wheel0 = pwheel[0].Si, wp0 = pwheel[0].Wp;
	uint wi0 = wp0 >> SI_BIT, offset0 = (wheel0 >> SI_BIT) - sieve_byte;

	uint& wheel1 = pwheel[1].Si, wp1 = pwheel[1].Wp;
	uint wi1 = wp1 >> SI_BIT, offset1 = (wheel1 >> SI_BIT) - sieve_byte;

	uint& wheel2 = pwheel[2].Si, wp2 = pwheel[2].Wp;
	uint wi2 = wp2 >> SI_BIT, offset2 = (wheel2 >> SI_BIT) - sieve_byte;

	WheelElement* wdata0 = WHEEL_MAP[wp0 % (1 << SI_BIT)];
	WheelElement* wdata1 = WHEEL_MAP[wp1 % (1 << SI_BIT)];
	WheelElement* wdata2 = WHEEL_MAP[wp2 % (1 << SI_BIT)];

	WheelElement we0, we1, we2;
	we0.WheelIndex = wheel0 % (1 << SI_BIT);
	we1.WheelIndex = wheel1 % (1 << SI_BIT);
	we2.WheelIndex = wheel2 % (1 << SI_BIT);

	while ((int)offset0 < 0) {
		SAFE_SET(0);
#if 1
		if ((int)offset1 >= 0) break;
		SAFE_SET(1);
		if ((int)offset2 >= 0) break;
		SAFE_SET(2);
#else
		if ((int)offset1 < 0) { SAFE_SET(1); }
		if ((int)offset2 < 0) { SAFE_SET(2); }
#endif
	}

	while ((int)offset0 < 0) { SAFE_SET(0); }
	while ((int)offset1 < 0) { SAFE_SET(1); }
	while ((int)offset2 < 0) { SAFE_SET(2); }

	wheel0 = offset0 << SI_BIT | we0.WheelIndex;
	wheel1 = offset1 << SI_BIT | we1.WheelIndex;
	wheel2 = offset2 << SI_BIT | we2.WheelIndex;

	return 3;
}

#if W30
static WheelPrime* sieveMediumW30(uchar bitarray[], const uint sieve_byte, const uint minwp, WheelPrime* pwheel)
{
	for (uint wp = pwheel->Wp; wp < minwp; wp = pwheel->Wp) {
		const uint wi = wp >> SI_BIT, pi = wp % (1 << SI_BIT);
		const uint p = wi * WHEEL + PATTERN[pi];
		WheelElement* wdata = WHEEL_MAP[pi];
		WheelElement* wheel = wdata + pwheel->Si % (1 << SI_BIT);

		uchar* ps = bitarray + (pwheel->Si >> SI_BIT);
		uchar* pend = bitarray + sieve_byte;
		uchar* pmin = pend + p;

		int nwi, nwa[8];
		for (int i = 0; i < 8; i += 2) {
			uchar* ps0 = ps;
			uchar masks0 = wheel->UnsetBit;
			ps += wheel->Correct + wheel->NextMultiple * wi;
			nwa[i + 0] = wheel->WheelIndex;
			wheel = wdata + wheel->WheelIndex;

			nwa[i + 1] = wheel->WheelIndex;
			uchar* ps1 = ps;
			uchar masks1 = wheel->UnsetBit;
			ps += wheel->Correct + wheel->NextMultiple * wi;
			wheel = wdata + wheel->WheelIndex;
			while (ps1 < pend) {
				*ps1 |= masks1, ps1 += p;
				*ps0 |= masks0, ps0 += p;
			}
			if (ps0 < pend) {
				*ps0 |= masks0, ps0 += p;
			}

			if (ps1 < pmin)
				pmin = ps1, nwi = i;
			if (ps0 < pmin)
				pmin = ps0, nwi = (i + 7) % 8;
		}
		pwheel ++->Si = ((pmin - pend) << SI_BIT) | nwa[nwi];
	}

	return pwheel;
}
#endif

//sieve prime multiples in [start, start + sieve_size) with medium algorithm
static void eratSieveMedium(uchar bitarray[], const uint64 start, const uint sieve_size, const uint wmi, uint maxp)
{
	if (start + sieve_size < (uint64)maxp * maxp) {
		maxp = isqrt(start + sieve_size) + 1;
	}

	const uint sieve_byte = sieve_size / WHEEL30 + WheelInit30[sieve_size % WHEEL30].WheelIndex;
	WheelPrime* pwheel = MediumSieve + wmi;

#if W30
	uint minwp = MIN(maxp, sieve_byte * 6 / 10);
	minwp = (minwp / WHEEL << SI_BIT); //  + WHEEL_INIT[minwp % WHEEL].PrimeIndex;
	pwheel = sieveMediumW30(bitarray, sieve_byte, minwp, pwheel);
#endif

	const uint maxwp = (maxp / WHEEL << SI_BIT) + WHEEL_INIT[maxp % WHEEL].PrimeIndex;
	bitarray += sieve_byte;

	while (pwheel[2].Wp < maxwp) {
		pwheel += sieveMedium3(bitarray, sieve_byte, pwheel);
	}
	if (pwheel[0].Wp < maxwp)
		sieveMedium1(bitarray, sieve_byte, pwheel ++);
	if (pwheel[0].Wp < maxwp)
		sieveMedium1(bitarray, sieve_byte, pwheel);
}

//sieve big prime from bucket
static void sieveBig1(uchar bitarray[], const uint sieve_size, const WheelPrime* cur_wheel)
{
	uint offset = cur_wheel->Si, wp = cur_wheel->Wp;
	const WheelElement* wh = Wheel210[wp % (1 << WP_BIT)];
	const WheelElement* wheel = wh + offset % (1 << SI_BIT);
	bitarray[offset >>= SI_BIT] |= wheel->UnsetBit;
	offset += wheel->Correct + wheel->NextMultiple * (wp >> WP_BIT);

	if (offset < sieve_size) {
		wheel = wh + wheel->WheelIndex;
		bitarray[offset] |= wheel->UnsetBit;
		offset += wheel->Correct + wheel->NextMultiple * (wp >> WP_BIT);
	}
	pushBucket(offset, wp, wheel->WheelIndex);
}

//sieve 2 big prime from bucket, 15% improvement
static void sieveBig2(uchar bitarray[], const uint sieve_size, const WheelPrime* cur_wheel)
{
	uint offset1 = cur_wheel[0].Si, wp1 = cur_wheel[0].Wp;
	uint offset2 = cur_wheel[1].Si, wp2 = cur_wheel[1].Wp;
	const WheelElement* wh1 = Wheel210[wp1 % (1 << WP_BIT)];
	const WheelElement* wh2 = Wheel210[wp2 % (1 << WP_BIT)];

	const WheelElement* wheel1 = wh1 + offset1 % (1 << SI_BIT);
	const WheelElement* wheel2 = wh2 + offset2 % (1 << SI_BIT);
	bitarray[offset1 >>= SI_BIT] |= wheel1->UnsetBit;
	offset1 += wheel1->Correct + wheel1->NextMultiple * (wp1 >> WP_BIT);

	bitarray[offset2 >>= SI_BIT] |= wheel2->UnsetBit;
	offset2 += wheel2->Correct + wheel2->NextMultiple * (wp2 >> WP_BIT);
#if ERAT_BIG > 2
	if (offset1 < sieve_size) {
		wheel1 = wh1 + wheel1->WheelIndex;
		bitarray[offset1] |= wheel1->UnsetBit;
		offset1 += wheel1->Correct + wheel1->NextMultiple * (wp1 >> WP_BIT);
	}
	if (offset2 < sieve_size) {
		wheel2 = wh2 + wheel2->WheelIndex;
		bitarray[offset2] |= wheel2->UnsetBit;
		offset2 += wheel2->Correct + wheel2->NextMultiple * (wp2 >> WP_BIT);
	}
#endif
	pushBucket(offset2, wp2, wheel2->WheelIndex);
	pushBucket(offset1, wp1, wheel1->WheelIndex);
}

//This implementation uses a sieve array with WHEEL210 numbers per byte and
//a modulo wheel that skips multiples of 2, 3, 5 and 7.
static void eratSieveBig(uchar bitarray[], const uint sieve_size)
{
	uint loops = (uint)(size_t)Bucket[0].Wheel % MEM_WHEEL / sizeof(WheelPrime);
	if (loops % 2) {
		sieveBig1(bitarray, sieve_size, Bucket[0].Wheel - 1), loops --;
		//*(Bucket[0].Wheel) = *(Bucket[0].Wheel - 1), loops ++;
	} else if (loops == 0) {
		loops = WHEEL_SIZE;
	}

	for (Stock* phead = Bucket[0].Head; phead != NULL; loops = WHEEL_SIZE) {
		WheelPrime* cur_wheel = phead->RWheel;

		while (loops) {
			sieveBig2(bitarray, sieve_size, cur_wheel);
			loops -= 2, cur_wheel += 2;
		}

		Stock* pnext = phead->Next; phead->Next = StockHead; StockHead = phead;
		phead = pnext;
		BucketInfo.CurStock ++;
	}

	BucketInfo.MaxBucket --;
	memmove(Bucket, Bucket + 1, BucketInfo.LoopSize * sizeof(Bucket[0]));
}

static int segmentProcessed(uchar bitarray[], const uint64 start, const uint bytes, Cmd* cmd)
{
	uint64& lastqw = *(uint64*)(bitarray + bytes);
	const uint64 tmp = lastqw;
	lastqw = ~0;
	int primes = countBit0sArray((uint64*)bitarray, bytes * sizeof(uint64));
	lastqw = tmp;
	return primes;
}

//sieve prime multiples in [start, start + sieve_size)
static int segmentedSieve(uchar bitarray[], const uint64 start, const uint sieve_size)
{
	const uint sqrtp = isqrt(start + sieve_size);
	const uint medium = MIN(sqrtp, Threshold.Medium) + 1;
	const bool bmsieve = sqrtp >= Threshold.L1Maxp;

	preSieve(bitarray, start, sieve_size);
	//copy last L1 seg to the first seg
	const uint copy_size = Threshold.L1Maxp, copy_index = Config.SieveSize;
	for (uint i = 0; i < copy_size; i += sizeof(uint64)) {
		uint64& c = *(uint64*)(bitarray + i + copy_index);
		*(uint64*)(bitarray + i) |= c, c = 0;
	}

	for (uint offset = 0, l2size = Threshold.L2Size * WHEEL30; offset < sieve_size; offset += l2size) {
		if (l2size + offset > sieve_size)
			l2size = sieve_size - offset;

		uchar* buffer = bitarray + offset / WHEEL30;
		eratSieveSmall(buffer, start + offset, l2size);
		if (bmsieve) {
			eratSieveMedium(buffer, start + offset, l2size, Threshold.L1Index, Threshold.L2Maxp);
		}
	}

	if (medium > Threshold.L2Maxp) {
		eratSieveMedium(bitarray, start, sieve_size, Threshold.L2Index, medium);
	}

	if (start >= Threshold.BucketStart) {
		eratSieveBig(bitarray, Config.SieveSize);
	}

	return sieve_size / WHEEL30 + (7 + WheelInit30[sieve_size % WHEEL30].WheelIndex) / 8;
}

static int segmentedSieve2(uchar bitarray[], uint start, uint sieve_size)
{
	const uint sqrtp = isqrt(start + sieve_size);
	const uint align30 = (uint)(start % WHEEL30);
	const uint l1_maxp = Threshold.L1Maxp;
	assert(align30 % WHEEL210 == 0);
	const uint bytes = sieve_size / WHEEL30 + (7 + WheelInit30[sieve_size % WHEEL30].WheelIndex) / 8;

	preSieve(bitarray, start, sieve_size);
	*(uint64*)(bitarray + bytes) = ~0;

	for (uint offset = 0, l1size = Threshold.L1Size * WHEEL30; offset < sieve_size; offset += l1size) {
		if (l1size + offset > sieve_size)
			l1size = sieve_size - offset;
		eratSieveL0(bitarray + offset / WHEEL30, start + offset, l1size, l1_maxp);
	}

//	assert(l1_maxp > start / l1_maxp);
	const uchar* pend = bitarray + sieve_size / WHEEL30;
	for (uint j = Threshold.L1Index, p = l1_maxp; p <= sqrtp; p = Prime[++j]) {
		uint offset = p - start % p;
		if (p * p > start)
			offset = p * p - start;
		const WheelFirst& wf = WheelFirst30[offset % WHEEL30][WheelInit30[p % WHEEL30].PrimeIndex];
		const uint multiples = (WHEEL_SKIP >> (wf.NextMultiple * 3)) | (WHEEL_SKIP << (24 - wf.NextMultiple * 3));
		sieveSmall1(bitarray, pend, p, offset + wf.Correct * p, multiples);
	}

	return bytes;
}

static void setThresholdL1()
{
	if (Threshold.L1Maxp > Prime[PRIME_SIZE - 100])
		Threshold.L1Maxp = Prime[PRIME_SIZE - 100];
	for (uint p = Prime[1], j = 1; ; p = Prime[++j]) {
		if (p >= Threshold.L1Maxp && p % WHEEL210 == 1) {
			Threshold.L1Index = j;
			Threshold.L1Maxp = p;
			return;
		}
	}
	//assert(false);
}

void setCpuCache(int level, uint cache)
{
//	cache = 1 << ilog(cache, 2);
	if (level == 1 && cache >= 16 && cache <= (Threshold.L2Size >> 10)) {
		Threshold.L1Size = cache << 10;
		Threshold.L1Maxp = Threshold.L1Size / Threshold.L1Segs;
		Threshold.L2Size = Threshold.L2Size / Threshold.L1Size * Threshold.L1Size;
		Threshold.L2Maxp = Threshold.L2Size / Threshold.L2Segs;
		setThresholdL1();
	} else if (level == 2 && cache >= (Threshold.L1Size >> 10) && cache <= MAX_SEGMENT) {
		Threshold.L2Size = (cache << 10) / Threshold.L1Size * Threshold.L1Size;
		Threshold.L2Maxp = Threshold.L2Size / Threshold.L2Segs;
	}
}

void setLevelSegs(uint level, uint segs)
{
	if (segs > 0 && segs <= 6) {
		if (level == 1) {
			Threshold.L1Segs = segs;
			Threshold.L1Maxp = Threshold.L1Size / segs;
			setThresholdL1();
		} else if(level == 2) {
			Threshold.L2Segs = segs;
			Threshold.L2Maxp = Threshold.L2Size /segs;
		} else if (level == 3 && segs > 1) {
			Threshold.Msegs = segs;
		}
	}
}

//sieve_size : 32k - 4M
int setSieveSize(uint sieve_size)
{
	const uint l1_size = Threshold.L1Size >> 10, l2_Size = Threshold.L2Size >> 10;
	if (sieve_size <= MAX_SEGMENT && sieve_size > l2_Size) {
		sieve_size = (sieve_size / l2_Size) * l2_Size << 10;//  1 << (ilog(sieve_size, 2) + 10);
	} else if (sieve_size <= l2_Size && sieve_size >= l1_size) {
		sieve_size = (sieve_size / l1_size) * l1_size << 10;
	} else if (sieve_size <= (MAX_SEGMENT >> 10) && sieve_size > 0) {
		sieve_size = sieve_size << 20;
	} else {
		sieve_size = L2_DCACHE_SIZE << 10;
	}

	return Config.SieveSize = sieve_size;
}

static int checkSmall(const uint64 start, const uint64 end, Cmd* cmd)
{
	int primes = 0;
	for (int i = 0, p = 2; p < 7; p = Prime[++i]) {
		if (start <= p && p <= end) {
			primes ++;
			if (cmd && cmd->Oper == PCALL_BACK) {
				(*(sieve_call)cmd->Data)(primes, p);
				cmd->Primes += 1;
			} else if (cmd && cmd->Oper == SAVE_PRIME) {
				((uint64*)cmd->Data)[primes - 1] = p;
				cmd->Primes += 1;
			}
		}
	}

	return primes;
}

static int64 pi(uchar* bitarray, uint64 start, uint64 end, Cmd* cmd)
{
	const int64 ts = getTime();
	uint64 primes = checkSmall(start, end, cmd);

	uint align210 = (uint)(start % WHEEL210);
	start -= align210;

	if (++end == 0) end --; //watch overflow if end = 2^64-1
	//pi(n) ~= n/log(n), complexty ~= n*log(log(n)), replaced by n*log(n) more accurate
#ifndef _DEBUG
	double lge = end * log((double)end), lgs = start * log(start + 10.0);
	double pie = PI(end, 1.2), pis = PI(start, 1.2);
#endif

	for (uint si = 0, sieve_size = Config.SieveSize * WHEEL30; start < end; start += sieve_size) {
		if (sieve_size > end - start) {
			sieve_size = (uint)(end - start);
		}

		const uint bytes = segmentedSieve(bitarray, start, sieve_size);
		if (align210 > 0) {
			memset(bitarray, ~0, align210 / WHEEL30);
			bitarray[align210 / WHEEL30] |= (1 << WheelInit30[align210 % WHEEL30].WheelIndex) - 1;
			align210 = 0;
		}

		primes += segmentProcessed(bitarray, start, bytes, cmd);
#ifndef _DEBUG
		if ((si ++ & Config.Progress) == 16) {
			const double cur = start + sieve_size;
			double lgc = cur * log(cur), pic = PI(cur, 1.2);
			double tratio = (lgc - lgs) / (lge - lgs) * 100;
			double pratio = (pic - pis) / (pie - pis);
			double timeuse = (getTime() - ts) / (10 * tratio);
			const uint64 picount = (int64)((int64)primes / pratio);
			if (timeuse < 3600)
				printf(">> %.2f%%, time ~= %.2f %s, primes ~= %llu\r", tratio, timeuse, "sec", picount);
			else
				printf(">> %.2f%%, time ~= %.3f %s, primes ~= %llu\r", tratio, timeuse / 3600, "hour", picount);
			fflush(stdout);
		}
#endif
	}

	return primes;
}

static void convertSci(uint64 n, char buff[32])
{
	const int logn = ilog(n, 10);
	const uint64 pown = ipow(10, logn);
	if (n % pown == 0)
		sprintf(buff, "%de%d", (int)(n / pown), logn);
	else if (n % (pown / 10) == 0)
		sprintf(buff, "%de%d", (int)(n / (pown / 10)), logn - 1);
	else if ((n & (n - 1)) == 0)
		sprintf(buff, "2^%d", ilog(n, 2));
	else if (n > 1000000000l) {
		uint64 r = n - (int)(n / pown) * pown;
		const int logr = ilog(r, 10);
		const uint64 powr = ipow(10, logr);
		if (r % powr == 0)
			sprintf(buff, "%de%d+%de%d", (int)(n / pown), logn, (int)(r / powr), logr);
	}
}

static void printResult(const uint64 start, const uint64 end, uint64 primes)
{
	char buff[128] = {0};
	char begbuff[40] = {0}, endbuff[40] = {0}, rangebuff[40] = {0};
	if (end > 10000) {
		convertSci(start, begbuff);
		convertSci(end, endbuff);
	}
	if (end + 1 == 0)
		sprintf(endbuff, "2^64");

	const uint64 range = end - start;
	if (range > 10000)
		convertSci(range, rangebuff);

	if (start == 0) {
		if (endbuff[0] == 0)
			sprintf(buff, "%llu", end);
		else
			sprintf(buff, "%s", endbuff);
	} else if (begbuff[0] && endbuff[0]) {
		if (rangebuff[0] && range <= start / 10)
			sprintf(buff, "%s,+%s", begbuff, rangebuff);
		else
			sprintf(buff, "%s,%s", begbuff, endbuff);
	} else if (rangebuff[0]) {
		if (begbuff[0])
			sprintf(buff, "%s,+%s", begbuff, rangebuff);
		else if (endbuff[0])
			sprintf(buff, "%s-%s,%s", endbuff, rangebuff, endbuff);
		else
			sprintf(buff, "%llu,+%s", start, rangebuff);
	} else {
		if (begbuff[0] == 0)
			sprintf(begbuff, "%llu", start);
		if (range < start / 10 && endbuff[0] == 0)
			sprintf(endbuff, "+%llu", range);
		else if (endbuff[0] == 0)
			sprintf(endbuff, "%llu", end);
		sprintf(buff, "%s,%s", begbuff, endbuff);
	}
	printf("\rpi(%s) = %llu", buff, primes);
}

static uint setBucketStart(const uint64 start, const uint sqrtp, const uint sieve_size)
{
	uint medium = sieve_size / Threshold.Msegs + 1;

	uint64 offset = (uint64)medium * medium;
	if (sqrtp > medium && offset > start) {
		offset += sieve_size - offset % sieve_size + start % sieve_size;
		while (offset % WHEEL210 != start % WHEEL210) offset += sieve_size;
		medium = isqrt(offset - offset % WHEEL210 + sieve_size) + 1;
	} else {
		offset = start;
	}

	medium = MIN(sqrtp, medium);
	medium += WHEEL210 - medium % WHEEL210;

	Threshold.BucketStart = offset - offset % WHEEL210;
	Threshold.Medium = medium;
	return medium;
}

uint64 doSieve(const uint64 start, const uint64 end, Cmd* cmd = NULL)
{
	const int64 ts = getTime();
	const uint sqrtp = isqrt(end);
	uint sieve_size = Config.SieveSize;

	const bool bwheel = sqrtp > (sieve_size * WHEEL30) / Threshold.Msegs;
	if (bwheel) {
		assert((Threshold.L1Size & (Threshold.L1Size - 1)) == 0);
		if (sieve_size < Threshold.L2Size || sieve_size >= Threshold.Msegs << 21)
			sieve_size = setSieveSize(SIEVE_SIZE);
		else if ((sieve_size & (sieve_size - 1)) != 0)
			sieve_size = setSieveSize(1 << ilog(sieve_size >> 10, 2));
	}

	//init medium sieve
	const int max_cache = MAX(sieve_size, L2_DCACHE_SIZE * 1024) + Threshold.L1Size;
	uchar* bitarray = (uchar*) malloc(max_cache + 1024);
	sieve_size *= WHEEL30;

	const uint medium = setBucketStart(start, sqrtp, sieve_size); //TODO
	const uint pix = PI(medium, 1.2);
	if (MediumSieve != NULL && (MediumSieve[0].Wp < pix || MediumSieve[0].Wp > pix * 3)) {
		free(MediumSieve); MediumSieve = NULL;
	}
	if (MediumSieve == NULL) {
		SmallSieve = MediumSieve = (WheelPrime*) malloc(sizeof(WheelPrime) * (pix + 100));
		MediumSieve[0].Wp = pix;
	}

	const bool bmsieve = sqrtp >= Threshold.L1Maxp;
	if (bmsieve) {
		setWheelMedium(bitarray, sieve_size, medium + 255, start - start % WHEEL210); //TODO
	}

	//init bucket sieve
	uint64& buckets = Threshold.BucketStart;
	if (sqrtp > medium) {
		setBucketInfo(sieve_size, sqrtp, end - buckets);
		setWheelBig(bitarray, medium, sqrtp, buckets, end - buckets);
		if (BucketInfo.CurStock < BucketInfo.LoopSize)
			allocWheelBlock(1);
	} else {
		buckets = end + 1;
	}

	//init small sieve
	setWheelSmall(start);
	memset(bitarray + sieve_size / WHEEL30, 0, Threshold.L1Maxp + 64);

	const int64 ti = getTime();
	const int64 primes = pi(bitarray, start, end, cmd);

	free(bitarray);
	if (BucketInfo.CurStock > MAX_STOCK / 5) {
		for (uint i = 0; i < BucketInfo.PoolSize; i ++) free(WheelPool[i]);
		memset(&BucketInfo, 0, sizeof(BucketInfo));
		memset(Bucket, 0, sizeof(Bucket));
		StockHead = NULL;
	}

	if (Config.Flag & PRINT_RET) {
		const int64 ta = getTime();
		printResult(start, end, primes);
		if (Config.Flag & PRINT_TIME)
			printf(" (%.2f + %.2f = %.2f sec, %d kb L%d%d%d)",
					(ti - ts) / 1000.0, (ta - ti) / 1000.0, (ta - ts) / 1000.0, Config.SieveSize >> 10, Threshold.L1Segs, Threshold.L2Segs, Threshold.Msegs);
		putchar('\n');
	}

#ifndef B_R
	assert(BucketInfo.StockSize == BucketInfo.CurStock);
#endif

	return primes;
}

//the simple sieve of Eratosthenes implementated by bit packing
static int eratoSimple()
{
	int primes = 1;
	const uint maxp = (1 << 16) + 32;
	uchar bitarray[(maxp >> 4) + 10] = {0};

	for (uint p = 3; p < maxp; p += 2) {
		if (0 == (bitarray[p >> 4] & (1 << (p / 2 & 7)))) {
			Prime[primes ++] = p;
			for (uint j = p * (p / 2) + p / 2; j <= maxp / 2; j += p)
				bitarray[j >> 3] |= 1 << (j & 7);
		}
	}

	return primes;
}

//The first presieved template, cross off the first 8th prime multiples
static void initPreSieved()
{
	for (int i = 3; PRIME_PRODUCT % Prime[i] == 0; i ++) {
		int p = Prime[i];
		for (int offset = p; offset < sizeof(PreSieved) * WHEEL30; offset += p * 2) {
			PreSieved[offset / WHEEL30] |= WheelInit30[offset % WHEEL30].UnsetBit;
		}
	}
}

static void initBitTable()
{
	int i = 0;
	const int nbitsize = sizeof(WordNumBit1) / sizeof (WordNumBit1[0]);
	for (i = 1; i < nbitsize; i ++)
		WordNumBit1[i] = WordNumBit1[i >> 1] + (i & 1);

	uchar pattern[] = {1, 7, 11, 13, 17, 19, 23, 29};
	for (i = 0; i < sizeof(Pattern30) / sizeof (Pattern30[0]); i ++)
		Pattern30[i] = pattern[i % 8] + WHEEL30 * (i / 8);

#if BIT_SCANF == 0
	const int nbitsize2 = sizeof(Lsb) / sizeof(Lsb[0]);
	for (i = 0; i < nbitsize2; i += 2)
		Lsb[i + 0] = Lsb[i >> 1] + 1;

	for (i = 0; i < nbitsize2; i += 2) {
		Lsb[i + 0] = Pattern30[Lsb[i]];
		Lsb[i + 1] = Pattern30[0];
	}
#endif
}

static void initWheel30()
{
	//how to calculate the magic number ?
	const uint nextMultiple[8] =
	{
		0x74561230, 0x67325401,
		0x51734062, 0x42170653,
		0x35607124, 0x26043715,
		0x10452376, 0x03216547
	};

	for (int j = 0, wj = 0; j < WHEEL30; j ++) {
		WheelInit30[j].UnsetBit = 0;
		WheelInit30[j].WheelIndex = wj;
		WheelInit30[j].PrimeIndex = wj;
		if (j == Pattern30[wj]) {
			WheelInit30[j].UnsetBit = 1 << (wj ++);
		}
	}

	for (int i = 0; i < WHEEL30; i ++) {
		for (int pi = 0; pi < 8; pi ++) {
			int multiples = 0, offset = i;
			if (i % 2 == 0) {
				multiples = 1;
				offset += Pattern30[pi];
			}
			while (WheelInit30[offset % WHEEL30].UnsetBit == 0) {
				offset += Pattern30[pi] * 2;
				multiples += 2;
			}
			int wi = WheelInit30[offset % WHEEL30].WheelIndex;
			WheelElement& wf = WheelFirst30[i][pi];
			wf.NextMultiple = (nextMultiple[wi] >> (pi * 4)) & 15;
			wf.WheelIndex = wi;
			wf.Correct = multiples;
			wf.UnsetBit = 1 << wi;
		}
	}

	for (int wi = 0; wi < 8; wi ++) {
		for (int pi = 0; pi < 8; pi ++) {
			int multiples = 2;
			int next = Pattern30[wi] + Pattern30[pi] * 2;
			while (WheelInit30[next % WHEEL30].UnsetBit == 0) {
				next += Pattern30[pi] * 2;
				multiples += 2;
			}

			WheelElement& we30 = Wheel30[pi][wi];
			we30.NextMultiple = multiples * (WHEEL30 / WHEEL30);
			we30.WheelIndex = WheelInit30[next % WHEEL30].WheelIndex;
			we30.Correct = next / WHEEL30 - Pattern30[wi] / WHEEL30;
			we30.UnsetBit = WheelInit30[Pattern30[wi]].UnsetBit;
		}
	}
}

static void initWheel210()
{
	int wi = 0, i = 0;
	const int psize = 48, wsize = 48;
	int wpattern[WHEEL210] = {0};

	for (i = 0; i < WHEEL210; i ++) {
		const uchar mask = WheelInit30[i % WHEEL30].UnsetBit;
		WheelInit210[i].UnsetBit = mask;
		WheelInit210[i].WheelIndex = wi;
		WheelInit210[i].PrimeIndex = wi;
		if (mask && i % (WHEEL210 / WHEEL30)) {
			Pattern210[wi] = wpattern[wi] = i;
			wi ++;
		}
	}

	for (wi = 0; wi < wsize; wi ++) {
		for (int pi = 0; pi < psize; pi ++) {
			int multiples = 2;
			int next = wpattern[wi] + wpattern[pi] * 2;

			while (WheelInit210[next % WHEEL210].WheelIndex == WheelInit210[(next + 1) % WHEEL210].WheelIndex) {
				next += wpattern[pi] * 2;
				multiples += 2;
			}

			WheelElement& we210 = Wheel210[pi][wi];
			we210.Correct = next / WHEEL30 - wpattern[wi] / WHEEL30;
			we210.UnsetBit = WheelInit210[wpattern[wi]].UnsetBit;
			we210.WheelIndex = WheelInit210[next % WHEEL210].WheelIndex;
			we210.NextMultiple = multiples * (WHEEL210 / WHEEL30);
		}
	}

	for (i = 0; i < WHEEL210; i ++) {
		for (int pi = 0; pi < psize; pi ++) {
			int multiples = 0, next = i;
			if (i % 2 == 0) {
				multiples = 1;
				next += wpattern[pi];
			}

			while (WheelInit210[next % WHEEL210].WheelIndex == WheelInit210[(next + 1) % WHEEL210].WheelIndex) {
				next += wpattern[pi] * 2;
				multiples += 2;
			}

			WheelFirst210[i][pi].WheelIndex = WheelInit210[next % WHEEL210].WheelIndex;
			WheelFirst210[i][pi].Correct = multiples;
		}
	}
}

static void fixRangeTest(uint64 lowerBound, const int64 range, uint64 Ret)
{
	const int llog10 = ilog(lowerBound, 10), rlog10 = ilog(range, 10);
	const uint64 maxrange = ipow(2, 32);
	uint64 primes = 0, upperBound = lowerBound + range;
	if (upperBound + 1 == 0)
		printf("Sieving Pi[2^64-10^%d, 2^64-1] with range 10^%d\n", rlog10, ilog(maxrange, 10));
	else
		printf("Sieving pi[10^%d, 10^%d+10^%d] with range 10^%d randomly\n", llog10, llog10, rlog10, ilog(maxrange, 10));

	while (lowerBound < upperBound) {
		uint64 rd = rand() * rand();
		uint64 end = lowerBound + (rd * rand() * rd) % maxrange;
		if (end > upperBound || end < lowerBound)
			end = upperBound;

//		setCpuCache(1, 32 * (rand() % 2 + 1));
		setCpuCache(2, 256 * (rand() % 2 + 1));
		setSieveSize(rand() % MAX_SEGMENT + 64);
		setLevelSegs(1, rand() % 5 + 1), setLevelSegs(2, rand() % 5 + 2), setLevelSegs(3, rand() % 5 + 2);

		primes += doSieve(lowerBound, end - 1);
		if (lowerBound % 4 == 0)
		printf("chunk: %.2f%%\r", 100 - (int64)(upperBound - lowerBound) * 100.0 / range), fflush(stdout);
		lowerBound = end;
	}
	if (lowerBound + 1 == 0)
		printf("Pi[2^64-10^%d, 2^64-1] = %llu\n", rlog10, primes);
	else
		printf("Pi[10^%d, 10^%d+10^%d] = %llu\n", llog10, llog10, rlog10, primes);

	assert(primes == Ret);
}

static void startBenchmark(int flag)
{
	const uint primeCounts[] =
	{
		4,         // pi(10^1)
		25,        // pi(10^2)
		168,       // pi(10^3)
		1229,      // pi(10^4)
		9592,      // pi(10^5)
		78498,     // pi(10^6)
		664579,    // pi(10^7)
		5761455,   // pi(10^8)
		50847534,  // pi(10^9)
		455052511, // pi(10^10)
		4118054813u,//pi(10^11)
		203280221, // pi(2^32)
		155428406, // pi(10^12, 10^12 + 2^32)
		143482916, // pi(10^13, 10^13 + 2^32)
		133235063, // pi(10^14, 10^14 + 2^32)
		124350420, // pi(10^15, 10^15 + 2^32)
		116578809, // pi(10^16, 10^16 + 2^32)
		109726486, // pi(10^17, 10^17 + 2^32)
		103626726, // pi(10^18, 10^18 + 2^32)
		98169972,  // pi(10^19, 10^19 + 2^32)
		0,
	};

	int64 ts = getTime();
	Config.Flag &= ~PRINT_RET;
	uint64 primes = 0;
	Config.Progress = 0;
	for (int i = 1; i <= 10; i ++) {
		primes = doSieve(0, ipow(10, i));
		printf("pi(10^%2d) = %llu\n", i, primes);
	}

	srand(time(0));
	for (int j = 12; primeCounts[j]; j ++) {
		uint64 start = ipow(10, j), end = start + ipow(2, 30);
		setLevelSegs(1, rand() % 6 + 2), setLevelSegs(2, rand() % 6 + 2); setLevelSegs(3, rand() % 6 + 2);
		primes = doSieve(start, end);
		if (primes == primeCounts[j])
			printf("pi(10^%d, 10^%d+2^32) = %llu                 \n", j, j, primes);
	}

	Config.Progress = 0;
	printf("Time elapsed %.f sec\n\n", (getTime() - ts) / 1000.0);
	puts("All Big tests passed SUCCESSFULLY!\nStart Rand Test");

	const uint64 pow11 = ipow(10, 11);
	const uint64 pow12 = pow11 * 10, pow9 = pow11 / 100;

	const uint64 rangeData[][3] =
	{
		{ipow(01, 11), pow11, 4118054813ul},
		{-1 - pow11,   pow11, 2254197466ul},
		{ipow(10, 17), pow11, 2554712095ul},
		{ipow(10, 12), pow11, 3612791400ul},
		{ipow(10, 13), pow11, 3340141707ul},
		{ipow(01, 12), pow12, pow9 * 37 + 607912018},
		{ipow(10, 14), pow12, pow9 * 31 + 16203073},
		{ipow(10, 15), pow12, pow9 * 28 + 952450479},
		{ipow(10, 16), pow12, pow9 * 27 + 143405794},
		{ipow(10, 18), pow12, pow9 * 24 + 127637783},
		{ipow(10, 19), pow12, pow9 * 22 + 857444126},
		{ipow(10, 19), pow12, pow9 * 22 + 857444126},
		{-1 - pow12,   pow12, pow9 * 22 + 542106206},
		{-1 - pow12*10, pow12*10, pow9 * 225 + 420940155},
	};

	for (int k = 0; k < sizeof(rangeData) / sizeof(rangeData[0]); k ++) {
		fixRangeTest(rangeData[k][0], rangeData[k][1], rangeData[k][2]);
	}

	Config.Flag |= PRINT_RET;
}

static void printInfo()
{
	const char* sepator =
		"-------------------------------------------------------------------------------------------";
	puts(sepator);
	puts("Fast implementation of the segmented sieve of Eratosthenes (2^64 - 1)\n"
	"Copyright (C) by 2010-2018 Huang YuanBing bailuzhou@163.com\n"
	"Code: https://github.com/ktprime/ktprime/blob/master/PrimeNumber.cpp\n"
	"Compile: g++ -DSIEVE_SIZE=2048 -DW30 -march=native -funroll-loops -O3 -s -pipe PrimeNumber.cpp -o prime\n");

	char buff[500] = {0};
	char* info = buff;
#ifdef __clang_version__
	info += sprintf(info, "Compiled by %s", __clang_version__);
#elif _MSC_VER
	info += sprintf(info, "Compiled by vc ++ %d", _MSC_VER);
#elif __GNUC__
	info += sprintf(info, "Compiled by gcc %s", __VERSION__);
#endif

#if X86_64
	info += sprintf(info, " x86-64");
#endif

	info += sprintf(info, " %s %s\n", __TIME__, __DATE__);
	info += sprintf(info, "[MARCO] : MEM_WHEEL = %dM, WHEEL_SIZE = %dk, SIEVE_SIZE = %dk, WHEEL = %d\n",
			MEM_BLOCK * WHEEL_SIZE >> 17 , WHEEL_SIZE >> 7, SIEVE_SIZE, WHEEL);
	info += sprintf(info, "[CACHE] : L1Size = %u, L2Size = %u, SieveSize = %u, Bucket = %u\n",
			Threshold.L1Size >> 10, Threshold.L2Size >> 10, Config.SieveSize >> 10, BucketInfo.LoopSize);
	info += sprintf(info, "[ARGS ] : L1Segs/L2Segs/Mseg = (%u,%u,%u)\n",
			Threshold.L1Segs, Threshold.L2Segs, Threshold.Msegs);
	info += sprintf(info, "[ARGS ] : L1Maxp/L2Maxp/Medium/Large/SieveSize = (%u,%u,%u,%u,%u)",
			Threshold.L1Maxp, Threshold.L2Maxp, Threshold.Medium, isqrt(Threshold.BucketStart + 1), Config.SieveSize);
	puts(buff);
	puts(sepator);
}

//get the first digit number index
static int parseCmd(const char params[][60])
{
	int cmdi = -1;
	int cdata = 0;

	for (int i = 0; params[i][0]; i ++) {
		unsigned char c = params[i][0], n = params[i][1];
		if (c == '-')
			c = params[i][1];
		if (c >= 'a' && c <= 'z')
			c += 'A' - 'a';
		if (isdigit(c) || c == 'E') {
			if (cmdi < 0)
				cmdi = i;
			continue;
		}
		if (isdigit(n))
			cdata = atoi(params[i] + 1);
		else if (isdigit(params[i][2]))
			cdata = atoi(params[i] + 2);

		switch (c)
		{
			case 'H' :
				puts(Benchmark);
				if (n)
				puts(Help);
				break;
			case 'S':
				setSieveSize(cdata);
				break;
			case 'C':
				setCpuCache(cdata % 10, cdata / 10);
				break;
			case 'L':
				setLevelSegs(cdata % 10, cdata / 10);
				break;
			case 'M':
				Config.Progress = (1 << cdata) - 1;
				break;
			case 'D':
				Config.Flag ^= 1 << (n - 'A');
				break;
			case 'I':
				printInfo();
				break;
			default:
				cmdi = i;
				break;
		}
	}

	return cmdi;
}

//split ccmd string to params array
static int splitCmd(const char* ccmd, char cmdparams[][60])
{
	int nwords = 0;

	for (int i = 0; i < 64; i ++) {
		while (isspace(*ccmd) || ',' == *ccmd) {
			ccmd ++;
		}
		char* pc = cmdparams[i];
		if (*ccmd == 0 || *ccmd == ';') {
			break;
		}
		unsigned char c = *ccmd;
		bool isvalid = false;
		while (isalnum(c) || c == '^' || c == '/' || c == '+' || c == '-' || c == '*') {
			*pc ++ = c;
			c = * ++ccmd;
			isvalid = true;
		}
		if (isvalid)
			nwords ++;
		else
			ccmd ++;
	}

	return nwords;
}

static bool executeCmd(const char* cmd)
{
	initPrime(SIEVE_SIZE);
	while (cmd) {

		char params[14][60] = {0};

		char* pcmd = (char*) strchr(cmd, ';');
		if (splitCmd(cmd, params) <= 0)
			return false;

		int cmdi = parseCmd(params);
		if (cmdi == -1) {
			return true;
		}

		unsigned char cmdc = toupper(params[cmdi][0]);
		uint64 start = atoint64(params[cmdi]);
		uint64 end = atoint64(params[cmdi + 1]);
		if (!isdigit(cmdc) && cmdc != 'E') {
			start = end;
			end = atoint64(params[cmdi + 2]);
		}

		if (cmdc == 'B') {
			if (isdigit(params[cmdi + 1][0])) {
				int powi = atoi(params[cmdi + 1]);
				uint64 range = powi > 12 ? ipow(2, powi) : ipow(10, powi);
				for (int j = 11; j < 20 && powi > 0; j ++) {
					uint64 start2 = ipow(10, j);
					doSieve(start2, start2 + range);
				}
			}
			if (isdigit(params[cmdi + 2][0])) {
				int powi = atoi(params[cmdi + 2]);
				uint64 range = powi > 12 ? ipow(2, powi) : ipow(10, powi);
				for (int i = 32; i < 64 && powi > 0; i ++) {
					uint64 start2 = ipow(2, i);
					doSieve(start2, start2 + range);
				}
			} else
				startBenchmark(0);
		} else if (cmdi >= 0) {
			if (end == 0)
				end = start, start = 0;
			else if (end < start) {
				end += start;
				if (end < start)
					end = -1;
			}
			doSieve(start, end);
		}

		if (pcmd) {
			cmd = pcmd + 1;
		} else {
			break;
		}
	}

	return true;
}

void initPrime(int sieve_size)
{
	if (Prime[2] == 0) {
		eratoSimple();
		initBitTable();
		initWheel30();
		initWheel210();
		initPreSieved();
		setThresholdL1();
		setSieveSize(sieve_size);
	}
}

int main(int argc, char* argv[])
{
	if (argc == 2)
		executeCmd(argv[1]);

#if 1
	else if (argc > 2)
	{
	srand(time(0));
	initPrime(SIEVE_SIZE);

	Config.Progress = 0;
	double tl[4][8] = {0};

	const uint maxn = 5 * 5;
	for (int n = 0 ; n < maxn; n++)
	{
		int l3 = n % 5 + 2, l2 = n / 5 % 5 + 2, l1 = n / 25 + 2;
		int64 start = atoint64(argv[1]);
		int64 range = atoint64(argv[2]);

		if (maxn == 25 && (start < 1e14) || range < 1e8) {
			int tmp = l1;
			l1 = l3;
			l3 = tmp;
		}
		setLevelSegs(1, l1); setLevelSegs(2, l2); setLevelSegs(3, l3);

		double ts = getTime();
		doSieve(start, start + range);
		double ta = getTime() - ts;
		tl[1][l1] += ta;
		tl[2][l2] += ta;
		tl[3][l3] += ta;
	}
	for (int l = 1; l <= 3; l++)
	{
		printf("L%d ", l);
		for (int i = 2; i <= 6; i++)
			printf("%.2lf, ", tl[l][i] / maxn / 5);
		printf("\n");
	}

	Config.Flag ^= PRINT_RET;
	Config.Flag &= ~PRINT_TIME;
	for (int i = 10; i > 0; i --)
	for (int j = 1; j < 10; j ++) {
			setCpuCache(1, 32 * (rand() % 2 + 1)); setCpuCache(2, 256 * (rand() % 2 + 1));
			uint64 start = ((uint64)(rand() * rand()) << i) + (uint64)rand() * rand() * rand();
			int sieve_size = setSieveSize(rand() % 2000 + 32);
			uint64 range = (rand() % 32) * Config.SieveSize * WHEEL30 + rand();
			uint64 rm1 = doSieve(start, start + range, NULL);
			if (setSieveSize(L2_DCACHE_SIZE << (rand() % 5 + 0)) == sieve_size)
				setSieveSize(sieve_size*2);

			Config.Flag ^= PRINT_RET;
			setLevelSegs(3, rand() % 5 + 2); setLevelSegs(2, rand() % 6 + 2), setLevelSegs(1, rand() % 5 + 1);
			uint64 rm2 = doSieve(start, start + range, NULL);
			if(rm1 != rm2) { printInfo(); printf("%llu %llu\n", start, range); system("pause"); }
			if (j % 100 == 0) printf("\r %2d progress ~= %d\n", i, j);
			Config.Flag ^= PRINT_RET;
	}

	Config.Flag &= ~PRINT_TIME;
	Config.Flag ^= PRINT_RET;

	for (int j = 1; j <= 1000; j ++) {
		char cmd[100] = {0};
#if 1
		setCpuCache(1, 32 * (rand() % 2 + 1)); setCpuCache(2, 256 * (rand() % 2 + 1));
		setLevelSegs(1, rand() % 6 + 1);
		uint64 start = (uint64)(rand() * rand()) * (uint64)(rand() * rand());
		setSieveSize(L2_DCACHE_SIZE * rand() % 8 + L2_DCACHE_SIZE);
		start -= start % 2;
		uint64 range = (rand() % 20 + 2) * Config.SieveSize * WHEEL30;
		uint64 rm1 = doSieve(start, start + range - 1, NULL);

		Config.Flag ^= PRINT_RET;

		setLevelSegs(2, rand() % 8 + 3);
		setSieveSize(L2_DCACHE_SIZE << ((rand() % 5) + 1));
		uint64 rm2 = doSieve(start, start + range, NULL);
		sprintf(cmd, "%lld != %lld [%d] --- %lld %lld s%d L%d1 L%d2 L%d3", rm1, rm2,
				j, start, range, Config.SieveSize >> 10, Threshold.L1Segs, Threshold.L2Segs, Threshold.Msegs);
		if (rm1 != rm2)
			puts(cmd);
#endif
		for (int i = 2; i <= 6; i ++) {
			uint sieve_size = Config.SieveSize * WHEEL30;
			uint64 medium = sieve_size / i + WHEEL210;
			medium -= medium % WHEEL210;
			uint64 start = medium * medium - sieve_size * (rand() % 6 + 1) + rand();
			uint64 end = start + sieve_size * (rand() % 18 + 1) + rand();

			setLevelSegs(3, i); setLevelSegs(2, rand() % 6 + 2), setLevelSegs(1, rand() % 6 + 1);
			const uint64 r1 = doSieve(start, end, NULL);
			char cmd2[100] = {0};
			sprintf(cmd2, "  %llu %llu %llu s%d L%d1 L%d2 L%d3",
					r1, start, end, (sieve_size / WHEEL30) >> 10, Threshold.L1Segs, Threshold.L2Segs, Threshold.Msegs);

			setSieveSize(2048);
			setLevelSegs(3, rand() % 5 + 2); setLevelSegs(2, rand() % 6 + 2), setLevelSegs(1, rand() % 5 + 1);
			const uint64 r2 = doSieve(start, end, NULL);
			if (r1 != r2) {
				puts(cmd);
				puts(cmd2);
				printf("   %llu != %llu, %llu %llu s%d L%d1 L%d2 L%d3\n",
						r1, r2, start, end, Config.SieveSize >> 10, Threshold.L1Segs, Threshold.L2Segs, Threshold.Msegs);
			}
		}
		if (j % 10 == 0) printf("\rprogress ~= %d\n", j), fflush(stdout);
	}}
#endif

	executeCmd("e12 e10 h; e14 e10 s2");
	executeCmd("e10 s2; e12 e9; e16 e9;e18 e9;");

	while (true) {
		char ccmd[257];
		printf(">> ");
#if _MSC_VER > 1600 || GETS_S > 0
		if (!gets_s(ccmd) || !executeCmd(ccmd))
#else
		if (!gets(ccmd) || !executeCmd(ccmd))
#endif
			break;
	}

	return 0;
}

//TODO
//1.multi-threading
//2.improve bucket algorithm for big range ex[1e14, 1e16]
//3.redesign by c++
