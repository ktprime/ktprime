/* Fast segmented sieving of prime number which is based on site http://sweet.ua.pt/tos/software/prime_sieve.html ***/
const char* Benchmark =
"Mingw: g++ 5.1.0 bailuzhou@163\n"
"Windows 10 x64                   i3 350M,i5 3470,i7 7500u\n"
"pi(1e11, 1e10)  =  394050419     5.20    2.91    2.30\n"
"pi(1e12, 1e10)  =  361840208     6.31    3.40    2.70\n"
"pi(1e13, 1e10)  =  334067230     7.56    3.94    3.21\n"
"pi(1e14, 1e10)  =  310208140     9.24    4.64    3.87\n"
"pi(1e15, 1e10)  =  289531946     11.2    5.50    4.64\n"
"pi(1e16, 1e10)  =  271425366     13.0    6.40    5.42\n"
"pi(1e17, 1e10)  =  255481287     15.1    7.45    6.15\n"
"pi(1e18, 1e10)  =  241272176     18.3    8.88    7.35\n"
"pi(1e19, 1e10)  =  228568014     25.2    11.9    9.95\n"
"pi(2^64-1e9,1e9)=  22537866      8.75    4.28    3.62\n"
"pi(1e18, 1e6)   =  24280         0.75    0.46    0.34\n"
"pi(1e18, 1e8)   =  2414886       1.50    0.81    0.70\n"
"pi(1e18, 1e9)   =  24217085      3.80    1.80    1.48\n"
"pi(1e14, 1e12)  =  31016203073   930     466     324\n"
"pi(1e16, 1e12)  =  27143405794   1300    640     462\n"
"pi(1e18, 1e12)  =  24127637783   1600    810     592\n"
"pi(1e19, 1e12)  =  22857444126   1800    900     622\n"
"pi(1e18, 1e13)  =  241274558866          8180    6000\n"
"pi(1e19, 1e13)  =  228575545410          8800    6210\n";

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
	WHEEL_SKIP    = 0x799b799b,
	PATTERN_GAP   = 10,//max prime wheel 210 gap
	PRIME_SIZE    = 6542 + 10,//pi(2^16)
	MAX_SEGMENT   = (1 << 10) * 4,//max sieve_size
};

//performance marco
enum CONFIG
{
	L1_DCACHE_SIZE = 64,
	SI_BIT         = 8, //8 -- 10
	WP_BIT         = 6, //48 < 2^WP_BIT < WHEEL210
#if (L2_DCACHE_SIZE < 256 || L2_DCACHE_SIZE > 4096)
	L2_DCACHE_SIZE =  256,
#endif
#ifndef SIEVE_SIZE
	SIEVE_SIZE     =  2048,
#endif
#ifndef UINT_MAX
	UINT_MAX       = 0-1u,
#endif

	ERAT_SMALL     = 6, //4 - 16
	ERAT_MEDIUM    = 2, //2 - 6
};

enum EBUCKET
{
	UINT_PIMAX = 203280221, //= pi(2^32)
	MAX_BUCKET = 0xffffffff / (256 * (1 << 10) * 3) + 4, //5465 (sqrtp * (8 + 2) / sieve_size + 2);
	WHEEL_SIZE = 1 << 12, //=4096 11: 16k, [10 - 13]
	MEM_WHEEL  = WHEEL_SIZE * 8, //=32768
	MEM_BLOCK  = (1 << 19) / WHEEL_SIZE, //=256 1 << 20:8 MB
	MAX_STOCK  = UINT_PIMAX / WHEEL_SIZE + MAX_BUCKET + MEM_BLOCK, //=55221
	MAX_POOL   = UINT_PIMAX / (MEM_BLOCK * WHEEL_SIZE) + 100, //=487
};

enum EFLAG
{
	SLOW_TEST  = 1 << ('A' - 'A'),
	PRINT_RET  = 1 << ('R' - 'A'),
	PRINT_TIME = 1 << ('T' - 'A'),
	SAVE_DATA  = 1 << ('F' - 'A'),
};

enum ECMD
{
	COUNT_PRIME,
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
	1e16 10^10 s1024 c321";

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
	32 * WHEEL30 << 10, 256 * WHEEL30 << 10,
	(32 << 10) / ERAT_SMALL, 0, ERAT_SMALL,
	(256 << 10) / ERAT_MEDIUM, 0, ERAT_MEDIUM,
	SIEVE_SIZE * (WHEEL30 << 10) / ERAT_BIG, ERAT_BIG,
};

struct _Config
{
	uint Flag;
	uint Progress;
	uint SieveSize;
};

_Config Config =
{
	PRINT_RET | PRINT_TIME,
	(1 << 5) - 1,
	SIEVE_SIZE * (WHEEL30 << 10),
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
	WheelPrime* RWheel; //read only
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
	uint MaxBucket;
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
static WheelPrime* MediumWheel;

//small prime
static uint Prime[PRIME_SIZE];

//presieved buffer with prime <= 19.
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
static WheelElement Wheel210[48][64];

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

//x^n < 2^64
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

#define poSet(s1,o1,b1,s2,o2,b2,s3,o3,b3,s4,o4,b4,s5,o5,b5,s6,o6,b6,s7,o7,b7)\
	while (p /*+ o * s7 + o7*/ <= pend) {\
		p[o * 00 + 00] |= BIT##0,  p[o * s1 + o1] |= BIT##b1; p[o * s2 + o2] |= BIT##b2, p[o * s3 + o3] |= BIT##b3;\
		p[o * s4 + o4] |= BIT##b4, p[o * s5 + o5] |= BIT##b5; p[o * s6 + o6] |= BIT##b6, p[o * s7 + o7] |= BIT##b7;\
		p += step;\
	}

#if 0
	if (p + o *  0 +  0 <= pend)  p[o * 0  +  0] |= BIT##0;\
	if (p + o * s1 + o1 <= pend)  p[o * s1 + o1] |= BIT##b1;\
	if (p + o * s2 + o2 <= pend)  p[o * s2 + o2] |= BIT##b2;\
	if (p + o * s3 + o3 <= pend)  p[o * s3 + o3] |= BIT##b3;\
	if (p + o * s4 + o4 <= pend)  p[o * s4 + o4] |= BIT##b4;\
	if (p + o * s5 + o5 <= pend)  p[o * s5 + o5] |= BIT##b5;\
	if (p + o * s6 + o6 <= pend)  p[o * s6 + o6] |= BIT##b6;
#endif

static void crossOffWheelFactor2(uchar* p, const uchar* pend, const uint step)
{
	const uint o = step / WHEEL30;
	switch (WheelInit30[step % WHEEL30].PrimeIndex)
	{
		case 0 :
			poSet(6,0,1, 10,0,2, 12,0,3,  16,0,4,  18,0,5,  22,0,6,  28,0,7);
			break;
		case 1 :
			poSet(4,0,7, 6,1,3,  10,2,2,  16,3,6,  18,4,1,  24,5,5,  28,6,4);
			break;
		case 2 :
			poSet(2,0,6, 6,2,1,  8,2,7,   12,4,3,  18,6,5,  20,7,2,  26,9,4);
			break;
		case 3 :
			poSet(4,1,6, 6,2,5,  10,4,2,  12,5,1,  16,6,7,  22,9,4,  24,10,3);
			break;
		case 4:
			poSet(6,3,3, 8,4,4,  14,7,7,  18,10,1, 20,11,2, 24,13,5, 26,14,6);
			break;
		case 5 :
			poSet(4,2,4, 10,6,2, 12,7,5,  18,11,3, 22,13,7, 24,15,1, 28,17,6);
			break;
		case 6 :
			poSet(2,1,4, 6,4,5,  12,9,1,  14,10,6, 20,15,2, 24,18,3, 26,19,7);
			break;
		case 7 :
			poSet(2,1,7, 8,7,6,  12,11,5, 14,13,4, 18,17,3, 20,19,2, 24,23,1);
			break;
	}
}

static void crossOff2Factor(uchar* ps0, uchar* ps1, const uchar* pend, const ushort smask, const uint p)
{
	const uchar masks1 = smask >> 8, masks0 = (uchar)smask;
	while (ps1 <= pend) {
		*ps1 |= masks1, ps1 += p;
		*ps0 |= masks0, ps0 += p;
	}
	if (ps0 <= pend)
		*ps0 |= masks0;
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
	BucketInfo.CurStock +=  MEM_BLOCK;
}

#define PI(x, r) x / log((double)x) * (1 + r / log((double)x))
static int setBucketInfo(const uint sieve_size, const uint sqrtp, const uint64 range)
{
	StockHead = NULL;
	memset(Bucket, 0, sizeof(Bucket));

	//assert (range >> 32 < sieve_size);
	BucketInfo.MaxBucket = range / sieve_size + 1;
	BucketInfo.LoopSize = sqrtp / (sieve_size / PATTERN_GAP) + 2; //add the first, last bucket
	BucketInfo.Log2Size = ilog(sieve_size / WHEEL30, 2);
	BucketInfo.SieveSize = (1 << BucketInfo.Log2Size) - 1;
	//assert(BucketInfo.Log2Size <= (32 - SI_BIT) && SI_BIT >= WP_BIT);

	return 0;
}

static int segmentedSieve2(uchar bitarray[], uint start, uint sieve_size);
static void setWheelMedium(const uint sieve_size, const uint maxp, const uint64 start)
{
	uint j = Threshold.L1Index;
	Threshold.L2Index = 0;
	const uint pix = PI(maxp, 1.2);
	MediumWheel = (WheelPrime*) malloc(sizeof(WheelPrime) * (pix + 100));
	uint msize = MIN(sieve_size, Threshold.L2Size);

	uchar bitarray[(L1_DCACHE_SIZE + L2_DCACHE_SIZE) << 10];
	//assert(Threshold.L1Maxp * Threshold.L1Maxp > Threshold.L2Size);
	//assert(Threshold.L1Maxp < Threshold.L2Maxp);

	for (uint l2size = L2_DCACHE_SIZE * WHEEL30 << 10, sieve_index = Threshold.L1Maxp; sieve_index < maxp; sieve_index += l2size) {
		if (l2size > maxp - sieve_index)
			l2size = maxp - sieve_index;

		const int bytes = segmentedSieve2(bitarray, sieve_index, l2size);

		stype mask = 0, *sbitarray = (stype*)bitarray;
		uint noffset = sieve_index - sieve_index % WHEEL30 - sizeof(mask) * WHEEL30;
		const uint pn = 2 + bytes / sizeof(mask);

		for (uint i = 0; i < pn; ) {
			if (mask == 0) {
				mask = ~sbitarray[i ++];
				noffset += sizeof(mask) * WHEEL30;
				continue;
			}

			const uint p = noffset + PRIME_OFFSET(mask); mask &= mask - 1;
			if (p >= Threshold.L2Maxp && Threshold.L2Index == 0) {
				Threshold.L2Maxp = p;
				Threshold.L2Index = j;
				msize = sieve_size;
			}

			uint64 offset = start;
			const uint64 p2 = (uint64)p * p;
			if (p2 - msize >= offset) { //overflow
				offset += (p2 - offset) / msize * msize;
			}

			//assert(p2 < offset + msize);
			uint sieve_index = p - (uint)(offset % p);
			if (p2 > offset) {
				sieve_index = (uint)(p2 - offset);
			}

			const int pi = WHEEL_INIT[p % WHEEL].PrimeIndex;
			const WheelFirst& wf = WHEEL_FIRST[(sieve_index + uint(offset % WHEEL)) % WHEEL][pi];
			sieve_index += wf.Correct * p;
			//assert(sieve_index / WHEEL30 < (-1u >> SI_BIT));
			MediumWheel[j].Wp = (p / WHEEL << SI_BIT) + pi;
			MediumWheel[j ++].Si = (sieve_index / WHEEL30 << SI_BIT) + wf.WheelIndex;
		}
	}

	MediumWheel[j + 0].Wp = UINT_MAX;
	MediumWheel[j + 1].Wp = UINT_MAX;
}

static void pushBucket(const uint sieve_index, const uint wp, const uchar wheel_index)
{
	const uint next_bucket = sieve_index >> BucketInfo.Log2Size;
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
	wheel->Si = (sieve_index & BucketInfo.SieveSize) << SI_BIT | wheel_index;
}

static void setWheelPrime(uint medium, uint sqrtp, const uint64 start, const uint64 range)
{
	uint nextp = 0; uint64 remp = 0;
	if (sqrtp < UINT_MAX) sqrtp ++; //watch overflow if sqrtp = 2^32 - 1

	const uint irange = (uint)((range >> 32) > WHEEL30 ? UINT_MAX : range / WHEEL30);
	uchar bitarray[(L1_DCACHE_SIZE + L2_DCACHE_SIZE) << 10];

	for (uint l2size = L2_DCACHE_SIZE * WHEEL30 << 10; medium < sqrtp; medium += l2size) {
		if (l2size > sqrtp - medium)
			l2size = sqrtp - medium;

		const int bytes = segmentedSieve2(bitarray, medium, l2size);
		stype mask = 0, *sbitarray = (stype*)bitarray; //little endian
		uint offset = medium - medium % WHEEL30 - sizeof(mask) * WHEEL30;
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

//			if (sieve_index > range) continue;
			const uint wp = (p / WHEEL210 << WP_BIT) + WheelInit210[p % WHEEL210].PrimeIndex;
			const uint module_wheel = sieve_index % WHEEL210;
			const WheelFirst& wf = WheelFirst210[module_wheel][wp % (1 << WP_BIT)];
#if X86_64
			sieve_index = (sieve_index + wf.Correct * (uint64)p) / WHEEL30;
#else
			sieve_index = (sieve_index / WHEEL210 + (wp >> WP_BIT) * wf.Correct) * (WHEEL210 / WHEEL30);
			sieve_index += (wf.Correct * (p % WHEEL210) + module_wheel) / WHEEL30;
#endif

			if (sieve_index > irange) continue;
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

static void sieveSmall0(uchar bitarray[], const uchar* pend, const uint p, uint sieve_index, ushort multiples)
{
	uchar* ppbeg[8];
	for (int i = 0; i < 8; i ++) {
		const uchar windex = WheelInit30[sieve_index % WHEEL30].WheelIndex;
		ppbeg[windex] = bitarray + sieve_index / WHEEL30;
		sieve_index += (multiples % 4) * 2 * p; multiples /= 4;
	}

	#define OR_ADD(n)   *ps##n |= BIT##n, ps##n += p
	#define CMPEQ_OR(n) if (ps##n <= pend) *ps##n |= BIT##n

	uchar* ps0 = ppbeg[0], *ps1 = ppbeg[1], *ps2 = ppbeg[2], *ps3 = ppbeg[3];
	while (ps3 <= pend) {
		OR_ADD(0); OR_ADD(1); OR_ADD(2); OR_ADD(3);
	}
	CMPEQ_OR(0); CMPEQ_OR(1); CMPEQ_OR(2);

	uchar* ps4 = ppbeg[4], *ps5 = ppbeg[5], *ps6 = ppbeg[6], *ps7 = ppbeg[7];
	while (ps7 <= pend) {
		OR_ADD(4); OR_ADD(5); OR_ADD(6); OR_ADD(7);
	}
	CMPEQ_OR(4); CMPEQ_OR(5);CMPEQ_OR(6);
}

static void sieveSmall1(uchar bitarray[], const uchar* pend, const uint p, uint sieve_index, ushort multiples)
{
	for (int i = 0; i < 4; i ++) {
		uchar* ps0 = bitarray + sieve_index / WHEEL30;
		const uchar masks0 = WheelInit30[sieve_index % WHEEL30].UnsetBit;
		sieve_index += (multiples % 4) * 2 * p; multiples /= 4;

		uchar* ps1 = bitarray + sieve_index / WHEEL30;
		const uchar masks1 = WheelInit30[sieve_index % WHEEL30].UnsetBit;
		sieve_index += (multiples % 4) * 2 * p; multiples /= 4;
		crossOff2Factor(ps0, ps1, pend, masks0 | (ushort)masks1 << 8, p);
	}
}

static void sieveSmall3(uchar bitarray[], const uchar* pend, const uint p, uint sieve_index, ushort multiples)
{
	for (int i = 0; i < 8; i ++) {
		const uchar mask = WheelInit30[sieve_index % WHEEL30].UnsetBit;
		if (mask == BIT0)
			break;
		bitarray[sieve_index / WHEEL30] |= mask;
		sieve_index += (multiples % 4) * 2 * p; multiples /= 4;
	}
	crossOffWheelFactor2(bitarray + sieve_index / WHEEL30, pend, p);
}

static void eratSieveL1(uchar bitarray[], const uint64 start, const uint sieve_size, uint maxp)
{
	const uchar* pend = bitarray + sieve_size / WHEEL30;
	if (start + sieve_size < ((uint64)maxp) * maxp) {
		maxp = isqrt(start + sieve_size) + 1;
	}

	for (uint p = Prime[FIRST_INDEX], j = FIRST_INDEX; p < maxp; p = Prime[++j]) {
		uint sieve_index = p - (uint)(start % p);
		if (start <= p) {
			sieve_index = p * p - (uint)start;
		}

		const WheelFirst& wf = WheelFirst30[sieve_index % WHEEL30][WheelInit30[p % WHEEL30].PrimeIndex];
		const ushort multiples = WHEEL_SKIP >> (wf.NextMultiple * 2);
		sieve_index += wf.Correct * p;
		sieveSmall0(bitarray, pend, p, sieve_index, multiples);
	}
}

//Pre-sieve multiples of small primes <= 19
static void preSieve(uchar bitarray[], const uint64 start, const uint sieve_size)
{
	const uint offset = (uint)(start % PRIME_PRODUCT) / WHEEL30;
	const uint bits = sieve_size / WHEEL30 * 8 + WheelInit30[sieve_size % WHEEL30].WheelIndex;
	const int bytes = (bits + 7) / 8, remains = sizeof(PreSieved) - offset;
	if (remains > bytes) {
		memcpy(bitarray, PreSieved + offset, bytes);
	} else {
		memcpy(bitarray, PreSieved + offset, remains);
		memcpy(bitarray + remains, PreSieved, bytes - remains);
	}

	//1 is not prime, pattern < WHEEL30 is prime
	if (start == 0) {
		bitarray[0] = BIT0;
	}
	//set the last byte bit 1
	bitarray[bits >> 3] |= ~((1 << (bits & 7)) - 1);
}

static void eratSieveSmall(uchar bitarray[], const uint64 start, const uint sieve_size)
{
//	#pragma omp parallel for if (sieve_size > 2 * Threshold.L1Size)
//	preSieve(bitarray, start, sieve_size);
	for (uint sieve_index = 0, l1size = Threshold.L1Size; sieve_index < sieve_size; sieve_index += l1size) {
		if (l1size + sieve_index > sieve_size)
			l1size = sieve_size - sieve_index;
		preSieve(bitarray + sieve_index / WHEEL30, start + sieve_index, l1size);
		eratSieveL1(bitarray + sieve_index / WHEEL30, start + sieve_index, l1size, Threshold.L1Maxp);
	}
}

#define SAFE_SET(n) \
	we##n = wdata##n[(int)we##n.WheelIndex]; \
	bitarray[(int)sieve_index##n] |= we##n.UnsetBit; \
	sieve_index##n += we##n.Correct + (we##n.NextMultiple) * wi##n

//sieve 1 medium prime from array
static void sieveMedium1(uchar bitarray[], const uint sieve_byte, WheelPrime* pwheel)
{
	uint& wheel = pwheel->Si;
	int sieve_index = (wheel >> SI_BIT) - sieve_byte;
	if (sieve_index >= 0) {
		wheel -= sieve_byte << SI_BIT;
		return;
	}

	const uint wi = pwheel->Wp >> SI_BIT;
	WheelElement* wdata = WHEEL_MAP[pwheel->Wp % (1 << SI_BIT)];
	WheelElement we; we.WheelIndex = wheel % (1 << SI_BIT);

	do {
		we = wdata[we.WheelIndex];
		bitarray[(int)sieve_index] |= we.UnsetBit;
		sieve_index += we.Correct + (we.NextMultiple) * wi;
	} while (sieve_index < 0);

	wheel = (sieve_index << SI_BIT) | we.WheelIndex;
}

//sieve 3 medium prime from array
inline static int sieveMedium3(uchar bitarray[], const uint sieve_byte, WheelPrime* pwheel)
{
	uint& wheel0 = pwheel[0].Si, wp0 = pwheel[0].Wp;
	uint wi0 = wp0 >> SI_BIT, sieve_index0 = (wheel0 >> SI_BIT) - sieve_byte;

	uint& wheel1 = pwheel[1].Si, wp1 = pwheel[1].Wp;
	uint wi1 = wp1 >> SI_BIT, sieve_index1 = (wheel1 >> SI_BIT) - sieve_byte;

	uint& wheel2 = pwheel[2].Si, wp2 = pwheel[2].Wp;
	uint wi2 = wp2 >> SI_BIT, sieve_index2 = (wheel2 >> SI_BIT) - sieve_byte;

	WheelElement* wdata0 = WHEEL_MAP[wp0 % (1 << SI_BIT)];
	WheelElement* wdata1 = WHEEL_MAP[wp1 % (1 << SI_BIT)];
	WheelElement* wdata2 = WHEEL_MAP[wp2 % (1 << SI_BIT)];

	WheelElement we0, we1, we2;
	we0.WheelIndex = wheel0 % (1 << SI_BIT);
	we1.WheelIndex = wheel1 % (1 << SI_BIT);
	we2.WheelIndex = wheel2 % (1 << SI_BIT);

	while ((int)sieve_index0 < 0) {
		SAFE_SET(0);
		if ((int)sieve_index1 < 0) { SAFE_SET(1); }
		if ((int)sieve_index2 < 0) { SAFE_SET(2); }
	}

//	while ((int)sieve_index0 < 0) { SAFE_SET(0); }
	while ((int)sieve_index1 < 0) { SAFE_SET(1); }
	while ((int)sieve_index2 < 0) { SAFE_SET(2); }

	wheel0 = sieve_index0 << SI_BIT | we0.WheelIndex;
	wheel1 = sieve_index1 << SI_BIT | we1.WheelIndex;
	wheel2 = sieve_index2 << SI_BIT | we2.WheelIndex;

	return 3;
}

static WheelPrime* sieveMediumW30(uchar bitarray[], const uint sieve_byte, const uint minp, WheelPrime* pwheel)
{
	for (uint wp = pwheel->Wp; wp < minp; wp = pwheel->Wp) {
		const uint wi = wp >> SI_BIT, pi = wp % (1 << SI_BIT);
		const uint p = wi * WHEEL + PATTERN[pi];
		uint sieve_index = pwheel->Si >> SI_BIT;
		WheelElement* wdata = WHEEL_MAP[pi], we;
		WheelElement* wheel = wdata + pwheel->Si % (1 << SI_BIT);

		uchar* pbitarray = bitarray + sieve_index;
		sieve_index += (sieve_byte - sieve_index) / p * p - sieve_byte;
		const uchar* pend = bitarray + sieve_byte;
		for (int i = 0; i < 4; i ++) {
			uchar* ps0 = pbitarray;
			ushort smask = wheel->UnsetBit;
			pbitarray += wheel->Correct + wheel->NextMultiple * wi;
			wheel = wdata + wheel->WheelIndex;

			uchar* ps1 = pbitarray;
			smask |= (ushort)wheel->UnsetBit << 8;
			pbitarray += wheel->Correct + wheel->NextMultiple * wi;
			we.WheelIndex = wheel->WheelIndex, wheel = wdata + wheel->WheelIndex;
			crossOff2Factor(ps0, ps1, pend, smask, p);
		}

		while ((int)sieve_index < 0) {
//			SAFE_SET();
			we = wdata[we.WheelIndex];
//			bitarray[sieve_index] |= we.UnsetBit;
			sieve_index += we.Correct + we.NextMultiple * wi;
		}
		pwheel ++->Si = sieve_index << SI_BIT | we.WheelIndex;
	}

	return pwheel;
}

//sieve prime multiples in [start, start + sieve_size) with medium algorithm
static void eratSieveMedium(uchar bitarray[], const uint64 start, const uint sieve_size, const uint whi, uint maxp)
{
	if (start + sieve_size < (uint64)maxp * maxp) {
		maxp = isqrt(start + sieve_size) + 1;
	}

	const uint sieve_byte = sieve_size / WHEEL30 + WheelInit30[sieve_size % WHEEL30].WheelIndex;
	WheelPrime* pwheel = MediumWheel + whi;

#if W30
	uint minp = MIN(maxp, sieve_byte * 45 / 100);
	minp = (minp / WHEEL << SI_BIT) + WHEEL_INIT[minp % WHEEL].PrimeIndex;
	pwheel = sieveMediumW30(bitarray, sieve_byte, minp, pwheel);
#endif

	maxp = (maxp / WHEEL << SI_BIT) + WHEEL_INIT[maxp % WHEEL].PrimeIndex;
	bitarray += sieve_byte;

	while (pwheel[2].Wp < maxp) {
		pwheel += sieveMedium3(bitarray, sieve_byte, pwheel);
	}
	if (pwheel[0].Wp < maxp)
		sieveMedium1(bitarray, sieve_byte, pwheel ++);
	if (pwheel[0].Wp < maxp)
		sieveMedium1(bitarray, sieve_byte, pwheel);
}

//sieve big prime from bucket
static void sieveBig1(uchar bitarray[], const uint sieve_size, const WheelPrime* cur_wheel)
{
	uint sieve_index = cur_wheel->Si, wp = cur_wheel->Wp;
	const WheelElement* wh = Wheel210[wp % (1 << WP_BIT)];
	const WheelElement* wheel = wh + sieve_index % (1 << SI_BIT);
	bitarray[sieve_index >>= SI_BIT] |= wheel->UnsetBit;
	sieve_index += wheel->Correct + wheel->NextMultiple * (wp >> WP_BIT);

	if (sieve_index < sieve_size) {
		wheel = wh + wheel->WheelIndex;
		bitarray[sieve_index] |= wheel->UnsetBit;
		sieve_index += wheel->Correct + wheel->NextMultiple * (wp >> WP_BIT);
	}
	pushBucket(sieve_index, wp, wheel->WheelIndex);
}

//sieve 2 big prime from bucket, 15% improvement
static void sieveBig2(uchar bitarray[], const uint sieve_size, const WheelPrime* cur_wheel)
{
	uint sieve_index1 = cur_wheel[0].Si, wp1 = cur_wheel[0].Wp;
	uint sieve_index2 = cur_wheel[1].Si, wp2 = cur_wheel[1].Wp;
	const WheelElement* wh1 = Wheel210[wp1 % (1 << WP_BIT)];
	const WheelElement* wh2 = Wheel210[wp2 % (1 << WP_BIT)];

	const WheelElement* wheel1 = wh1 + sieve_index1 % (1 << SI_BIT);
	const WheelElement* wheel2 = wh2 + sieve_index2 % (1 << SI_BIT);
	bitarray[sieve_index1 >>= SI_BIT] |= wheel1->UnsetBit;
	sieve_index1 += wheel1->Correct + wheel1->NextMultiple * (wp1 >> WP_BIT);

	bitarray[sieve_index2 >>= SI_BIT] |= wheel2->UnsetBit;
	sieve_index2 += wheel2->Correct + wheel2->NextMultiple * (wp2 >> WP_BIT);

	if (sieve_index1 < sieve_size) {
		wheel1 = wh1 + wheel1->WheelIndex;
		bitarray[sieve_index1] |= wheel1->UnsetBit;
		sieve_index1 += wheel1->Correct + wheel1->NextMultiple * (wp1 >> WP_BIT);
	}
	if (sieve_index2 < sieve_size) {
		wheel2 = wh2 + wheel2->WheelIndex;
		bitarray[sieve_index2] |= wheel2->UnsetBit;
		sieve_index2 += wheel2->Correct + wheel2->NextMultiple * (wp2 >> WP_BIT);
	}
	pushBucket(sieve_index1, wp1, wheel1->WheelIndex);
	pushBucket(sieve_index2, wp2, wheel2->WheelIndex);
}

//This implementation uses a sieve array with WHEEL210 numbers per byte and
//a modulo 210 wheel that skips multiples of 2, 3, 5 and 7.
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
	int primes = 0;
	const int oper = cmd ? cmd->Oper : COUNT_PRIME;
//	if (COUNT_PRIME == oper) {
		primes = countBit0sArray((uint64*)bitarray, bytes * sizeof(uint64));
//	} else if (COPY_BITS == oper) {
//		primes = bytes;
//		cmd->Data = bitarray; // use stack buffer !!!
//	}

	return primes;
}

//sieve prime multiples in [start, start + sieve_size)
static int segmentedSieve(uchar bitarray[], const uint64 start, const uint sieve_size)
{
	const uint sqrtp = isqrt(start + sieve_size);
	const uint medium = MIN(sqrtp, Threshold.Medium) + 1;
	const bool bsieve = sqrtp >= Threshold.L1Maxp;

	for (uint sieve_index = 0, l2size = Threshold.L2Size; sieve_index < sieve_size; sieve_index += l2size) {
		if (l2size + sieve_index > sieve_size)
			l2size = sieve_size - sieve_index;

		uchar* buffer = bitarray + sieve_index / WHEEL30;
		eratSieveSmall(buffer, start + sieve_index, l2size);
		if (bsieve) {
			eratSieveMedium(buffer, start + sieve_index, l2size, Threshold.L1Index, Threshold.L2Maxp);
		}
	}

	if (medium > Threshold.L2Maxp) {
		eratSieveMedium(bitarray, start, sieve_size, Threshold.L2Index, medium);
	}

	if (start >= Threshold.BucketStart) {
		eratSieveBig(bitarray, Config.SieveSize / WHEEL30);
	}

	return sieve_size / WHEEL30 + (7 + WheelInit30[sieve_size % WHEEL30].WheelIndex) / 8;
}

static int segmentedSieve2(uchar bitarray[], uint start, uint sieve_size)
{
	const uint sqrtp = isqrt(start + sieve_size);
	const uint start_align = (uint)(start % WHEEL30);
	start -= start_align, sieve_size += start_align;
	const uint bytes = sieve_size / WHEEL30 + (7 + WheelInit30[sieve_size % WHEEL30].WheelIndex) / 8;

	eratSieveSmall(bitarray, start, sieve_size);
	const uchar* pend = bitarray + sieve_size / WHEEL30;
	for (uint j = Threshold.L1Index, p = Threshold.L1Maxp; p <= sqrtp; p = Prime[++j]) {
		const uint sieve_index = p - start % p;
		const WheelFirst& wf = WheelFirst30[sieve_index % WHEEL30][WheelInit30[p % WHEEL30].PrimeIndex];
		const ushort multiples = WHEEL_SKIP >> (wf.NextMultiple * 2);
		sieveSmall0(bitarray, pend, p, sieve_index + wf.Correct * p, multiples);
	}

	bitarray[0] |= (1 << WheelInit30[start_align].WheelIndex) - 1;
	*(uint64*)(bitarray + bytes) = ~0;
	return bytes;
}

static void setThresholdL1()
{
	for (uint p = Prime[1], j = 1; ; p = Prime[++j]) {
		if (p >= Threshold.L1Maxp) {
			Threshold.L1Index = j;
			Threshold.L1Maxp = p;
			break;
		}
	}
}

void setCpuCache(int level, uint cache)
{
	cache = 1 << ilog(cache, 2);

	if (level == 1 && cache >= 16 && cache < L2_DCACHE_SIZE) {
		Threshold.L1Size = cache * (WHEEL30 << 10);
		Threshold.L1Maxp = Threshold.L1Size / (WHEEL30 * Threshold.L1Segs);
		setThresholdL1();
	} else if (level == 2 && cache >= L2_DCACHE_SIZE && cache <= MAX_SEGMENT) {
		Threshold.L2Size = cache * (WHEEL30 << 10);
		Threshold.L2Maxp = Threshold.L2Size / (WHEEL30 * Threshold.L2Segs);
	}
}

void setLevelSegs(uint level, uint segs)
{
	if (segs > 1) {
		if (level == 1) {
			Threshold.L1Segs = segs;
			Threshold.L1Maxp = Threshold.L1Size / (WHEEL30 * segs);
			setThresholdL1();
		} else if(level == 2) {
			Threshold.L2Segs = segs;
			Threshold.L2Maxp = Threshold.L2Size / (WHEEL30 * segs);
		} else if (ERAT_BIG > 2 && segs <= 6) {
			Threshold.Msegs = segs;
		}
	}
}

int setSieveSize(uint sieve_size)
{
	if (sieve_size <= MAX_SEGMENT && sieve_size >= 32) {
		sieve_size = WHEEL30 << (ilog(sieve_size, 2) + 10);
	} else {
		sieve_size = L2_DCACHE_SIZE * WHEEL30 << 10;
	}

	return Config.SieveSize = sieve_size;
}

static int checkSmall(const uint64 start, const uint64 end, Cmd* cmd)
{
	int primes = 0;
	for (int i = 1, p = 2; i < 4; p = Prime[++i]) {
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

static int64 pi(uint64 start, uint64 end, Cmd* cmd)
{
	const int64 ts = getTime();
	uint64 primes = checkSmall(start, end, cmd);

	uint start_align = (uint)(start % WHEEL210);
	start -= start_align;

	if (++end == 0) end --; //watch overflow if end = 2^64-1
	const int64 range = (int64)(end - start);

	uchar* bitarray = (uchar*) malloc((Config.SieveSize + Threshold.L1Size) / WHEEL30);
	for (uint si = 0, sieve_size = Config.SieveSize; start < end; start += sieve_size) {
		if (sieve_size > end - start) {
			sieve_size = (uint)(end - start);
		}

		const uint bytes = segmentedSieve(bitarray, start, sieve_size);
		if (start_align > 0) {
			memset(bitarray, ~0, start_align / WHEEL30);
			bitarray[start_align / WHEEL30] |= (1 << WheelInit30[start_align % WHEEL30].WheelIndex) - 1;
			start_align = 0;
		}

		*(uint64*)(bitarray + bytes) = ~0;
		primes += segmentProcessed(bitarray, start, bytes, cmd);
#if 1
		if ((si ++ & Config.Progress) == 15) {
			double ratio = 100 - 100.0 * (int64)(end - start - sieve_size) / range;
			double timeuse = (getTime() - ts) / (10 * ratio);
			const uint64 picount = (int64)((100 * (int64)primes) / ratio);
			if (timeuse > 3600)
				printf("\r>> %.2f%%, sieve time ~= %.3f %s, primes ~= %llu\r", ratio, timeuse / 3600, "hour", picount);
			else
				printf("\r>> %.2f%%, sieve time ~= %.2f %s, primes ~= %llu\r", ratio, timeuse, "sec", picount);
		}
#endif
	}

	free(bitarray);
	return primes;
}

static void convertSci(uint64 n, char buff[32])
{
	const int logn = ilog(n, 10);
	const uint64 pown = ipow(10, logn);
	if (n % pown == 0)
		sprintf(buff, "%de%d", (int)(n / pown), logn);
	else if ((n & (n - 1)) == 0)
		sprintf(buff, "2^%d", ilog(n, 2));
}

static void printResult(const uint64 start, const uint64 end, uint64 primes)
{
	char buff[128] = {0};
	char begbuff[32] = {0}, endbuff[32] = {0}, rangebuff[32] = {0};
	if (end > 10000) {
		convertSci(start, begbuff);
		convertSci(end, endbuff);
	}
	if (end + 1 == 0)
		sprintf(endbuff, "2^64");

	const int64 range = end - start;
	if (range > 10000)
		convertSci(range, rangebuff);

	if (start == 0) {
		if (endbuff[0] == 0)
			sprintf(buff, "%llu", end);
		else
			sprintf(buff, "%s", endbuff);
	}
	else if (begbuff[0] && endbuff[0])
		sprintf(buff, "%s,%s", begbuff, endbuff);
	else if (rangebuff[0]) {
		if (begbuff[0])
			sprintf(buff, "%s,%s+%s", begbuff, begbuff, rangebuff);
		else if (endbuff[0])
			sprintf(buff, "%s-%s,%s", endbuff, rangebuff, endbuff);
		else
			sprintf(buff, "%llu,%llu+%s", start, start, rangebuff);
	} else {
		if (begbuff[0] == 0)
			sprintf(begbuff, "%llu", start);
		if (endbuff[0] == 0)
			sprintf(endbuff, "%llu", end);
		sprintf(buff, "%s,%s", begbuff, endbuff);
	}
	printf("\rpi(%s) = %llu", buff, primes);
}

static uint setBucketStart(const uint64 start, const uint sqrtp, const uint sieve_size)
{
	uint medium = sieve_size / Threshold.Msegs + 1;
	uint64 offset = (uint64)medium * medium;
	//adjust medium
	if (sqrtp > medium && offset > start) {
		offset += sieve_size - offset % sieve_size + start % sieve_size;
		while (offset % WHEEL210 != start % WHEEL210) offset += sieve_size;
		medium = isqrt(offset - offset % WHEEL210 + sieve_size) + 1;
	} else {
		offset = start;
	}

	medium = MIN(sqrtp, medium);

	Threshold.BucketStart = offset - offset % WHEEL210;
	Threshold.Medium = medium;

	return medium;
}

uint64 doSieve(const uint64 start, const uint64 end, Cmd* cmd = NULL)
{
	const int64 ts = getTime();
	const uint sqrtp = isqrt(end);
	uint& sieve_size = Config.SieveSize;
	if (sieve_size < Threshold.L2Size && sqrtp > 10000000) {
		setSieveSize(SIEVE_SIZE);
	} else if (sieve_size >= Threshold.Msegs * (WHEEL30 << 21)) {
		setSieveSize(SIEVE_SIZE);
	}

	//init medium sieve
	const uint medium = setBucketStart(start, sqrtp, sieve_size);
	if (sqrtp >= Threshold.L1Maxp) {
		setWheelMedium(sieve_size, medium + 255, start - start % WHEEL210);
	}

	memset(&BucketInfo, 0, sizeof(BucketInfo));
	//init bucket sieve
	uint64& buckets = Threshold.BucketStart;
	if (sqrtp > medium) {
		setBucketInfo(sieve_size, sqrtp, end - buckets);
		setWheelPrime(medium, sqrtp, buckets, end - buckets);
		allocWheelBlock(1);
	} else {
		buckets = end + 1;
	}

	const int64 ti = getTime();
	const int64 primes = pi(start, end, cmd);

	if (MediumWheel)
		free(MediumWheel), MediumWheel = NULL;
	for (uint i = 0; i < BucketInfo.PoolSize; i ++)
		free(WheelPool[i]);

	if (Config.Flag & PRINT_RET) {
		const int64 ta = getTime();
		//printf("\rpi(%llu, %llu) = %llu", start, end, primes);
		printResult(start, end, primes);
		if (Config.Flag & PRINT_TIME)
			printf(" (%.2f + %.2f(init) = %.2f sec %d)", (ta - ti) / 1000.0, (ti - ts) / 1000.0, (ta - ts) / 1000.0, Config.SieveSize / (30 << 10));
		putchar('\n');
	}
#ifndef B_R
	assert (BucketInfo.StockSize == BucketInfo.CurStock);
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
		for (int sieve_index = p; sieve_index < PRIME_PRODUCT; sieve_index += p * 2) {
			PreSieved[sieve_index / WHEEL30] |= WheelInit30[sieve_index % WHEEL30].UnsetBit;
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
			int multiples = 0, sieve_index = i;
			if (i % 2 == 0) {
				multiples = 1;
				sieve_index += Pattern30[pi];
			}
			while (WheelInit30[sieve_index % WHEEL30].UnsetBit == 0) {
				sieve_index += Pattern30[pi] * 2;
				multiples += 2;
			}
			int wi = WheelInit30[sieve_index % WHEEL30].WheelIndex;
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
	const uint64 maxrange = ipow(2, 30);
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

		setSieveSize(rand() % MAX_SEGMENT + 64);
		setLevelSegs(1, rand() % 6 + 2), setLevelSegs(2, rand() % 6 + 2), setLevelSegs(3, rand() % 5 + 2);

		primes += doSieve(lowerBound, end - 1);
		if (lowerBound % 8 == 0)
		printf("chunk: %.2f%%\r", 100 - (int64)(upperBound - lowerBound) * 100.0 / range);
		lowerBound = end;
	}
	if (lowerBound + 1 == 0)
		printf("Pi[2^64-10^%d, 2^64-1] = %llu\n", rlog10, primes);
	else
		printf("Pi[10^%d, 10^%d+10^%d] = %llu\n", llog10, llog10, rlog10, primes);

	assert(primes == Ret);
}

static void startBenchmark()
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
	puts("All Big tests passed SUCCESSFULLY!");
	const uint64 pow11 = ipow(10, 11);
	const uint64 pow12 = pow11 * 10, pow9 = pow11 / 100;

	const uint64 RangeData[][3] =
	{
		{-1 - pow11,   pow11, 2254197466ul},
		{ipow(01, 11), pow11, 4118054813ul},
		{ipow(10, 12), pow11, 3612791400ul},
		{ipow(10, 13), pow11, 3340141707ul},
		{ipow(10, 17), pow11, 2554712095ul},

		{-1 - pow12,   pow12, pow9 * 22 + 542106206},
		{ipow(10, 14), pow12, pow9 * 31 + 16203073},
		{ipow(10, 15), pow12, pow9 * 28 + 952450479},
		{ipow(10, 16), pow12, pow9 * 27 + 143405794},
		{ipow(10, 18), pow12, pow9 * 24 + 127637783},
		{ipow(10, 19), pow12, pow9 * 22 + 857444126},
		{ipow(10, 19), pow12, pow9 * 22 + 857444126},
		{-1 - pow12*10, pow12*10, pow9 * 225 + 420940155},
	};

	for (int k = 0; k < sizeof(RangeData) / sizeof(RangeData[0]); k ++) {
		fixRangeTest(RangeData[k][0], RangeData[k][1], RangeData[k][2]);
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
	"g++ -DSIEVE_SIZE=2048 -DW30 -march=native -funroll-loops -O3 -s -pipe\n");

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
	info += sprintf(info, "[ARGS ] : L1Size = %d, L2Size = %d, SieveSize = %d, Loop = %d\n",
			Threshold.L1Size / WHEEL30 >> 10, Threshold.L2Size / WHEEL30 >> 10, Config.SieveSize / WHEEL30 >> 10, BucketInfo.LoopSize);
	info += sprintf(info, "[ARGS ] : L1Seg/L2Seg/Bseg/Mseg = (%d,%d,%d,%d)\n",
		Threshold.L1Segs, Threshold.L2Segs, 2, Threshold.Msegs);
	info += sprintf(info, "[ARGS ] : L1Maxp/L2Maxp/Medium/Large/SieveSize = (%d,%d,%d,%u,%d)",
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
		char c = params[i][0];
		if (c >= 'a' && c <= 'z')
			c += 'A' - 'a';
		if (isdigit(c) || c == 'E') {
			if (cmdi < 0)
				cmdi = i;
			continue;
		}
		if (isdigit(params[i][1]))
			cdata = atoi(params[i] + 1);

		switch (c)
		{
			case 'H' :
				puts(Benchmark);
				puts(Help);
				break;
			case 'S':
				setSieveSize(cdata);
				break;
			case 'C':
				setCpuCache(cdata % 10, cdata / 10);
				break;
			case 'L':
				setLevelSegs(cdata / 10, cdata % 10);
				break;
			case 'M':
				Config.Progress = (1 << cdata) - 1;
				break;
			case 'D':
				Config.Flag ^= (1 << (toupper(params[i][1]) - 'A'));
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
		char c = *ccmd;
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

		char cmdc = toupper(params[cmdi][0]);
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
					uint64 start = ipow(10, j);
					doSieve(start, start + range);
				}
			}
			if (isdigit(params[cmdi + 2][0])) {
				int powi = atoi(params[cmdi + 2]);
				uint64 range = powi > 12 ? ipow(2, powi) : ipow(10, powi);
				for (int i = 32; i < 64 && powi > 0; i ++) {
					uint64 start = ipow(2, i);
					doSieve(start, start + range);
				}
			} else
				startBenchmark();
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
	if (argc > 1)
		executeCmd(argv[1]);

#if 0
	if (argc == 2)
	{
	srand(time(0));
	initPrime(SIEVE_SIZE);

	Config.Progress = 0;
//	Config.Flag &= ~PRINT_TIME;

	for (int i = 20; i > 0; i --)
	for (int j = 1; j < 10; j ++) {
		uint64 start = ((uint64)(rand() * rand()) << i) + (uint64)rand() * rand() * rand();
		int sieve_size = setSieveSize(MAX_SEGMENT / (rand() % 8 + 1));
		uint64 range = (rand() % 256) * Config.SieveSize + rand();
		uint64 rm1 = doSieve(start, start + range, NULL);
		if (setSieveSize(MAX_SEGMENT / (rand() % 4 + 1)) == sieve_size)
		{
			setSieveSize(sieve_size*2);
		}

		Config.Flag ^= PRINT_RET;
		setLevelSegs(3, rand() % 6 + 2); setLevelSegs(2, rand() % 6 + 2), setLevelSegs(1, rand() % 6 + 2);
		uint64 rm2 = doSieve(start, start + range, NULL);
		if(rm1 != rm2) { printInfo(); printf("%llu %llu\n", start, range); system("pause"); }
		if (j % 100 == 0) printf("\r %2d progress ~= %d\n", i, j);
		Config.Flag ^= PRINT_RET;
	}

	Config.Flag &= ~PRINT_TIME;
	Config.Flag ^= PRINT_RET;

	for (int j = 1; j <= 1000; j ++) {
		uint64 start = (uint64)(rand() * rand()) * (uint64)(rand() * rand());
		setSieveSize(SIEVE_SIZE);
		start -= start % WHEEL210;
		uint64 range = (rand() % 10 + 10) * Config.SieveSize;
		uint64 rm1 = doSieve(start, start + range - 1, NULL);

		setSieveSize(MAX_SEGMENT / (rand() % 16 + 1));
		uint64 rm2 = doSieve(start, start + range, NULL);
		assert(rm1 == rm2);
		for (int i = 2; i <= 6; i ++) {
			const uint64 medium = Config.SieveSize / i + 1;
			//const uint64 start = (uint64)(rand() * rand()) * (uint64)(rand() * rand());
			//const uint64 range = ((uint)rand() * rand()) % (1 << 27);
			uint64 start = medium * medium - Config.SieveSize * (rand() % 8 + 1) + rand();
			uint64 end = medium * medium + Config.SieveSize * (rand() % 8 + 1) + rand();

			setLevelSegs(3, i+1); setLevelSegs(2, rand() % 6 + 2), setLevelSegs(1, rand() % 6 + 2);
			const uint64 r1 = doSieve(start, end, NULL);

			Config.Flag ^= PRINT_RET;
			setSieveSize(Config.SieveSize * 2);
			setLevelSegs(3, i); setLevelSegs(2, rand() % 6 + 2), setLevelSegs(1, rand() % 6 + 2);
			const uint64 r2 = doSieve(start, end, NULL);
			if (r1 != r2) {
				printf(" --- %lld != %lld, %lld %lld c%d L1%d L2%d L3%d\n",
					r1, r2, start, end, Config.SieveSize, Threshold.L1Segs, Threshold.L2Segs, Threshold.Msegs);
			}
		}
		if (j % 10 == 0) printf("\r progress ~= %d\n", j);
	}}
#endif

	executeCmd("e9; e16 e10; e14 e12");

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
//2.improve bucket algorithm according to current segment, for ex[1e14, 1e16]
