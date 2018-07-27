/*Fast Single-thread segmented sieve of Eratosthenes prime number n < 2^64 ***/
static const char* Benchmark =
"g++ -DSIEVE_SIZE=2048 -DFSL1 -DFDIV -march=native -funroll-loops -O3 -s -pipe\n"
"Windows 10 x64               i3-350M,i5-3470,i7-7500u,i7-6700,r7-1700\n"
"Pi(0,    1e10) = 455052511    3.10   1.84    1.55     1.40    1.65\n"
"Pi(1e11, 1e10) = 394050419    4.35   2.50    2.02     1.83    2.00\n"
"Pi(1e12, 1e10) = 361840208    5.40   3.00    2.40     2.14    2.30\n"
"Pi(1e13, 1e10) = 334067230    6.60   3.50    2.85     2.50    2.70\n"
"Pi(1e14, 1e10) = 310208140    7.90   4.20    3.50     3.00    3.20\n"
"Pi(1e15, 1e10) = 289531946    10.2   5.10    4.32     3.67    3.91\n"
"Pi(1e16, 1e10) = 271425366    12.1   6.10    5.11     4.34    4.75\n"
"Pi(1e17, 1e10) = 255481287    14.4   7.09    5.95     5.17    5.75\n"
"Pi(1e18, 1e10) = 241272176    17.5   8.58    7.17     6.25    7.16\n"
"Pi(1e19, 1e10) = 228568014    24.0   11.6    9.86     8.64    9.82\n"
"Pi(0-1e9,10^9) = 22537866     8.15   4.28    3.92     3.42    3.82\n"
"Pi(1e18, 10^6) = 24280        0.65   0.46    0.34     0.52    0.60\n"
"Pi(1e18, 10^8) = 2414886      1.30   0.81    0.70     0.64    0.80\n"
"Pi(1e18, 10^9) = 24217085     3.58   1.80    1.58     1.40    1.55\n"
"Pi(0,    1e12) = 37607912018  500    270     224      200     220\n"
"Pi(1e14, 1e12) = 31016203073  790    420     354      295     320\n"
"Pi(1e16, 1e12) = 27143405794  1200   600     512      430     485\n"
"Pi(1e18, 1e12) = 24127637783  1505   760     622      560     640\n"
"Pi(1e19, 1e12) = 22857444126  1700   830     702      610     665\n";

static const char* Help = "\
	[B: Benchmark (0 - 12, 0 - 40)]\n\
	[D: D[T, R] dump time and result]\n\
	[M: Progress of calculating (0 - 20)]\n\
	[C: Cpu L1/L2 data cache size (L1:16-64, L2:256-4096)k]\n\
	[S: Set the sieve size (32 - 4096)k]\n\
	[L: Set sieve cache segs L(2-12)1, L(2-8)2 L(2-6)3]\n\
	[I: Info of programming]\n\
	[P: Print prime in [start, end]]\n\n\
Example:\n\
	1e16 10^10 s1024 c321 c2562";

#include <ctype.h>
#include <math.h>
#include <memory.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#if _MSC_VER >= 1400
#include <intrin.h>
#endif

#if __x86_64__ || __amd64__ || _M_X64 || __amd64 || __x86_64
	# define X86_64       1
#endif

#if _M_IX86 | __i386 | __i386__ | _X86_
	# define X86          1
#endif

#if _WIN32 && _MSC_VER < 1500
	typedef unsigned __int64 uint64;
	typedef __int64 int64;
#else
	typedef unsigned long long uint64;
	typedef long long int64;
#endif

typedef unsigned int uint;
typedef unsigned char uchar;
typedef unsigned short ushort;

#if X86_64 && _MSC_VER
	# define ASM_X86      0
//	#ifndef __clang_version__
	# define BIT_SCANF    1
//	#endif
#elif X86_64 || X86
//#ifndef __clang_version__
	# define ASM_X86      1
	# define BIT_SCANF    1
//#endif
#else
	# define ASM_X86      0
	# define BIT_SCANF    0
#endif

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

#ifndef W210 //big sieve wheel
	# define WHEEL        WHEEL30
	# define WHEEL_DATA   WheelData30
	# define WHEEL_INIT   WheelInit30
	# define WHEEL_FIRST  WheelFirst30
#else
	# define WHEEL        WHEEL210
	# define WHEEL_DATA   WheelData210
	# define WHEEL_INIT   WheelInit210
	# define WHEEL_FIRST  WheelFirst210
#endif

#ifndef __cplusplus
	#define bool   int
	#define true   1
	#define false  0
	typedef struct WheelElement  WheelElement;
	typedef struct WheelInit     WheelInit;
	typedef struct SievePrime    SievePrime;
	typedef struct BucketInfo_   BucketInfo_;
	typedef struct Stock         Stock;
	typedef struct Bucket_       Bucket_;
	typedef struct SmallSieve    SmallSieve;
	typedef struct PrimeCall     PrimeCall;
#endif

#define MIN(a, b)         (a < b ? a : b)
#define MAX(a, b)         (a > b ? a : b)

enum ECONST
{
	ERAT_SMALL    = 2, //4 - 16
	ERAT_MEDIUM   = 2, //2 - 6

#ifndef ERAT_BIG
# define ERAT_BIG   5 //2 - 6
#endif

	WHEEL30       = 30,
	WHEEL210      = 210,
	PRIME_PRODUCT = 210 * 11 * 13 * 17 * 19,
	FIRST_INDEX   = PRIME_PRODUCT / 9699690 + 7,
	NEXT_MULTIPLE = 0x5A28A6, //magic number from 0x799b799b. 4 bit shrift into 3 bit
	MAX_WHEEL_GAP = 11 - 1,   //max prime wheel 210 gap
	PI_65536      = 6542 + 1, //pi(2^16) + 1
	MAX_SEGMENT   = 4 << 10,  //4Mkb
#ifndef UINT_MAX
	UINT_MAX      = 0-1u,
#endif

	SP_BIT        = 6, //48 < 2^SP_BIT < WHEEL210
	SI_BIT        = 8, //[8 - 10]

#if (L2_DCACHE_SIZE < 256 || L2_DCACHE_SIZE > 4096)
	L2_DCACHE_SIZE= 256,
#endif
#ifndef SIEVE_SIZE
	SIEVE_SIZE    = 2048,
#endif
};

enum EBUCKET
{
	UINT_PIMAX = 203280221, //= pi(2^32)
	MAX_BUCKET = 5465,      //= 0xffffffff / (256 * (1 << 10) * 3) + 4, 5465 (sqrtp * (8 + 2) / sieve_size + 2);
	WHEEL_SIZE = 1 << 11,   //= 2048 11: 16k, [10 - 13]
	MEM_WHEEL  = WHEEL_SIZE * sizeof(int) * 2, //=32768
	MEM_BLOCK  = (1 << 19) / WHEEL_SIZE, //19:4 MB
	MAX_STOCK  = UINT_PIMAX / WHEEL_SIZE + MAX_BUCKET + MEM_BLOCK, //=55221
	MAX_POOL   = UINT_PIMAX / (MEM_BLOCK * WHEEL_SIZE) + 100 //=487
};

enum EFLAG
{
	PRINT_TIME = 1 << ('T' - 'A'),
	SLOW_TEST  = 1 << ('A' - 'A'),
	SAVE_DATA  = 1 << ('F' - 'A'),
	PRINT_RET  = 1 << ('R' - 'A'),
};

enum EBITMASK
{
	BIT0 = 1 << 0, BIT1 = 1 << 1,
	BIT2 = 1 << 2, BIT3 = 1 << 3,
	BIT4 = 1 << 4, BIT5 = 1 << 5,
	BIT6 = 1 << 6, BIT7 = 1 << 7,
};

struct Threshold_
{
	uint L1Size;
	uint L2Size;

	uint L1Maxp;
	uint L1Index;

	uint L2Maxp;
	uint L2Index;

	uint Medium;

	uint L3Maxp;
	uint L3Index;
	uint64 BucketStart; //min bucket start
};

struct Config_
{
	uint SieveSize;
	uint Progress;
	uint Flag;

	uint L1Segs;
	uint L2Segs;
	uint Msegs;
};

struct WheelElement
{
	char WheelIndex;
	uchar MaskBit;
	uchar Correct;
	uchar Multiple;
};

struct WheelInit
{
	char WheelIndex;
	uchar MaskBit;
};

struct SmallSieve
{
	uint Prime;
	uint Multiple;
	uint Si;
};

struct SievePrime
{
	//[0 - 7]: sieve_index % wheel, [8 - 31]: sieve_index / WHEEL30
	uint Si;
	//[0 - 6]: p % wheel, [7 - 31]: p / wheel
	uint Sp;
};

struct BucketInfo_
{
	uint Log2Size;
	uint MaxBucket;
	uint LoopSize;
	uint ModuloSize;

	uint CurStock;
	uint StockSize;
	uint PoolSize;
};

struct Stock
{
	SievePrime* Sprime; //read only
	Stock* Next;
};

struct Bucket_
{
	SievePrime* Sprime; //write only
	Stock* Head; //Head->s1->s2->....->sn
};

//api start
typedef void (*callback)(void* data, const uint64 prime);
struct PrimeCall
{
	void* data;
	callback func;
	uint64 primes;
};

#if 0
void printInfo();
void initPrime(int sieve_size);
//uint64 doSieve(const uint64 start, const uint64 end, PrimeCall* pcall);
uint setSieveSize(uint sieve_size);
void setCacheSegs(uint level, uint cachesegs);
void setCacheSize(uint level, uint cachecpu);
#endif

//common cache
typedef WheelElement WheelFirst;
static uchar Pattern30[64], Pattern210[48];
static WheelFirst WheelFirst30[WHEEL30][8];
static WheelFirst WheelFirst210[WHEEL210][48];

static WheelInit WheelInit30[WHEEL30];
static WheelInit WheelInit210[WHEEL210];

static WheelElement WheelData30[8][8];
static WheelElement WheelData210[48][48];
//presieved with prime <= 19
static uchar PreSieved[PRIME_PRODUCT / WHEEL30];

#if BIT_SCANF == 0
static uchar Lsb[1 << 16];
#endif

#if POPCNT == 0
//number of bits 1 binary representation table in Range[0-2^16)
static uchar WordNumBit1[1 << 16];
#else
# include <popcntintrin.h>
#endif

////////////////thread/task data //////////////////////////////
/**********************************
stock pool for big sieve
bucket.array[n]
|----
|3|.|----     ----------------
|2|2| | |---> |s1->s2->...->sn|---> wheel.array(4k size)
|1|1|.|n|     ----------------
---------       stock.list
b1.b2.bn
****************************/

static BucketInfo_ BucketInfo;
static Stock* StockHead;
static Stock  StockCache [MAX_STOCK];

//big/bucket wheel pool
static Bucket_   Bucket [MAX_BUCKET];
static SievePrime* WheelPool [MAX_POOL];

//medium pool
static SievePrime* MediumPrime;

//small  pool
static SmallSieve SmallPrime[PI_65536];

//config
static struct Threshold_ Threshold =
{
	32 << 10, L2_DCACHE_SIZE << 10,
	(32 << 10) / ERAT_SMALL, 0,
	(L2_DCACHE_SIZE << 10) / ERAT_MEDIUM, 0,
	SIEVE_SIZE * (WHEEL30 << 10) / ERAT_BIG,
	100000000,
};

static struct Config_ Config =
{
	SIEVE_SIZE << 10,
	(1 << 6) - 1,
	PRINT_RET | PRINT_TIME,
	ERAT_SMALL,
	ERAT_MEDIUM,
	ERAT_BIG,
};

////////////////basic func //////////////////////////////

//the simple sieve of Eratosthenes implementated by bit packing
static int eratoSimple()
{
	int primes = 1;
	const uint maxp = (1 << 16) + 1;
	uchar bitarray[1 << 12] = {0};

	for (uint p = 3; p < maxp; p += 2) {
		if (0 == (bitarray[p >> 4] & (1 << (p / 2 & 7)))) {
			SmallPrime[primes ++].Prime = p;
			for (uint j = p * p; j <= maxp; j += p * 2)
				bitarray[j >> 4] |= 1 << (j / 2 & 7);
		}
	}

	SmallPrime[primes ++].Prime = maxp;
	return primes;
}

//The first presieved template, cross off the first 8th prime multiples
static void initBitTable()
{
	uint i = 0;
#if POPCNT == 0
	const uint nbitsize = sizeof(WordNumBit1) / sizeof (WordNumBit1[0]);
	for (i = 1; i < nbitsize; i ++)
		WordNumBit1[i] = WordNumBit1[i >> 1] + (i & 1);
#endif

	const uchar pattern[] = {1, 7, 11, 13, 17, 19, 23, 29};
	for (i = 0; i < sizeof(Pattern30) / sizeof (Pattern30[0]); i ++)
		Pattern30[i] = pattern[i % 8] + WHEEL30 * (i / 8);

#if BIT_SCANF == 0
	const uint nbitsize2 = sizeof(Lsb) / sizeof(Lsb[0]);
	for (i = 0; i < nbitsize2; i += 2)
		Lsb[i + 0] = Lsb[i >> 1] + 1;

	for (i = 0; i < nbitsize2; i += 2) {
		Lsb[i + 0] = Pattern30[Lsb[i]];
		Lsb[i + 1] = 1;
	}
#endif

	for (int j = 0, wj = 0; j < WHEEL30; j ++) {
		WheelInit30[j].MaskBit = 0;
		WheelInit30[j].WheelIndex = wj;
		if (j == Pattern30[wj])
			WheelInit30[j].MaskBit = 1 << (wj ++);
	}

	for (i = 1; PRIME_PRODUCT % pattern[i] == 0; i ++) {
		const int p = pattern[i];
		for (uint offset = p; offset < sizeof(PreSieved) * WHEEL30; offset += p * 2)
			PreSieved[offset / WHEEL30] |= WheelInit30[offset % WHEEL30].MaskBit;
	}
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

	const int psize = 8;
	for (int i = 0; i < WHEEL30; i ++) {
		for (int pi = 0; pi < psize; pi ++) {
			int multiples = 1 - (i % 2);
			int next = i + multiples * Pattern30[pi];

			while (WheelInit30[next % WHEEL30].MaskBit == 0)
				next = i + Pattern30[pi] * (multiples += 2);

			int wi = WheelInit30[next % WHEEL30].WheelIndex;
			WheelElement* wf = &WheelFirst30[i][pi];
			wf->Multiple = (nextMultiple[wi] >> (pi * 4)) & 15;
			wf->WheelIndex = wi;
			wf->Correct = multiples;
//			wf.MaskBit = 1 << wi;
		}
	}

	for (int wi = 0; wi < psize; wi ++) {
		for (int pi = 0; pi < psize; pi ++) {
			int multiples = 2;
			int next = Pattern30[wi] + Pattern30[pi] * multiples;
			while (WheelInit30[next % WHEEL30].MaskBit == 0)
				next = Pattern30[wi] + Pattern30[pi] * (multiples += 2);

			WheelElement* we30 = &WheelData30[pi][wi];
			we30->Multiple = multiples;
			we30->Correct = next / WHEEL30 - Pattern30[wi] / WHEEL30;
			we30->WheelIndex = WheelInit30[next % WHEEL30].WheelIndex;
			we30->MaskBit = WheelInit30[Pattern30[wi]].MaskBit;
		}
	}
}

static void initWheel210()
{
	int wi = 0, i = 0;
	const int psize = sizeof(Pattern210) /sizeof(Pattern210[0]);
	uchar wpattern[WHEEL210];

	for (i = 0; i < WHEEL210; i ++) {
		const uchar mask = WheelInit30[i % WHEEL30].MaskBit;
		WheelInit210[i].MaskBit = mask;
		WheelInit210[i].WheelIndex = wi;
		if (mask && i % (WHEEL210 / WHEEL30)) {
			Pattern210[wi] = wpattern[wi] = i;
			wi ++;
		}
	}

	for (i = 0; i < WHEEL210; i ++) {
		for (int pi = 0; pi < psize; pi ++) {
			int multiples = 1 - (i % 2);
			int next = i + wpattern[pi] * multiples;

			WheelElement* wf = &WheelFirst210[i][pi];
			while (WheelInit210[next % WHEEL210].WheelIndex == WheelInit210[(next + 1) % WHEEL210].WheelIndex)
				next = i + wpattern[pi] * (multiples += 2);

			wf->WheelIndex = WheelInit210[next % WHEEL210].WheelIndex;
			wf->Correct = multiples;
		}
	}

	for (int pi = 0; pi < psize; pi ++) {
		for (wi = 0; wi < psize; wi ++) {
			int multiples = 2;
			int next = wpattern[wi] + wpattern[pi] * 2;

			while (WheelInit210[next % WHEEL210].WheelIndex == WheelInit210[(next + 1) % WHEEL210].WheelIndex) {
				next += wpattern[pi] * 2;
				multiples += 2;
			}

			WheelElement* we210 = &WheelData210[pi][wi];
			we210->Correct = next / WHEEL30 - wpattern[wi] / WHEEL30;
			we210->MaskBit = WheelInit210[wpattern[wi]].MaskBit;
			we210->WheelIndex = WheelInit210[next % WHEEL210].WheelIndex;
			we210->Multiple = multiples * (WHEEL210 / WHEEL30);
		}
	}
}

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
	while (n --)
		pown *= x;

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
	 stype index;
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
#elif __GNUC__
	return __builtin_ffsll(n) - 1;
	#if X86_64
	#else
	return __builtin_ffsl(n) - 1;
	#endif
#else
	stype pos = 0;
	#if X86_64
	__asm__ ("bsfq %1, %0" : "=r" (pos) : "rm" (n));
	return (uint)pos;
#else
	__asm__ ("bsfl %1, %0\n" : "=r" (pos) : "rm" (n) : "cc");
	return pos;
	#endif
#endif
}
#endif

inline static uint fastDiv(const uint64 n, uint p)
{
#if ASM_X86 == 0
	p = (uint)(n % p);
#elif __GNUC__ || __TINYC__
	//(n / p) < 2^32 is more efficient
	const uint loww = (uint)n, higw = (uint)(n >> 32);
	__asm__ ( "divl %%ecx\n" : "=d" (p) : "d"(higw), "a"(loww), "c"(p));
#else
	const uint loww = (uint)n, higw = (uint)(n >> 32);
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

//valid is [number][e][+-*^][number]
//ex: 1000, 2e9, 300-2e1, 2^32*2-e5 2^30-1E2, 2e9+2^20
static uint64 atoint64(const char* str)
{
	uint64 ret = 0;
	while (isdigit(*str))
		ret = ret * 10 + *str ++ - '0';

	if (*str && isdigit(str[1])) {
		if (str[0] == '^') {
			ret = ipow((uint)ret, atoi(str + 1));
		} else if (str[0] == 'e' || str[0] == 'E') {
			if (ret == 0)
				ret = 1;
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

//count number of bit 0 in binary representation of array, bytes % 8 = 0
static uint countBit0sArray(const uint64 bitarray[], const uint bytes)
{
	uint bit1s = 0, bit2s = 0;
	uint loops = bytes / sizeof(uint64);
	while (loops -- > 0) {
		const uint64 qw = *bitarray++;
#if POPCNT
		bit2s += _mm_popcnt_u64(qw);
#else
		const uint hig = (uint)(qw >> 32);
		const uint low = (uint)(qw);
		bit1s += WordNumBit1[(ushort)hig] + WordNumBit1[hig >> 16];
		bit2s += WordNumBit1[(ushort)low] + WordNumBit1[low >> 16];
#endif
	}

	return bytes * sizeof(uint64) - bit1s - bit2s;
}

//////////////////////// core code  //////////////////////////////////
static void allocWheelBlock(const uint blocks)
{
	SievePrime *pSprime = (SievePrime*) malloc((MEM_BLOCK + 1) * MEM_WHEEL);
	WheelPool[BucketInfo.PoolSize ++] = pSprime;
	//assert (BucketInfo.PoolSize < sizeof(WheelPool) / sizeof(WheelPool[0]));
	//assert (BucketInfo.StockSize + blocks < sizeof(StockCache) / sizeof(StockCache[0]));

	//align by MEM_WHEEL
	pSprime = (SievePrime*)((size_t)pSprime + MEM_WHEEL - (size_t)pSprime % MEM_WHEEL);
	Stock* pStock = StockCache + BucketInfo.StockSize;
	for (uint i = 0; i < MEM_BLOCK; i ++) {
		pStock->Sprime = pSprime + WHEEL_SIZE * i;
		pStock->Next = pStock + 1;
		pStock ++;
	}
	pStock[-1].Next = StockHead;
	StockHead = pStock - MEM_BLOCK;

	BucketInfo.StockSize += MEM_BLOCK;
	BucketInfo.CurStock  += MEM_BLOCK;
}

#define PI(x, r) x / log((double)x) * (1 + r / log((double)x))
static int setBucketInfo(const uint sieve_size, const uint sqrtp, const uint64 range)
{
	//assert (range >> 32 < sieve_size);
	BucketInfo.MaxBucket = range / sieve_size + 1;//overflow for big range
	BucketInfo.LoopSize = sqrtp / (sieve_size / MAX_WHEEL_GAP) + 2; //add the first, last bucket
	BucketInfo.Log2Size = ilog(sieve_size / WHEEL30, 2);
	BucketInfo.ModuloSize = (1 << BucketInfo.Log2Size) - 1;
	//assert (BucketInfo.Log2Size <= (32 - SI_BIT) && SI_BIT >= SP_BIT);
	return 0;
}

static void setWheelSmall(const uint64 start, const uint maxp)
{
	const uint l1_size = Threshold.L1Size * WHEEL30;
	uint64 sstart = start - start % WHEEL210;
	uint p = 0, j = Threshold.L1Index;

	for (; p < maxp; p = SmallPrime[++j].Prime) {
		SmallPrime[j].Si = 0;
	}

	for (p = SmallPrime[FIRST_INDEX].Prime, j = FIRST_INDEX; p < maxp; p = SmallPrime[++j].Prime) {
		const uint64 p2 = (uint64)p * p;
		if (p2 > sstart && p2 > sstart + l1_size) //overflow
			sstart += (p2 - sstart) / l1_size * l1_size;

		//assert(p2 < sstart);
		uint sieve_index = p - (uint)(sstart % p);
		if (p2 > sstart)
			sieve_index = (uint)(p2 - sstart);

		const WheelFirst wf = WheelFirst30[sieve_index % WHEEL30][WheelInit30[p % WHEEL30].WheelIndex];
		const uint multiples = (NEXT_MULTIPLE >> (wf.Multiple * 3)) | (NEXT_MULTIPLE << (24 - wf.Multiple * 3));
		sieve_index += wf.Correct * p;

		SmallPrime[j].Multiple = multiples;
#ifdef SM0
		SmallPrime[j].Si = sieve_index;
#else
		SmallPrime[j].Si = 0 - sieve_index;
#endif
	}
}

static int segmentedSieve2(uchar bitarray[], const uint start, const uint sieve_size, bool bcopy);
static void setWheelMedium(uchar* bitarray, const uint sieve_size, const uint medium, const uint64 start)
{
	uint j = Threshold.L1Index, l1_maxp = Threshold.L1Maxp - Threshold.L1Maxp % WHEEL210;
	Threshold.L2Index = Threshold.L3Index = 0;
	Threshold.L3Maxp = sieve_size / 6 + 1;

	uint msize = MIN(sieve_size, Threshold.L2Size * WHEEL30);
	const int bytes = segmentedSieve2(bitarray, l1_maxp, medium - l1_maxp, false);

	stype mask = 0, *sbitarray = (stype*)bitarray;
	uint noffset = l1_maxp - sizeof(mask) * WHEEL30;
	const uint pn = 1 + bytes / sizeof(mask);

	for (uint i = 0; i < pn; ) {
		if (mask == 0) {
			mask = ~sbitarray[i ++];
			noffset += sizeof(mask) * WHEEL30;
			continue;
		}

		const uint p = noffset + PRIME_OFFSET(mask); mask &= mask - 1;
		if (Threshold.L2Index == 0 && p >= Threshold.L2Maxp) {
			Threshold.L2Maxp = p;
			Threshold.L2Index = j;
			msize = sieve_size;
		}
		if (Threshold.L3Index == 0 && p > Threshold.L3Maxp) {
			Threshold.L3Maxp = p;
			Threshold.L3Index = j;
		}

		uint64 offset = start;
		const uint64 p2 = (uint64)p * p;
		if (p2 - msize >= offset) //overflow
			offset += (p2 - offset) / msize * msize;

		//assert (p2 < offset + msize);
		uint sieve_index = p - (uint)(offset % p);
		if (p2 > offset)
			sieve_index = (uint)(p2 - offset);

		const int pi = WHEEL_INIT[p % WHEEL].WheelIndex;
		const WheelFirst wf = WHEEL_FIRST[(sieve_index + (uint)(offset % WHEEL)) % WHEEL][pi];
		sieve_index += wf.Correct * p;
		MediumPrime[j + 0].Sp = (p / WHEEL << SI_BIT) + pi;
		MediumPrime[j ++ ].Si = (sieve_index / WHEEL30 << SI_BIT) + wf.WheelIndex;
	}
	MediumPrime[j + 0].Sp = MediumPrime[j + 1].Sp = MediumPrime[j + 2].Sp = UINT_MAX;
}

static void pushBucket(const uint offset, const uint sp, const uint wi)
{
	const uint next_bucket = offset >> BucketInfo.Log2Size;
#ifndef B_R
	if (next_bucket >= BucketInfo.MaxBucket)
		return;
#endif

	Bucket_* pbucket = Bucket + next_bucket;
	SievePrime* pSprime = pbucket->Sprime ++;
	if ((size_t)(pSprime) % MEM_WHEEL == 0) {
		BucketInfo.CurStock --;
//		if (BucketInfo.CurStock -- == 0) allocWheelBlock(1);
		//pop from free stock list
		Stock* psfree = StockHead; StockHead = psfree->Next;
		psfree->Next = pbucket->Head;
		pbucket->Head = psfree;
		pbucket->Sprime = (pSprime = psfree->Sprime) + 1;
	}

	pSprime->Sp = sp;
	pSprime->Si = (offset & BucketInfo.ModuloSize) << SI_BIT | wi;
}

static void setWheelBig(uchar* bitarray, uint medium, uint sqrtp, const uint64 start, const uint64 range)
{
	uint nextp = 0; uint64 remp = 0;
	if (sqrtp < UINT_MAX) sqrtp ++; //watch overflow if sqrtp = 2^32 - 1

	const uint irange = (uint)((range >> 32) > WHEEL30 ? UINT_MAX : range / WHEEL30);
	for (uint l2size = L2_DCACHE_SIZE * WHEEL30 << 10; medium < sqrtp; medium += l2size) {
		if (l2size > sqrtp - medium)
			l2size = sqrtp - medium;

		const int bytes = segmentedSieve2(bitarray, medium, l2size, nextp > 0);
		stype mask = 0, *sbitarray = (stype*)bitarray; //little endian
		uint offset = medium - sizeof(mask) * WHEEL30;
		const uint pn = 1 + bytes / sizeof(mask);

		for (uint j = 0; j < pn; ) {
			if (mask == 0) {
				mask = ~sbitarray[j ++];
				offset += sizeof(mask) * WHEEL30;
				continue;
			}

			const uint p = offset + PRIME_OFFSET(mask); mask &= mask - 1;
#if ASM_X86 && FDIV == 1
			if (p > nextp) { //ugly & difficult to understand but efficient
				remp = start / (nextp = p + (uint64)p * p / (uint)(start >> 32));
				if (p > nextp) //overflow
					remp = start >> 32, nextp = UINT_MAX;
			}
			uint sieve_index = p - fastDiv(start - remp * p, p);
#else
			uint sieve_index = p - (uint)(start % p);
#endif

#if 0
			if (sieve_index > range) continue;
#endif
			const uint sp = (p / WHEEL210 << SP_BIT) + WheelInit210[p % WHEEL210].WheelIndex;
			const uint modulo_210 = sieve_index % WHEEL210;
			const WheelFirst* wf = &WheelFirst210[modulo_210][sp % (1 << SP_BIT)];
#if X86_64
			sieve_index = (sieve_index + wf->Correct * (uint64)p) / WHEEL30;
#else
			sieve_index = (sieve_index / WHEEL210 + (sp >> SP_BIT) * wf->Correct) * (WHEEL210 / WHEEL30);
			sieve_index += (wf->Correct * (p % WHEEL210) + modulo_210) / WHEEL30;
#endif

#ifndef B_R
			if (sieve_index > irange) continue;
#endif
	//		pushBucket(sieve_index, sp, wf.WheelIndex);

			Bucket_* pbucket = Bucket + (sieve_index >> BucketInfo.Log2Size);
			SievePrime* pSprime = pbucket->Sprime ++;
			if ((size_t)(pSprime) % MEM_WHEEL == 0) {
				if (BucketInfo.CurStock -- == 0) allocWheelBlock(1);
				//pop from free stock list
				Stock* psfree = StockHead; StockHead = psfree->Next;
				psfree->Next = pbucket->Head;
				pbucket->Head = psfree;
				pSprime = psfree->Sprime;
				pbucket->Sprime = pSprime + 1;
			}
			pSprime->Sp = sp;
			pSprime->Si = (sieve_index & BucketInfo.ModuloSize) << SI_BIT | wf->WheelIndex;
		}
	}
}

static int crossSmall0(uchar bitarray[], const uchar* pend, const uint p, uint offset, uint multiples)
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
		pmax += p; //05143627
		OR_ADD(0); OR_ADD(5); OR_ADD(1); OR_ADD(4); OR_ADD(3); OR_ADD(6); OR_ADD(2); OR_ADD(7);
	}
	return (int)(pmax - ppbeg[wi]) * WHEEL30;
}

static int crossSmall1(uchar bitarray[], const uchar* pend, const uint p, uint offset, uint multiples)
{
	int mini = 1 << 31;
	for (int i = 0; i < 4; i ++) {
		uchar* ps1 = bitarray + offset / WHEEL30;
		const uchar masks1 = WheelInit30[offset % WHEEL30].MaskBit;
		offset += (multiples % 8) * p; multiples /= 8;

		uchar* ps2 = bitarray + offset / WHEEL30;
		const uchar masks2 = WheelInit30[offset % WHEEL30].MaskBit;
		offset += (multiples % 8) * p; multiples /= 8;

		while (ps2 <= pend) {
			*ps2 |= masks2, ps2 += p;
			*ps1 |= masks1, ps1 += p;
		}
		if (ps1 <= pend)
			*ps1 |= masks1, ps1 += p;

		int next;
		if (ps1 > ps2)
			next = (ps2 - pend) * WHEEL30 + PRIME_OFFSET(masks2);
		else
			next = (ps1 - pend) * WHEEL30 + PRIME_OFFSET(masks1);
		if (next < mini)
			mini = next;
	}

	return mini;
}

static int crossSmall3(uchar* ps, const uchar* pend, const uint p)
{
	#define POSET(s1,o1,b1,s2,o2,b2,s3,o3,b3,s4,o4,b4,s5,o5,b5,s6,o6,b6,s7,o7,b7)\
	while (ps <= pend) {\
		ps[o * 00 + 00] |= 1 <<0x0, ps[o * s1 + o1] |= BIT##b1;\
		ps[o * s2 + o2] |= BIT##b2, ps[o * s3 + o3] |= BIT##b3;\
		ps[o * s4 + o4] |= BIT##b4, ps[o * s5 + o5] |= BIT##b5;\
		ps[o * s6 + o6] |= BIT##b6, ps[o * s7 + o7] |= BIT##b7;\
		ps += p;\
	}

	const uchar* pbeg = ps;
	const uint o = p / WHEEL30;
	switch (WheelInit30[p % WHEEL30].WheelIndex) //37216405
	{
		case 3: POSET(4,1,6, 6,2,5,  10,4,2,  12,5,1,  16,6,7,  22,9,4,  24,10,3) break;
		case 7: POSET(2,1,7, 8,7,6,  12,11,5, 14,13,4, 18,17,3, 20,19,2, 24,23,1) break;
		case 2: POSET(2,0,6, 6,2,1,  8,2,7,   12,4,3,  18,6,5,  20,7,2,  26,9, 4) break;
		case 1: POSET(4,0,7, 6,1,3,  10,2,2,  16,3,6,  18,4,1,  24,5,5,  28,6, 4) break;
		case 6: POSET(2,1,4, 6,4,5,  12,9,1,  14,10,6, 20,15,2, 24,18,3, 26,19,7) break;
		case 4: POSET(6,3,3, 8,4,4,  14,7,7,  18,10,1, 20,11,2, 24,13,5, 26,14,6) break;
		case 0: POSET(6,0,1, 10,0,2, 12,0,3,  16,0,4,  18,0,5,  22,0,6,  28,0, 7) break;
		case 5: POSET(4,2,4, 10,6,2, 12,7,5,  18,11,3, 22,13,7, 24,15,1, 28,17,6) break;
	}

	return (int)(ps - pbeg) * WHEEL30;
}

#define MEDIUM_SET(n) \
	we##n = wd##n[(int)we##n.WheelIndex]; \
	bitarray[(int)offset##n] |= we##n.MaskBit; \
	offset##n += we##n.Correct + (we##n.Multiple) * wi##n

//cross out one medium prime from array
static int crossMedium1(uchar bitarray[], const uint sieve_byte, SievePrime* pSprime)
{
	uint si = pSprime->Si;
	int offset = (si >> SI_BIT) - sieve_byte;
	if (offset >= 0) {
		pSprime->Si -= sieve_byte << SI_BIT;
		return 1;
	}

	const uint wi = pSprime->Sp >> SI_BIT;
	WheelElement* wd = WHEEL_DATA[pSprime->Sp % (1 << SI_BIT)];
	WheelElement we; we.WheelIndex = si % (1 << SI_BIT);

	do {
		MEDIUM_SET();
	} while (offset < 0);

	pSprime->Si = (offset << SI_BIT) | we.WheelIndex;
	return 1;
}

//cross out 2 medium prime from array
static int crossMedium2(uchar bitarray[], const uint sieve_byte, SievePrime* pSprime)
{
	uint si1 = pSprime[0].Si, sp1 = pSprime[0].Sp;
	uint wi1 = sp1 >> SI_BIT, offset1 = (si1 >> SI_BIT) - sieve_byte;
	WheelElement* wd1 = WHEEL_DATA[sp1 % (1 << SI_BIT)];
	WheelElement we1; we1.WheelIndex = si1 % (1 << SI_BIT);

	uint si2 = pSprime[1].Si, sp2 = pSprime[1].Sp;
	uint wi2 = sp2 >> SI_BIT, offset2 = (si2 >> SI_BIT) - sieve_byte;
	WheelElement* wd2 = WHEEL_DATA[sp2 % (1 << SI_BIT)];
	WheelElement we2; we2.WheelIndex = si2 % (1 << SI_BIT);

	while ((int)offset1 < 0) {
		MEDIUM_SET(1);
		if ((int)offset2 >= 0) break;
		MEDIUM_SET(2);
	}

	while ((int)offset1 < 0) { MEDIUM_SET(1); } pSprime[0].Si = offset1 << SI_BIT | we1.WheelIndex;
	while ((int)offset2 < 0) { MEDIUM_SET(2); } pSprime[1].Si = offset2 << SI_BIT | we2.WheelIndex;

	return 2;
}

//more efficient for medium if prime < sieve_size / 30
static SievePrime* crossMediumW30(uchar bitarray[], const uint sieve_byte, const uint minsp, SievePrime* pSprime)
{
	for (uint sp = pSprime->Sp; sp < minsp; sp = pSprime->Sp) {
		const uint wi = sp >> SI_BIT, pi = sp % (1 << SI_BIT);
		const uint p = wi * WHEEL30 + Pattern30[pi];
		WheelElement* wd = WheelData30[pi];
		WheelElement* we = wd + pSprime->Si % (1 << SI_BIT);

		uchar* ps = bitarray + (pSprime->Si >> SI_BIT);
		uchar* pend = bitarray + sieve_byte;
		uchar* pmin = pend + p;

		int nwi = 0, lw1 = -1;
		int cw0, cw1;
		for (int i = 0; i < 8; i += 2) {
			uchar* ps1 = ps;
			const uchar masks1 = we->MaskBit;
			ps += we->Correct + we->Multiple * wi;
			cw0 = we->WheelIndex;
			we = wd + we->WheelIndex;

			cw1 = we->WheelIndex;
			uchar* ps2 = ps;
			const uchar masks2 = we->MaskBit;
			ps += we->Correct + we->Multiple * wi;
			we = wd + we->WheelIndex;
			while (ps2 < pend) {
				*ps2 |= masks2, ps2 += p;
				*ps1 |= masks1, ps1 += p;
			}
			if (ps1 < pend)
				*ps1 |= masks1, ps1 += p;

			if (ps2 < pmin)
				pmin = ps2, nwi = cw0;
			if (ps1 < pmin)
				pmin = ps1, nwi = lw1;
			lw1 = cw1;
		}

		if (nwi < 0) nwi = cw1;
		pSprime ++->Si = ((uint)(pmin - pend) << SI_BIT) | nwi;
	}

	return pSprime;
}

//cross out big prime's multiples from bucket
static void crossBig1(uchar bitarray[], const uint sieve_size, const SievePrime* pSprime)
{
	uint offset = pSprime->Si, sp = pSprime->Sp;
	const WheelElement* wd = WheelData210[sp % (1 << SP_BIT)];
	const WheelElement* we = wd + offset % (1 << SI_BIT);
	bitarray[offset >>= SI_BIT] |= we->MaskBit;
	offset += we->Correct + we->Multiple * (sp >> SP_BIT);

#if ERAT_BIG > 2
	if (offset < sieve_size) {
		we = wd + we->WheelIndex;
		bitarray[offset] |= we->MaskBit;
		offset += we->Correct + we->Multiple * (sp >> SP_BIT);
	}
#endif
	pushBucket(offset, sp, we->WheelIndex);
}

//cross out 2 big prime's multiples from bucket, 15% improvement
static int crossBig2(uchar bitarray[], const uint sieve_size, const SievePrime* pSprime)
{
	uint offset1 = pSprime[0].Si, sp1 = pSprime[0].Sp;
	uint offset2 = pSprime[1].Si, sp2 = pSprime[1].Sp;
	const WheelElement* wd1 = WheelData210[sp1 % (1 << SP_BIT)];
	const WheelElement* wd2 = WheelData210[sp2 % (1 << SP_BIT)];

	const WheelElement* we1 = wd1 + offset1 % (1 << SI_BIT);
	const WheelElement* we2 = wd2 + offset2 % (1 << SI_BIT);
	bitarray[offset1 >>= SI_BIT] |= we1->MaskBit;
	offset1 += we1->Correct + we1->Multiple * (sp1 >> SP_BIT);

	bitarray[offset2 >>= SI_BIT] |= we2->MaskBit;
	offset2 += we2->Correct + we2->Multiple * (sp2 >> SP_BIT);

#if ERAT_BIG > 2
	if (offset1 < sieve_size) {
		we1 = wd1 + we1->WheelIndex;
		bitarray[offset1] |= we1->MaskBit;
		offset1 += we1->Correct + we1->Multiple * (sp1 >> SP_BIT);
	}
	if (offset2 < sieve_size) {
		we2 = wd2 + we2->WheelIndex;
		bitarray[offset2] |= we2->MaskBit;
		offset2 += we2->Correct + we2->Multiple * (sp2 >> SP_BIT);
	}
#endif

	pushBucket(offset1, sp1, we1->WheelIndex);
	pushBucket(offset2, sp2, we2->WheelIndex);

	return 2;
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

	//wheel pattern < WHEEL30 is prime except 1
	if (start == 0)
		bitarray[0] = BIT0;

	//pack the last byte with bit 1
	if (bits % 8 != 0)
		bitarray[bits / 8] |= ~((1 << (bits & 7)) - 1);
}

//sieve prime multiples in [start, start + sieve_size) with small algorithm
static void eratSieveSmall(uchar bitarray[], const uint64 start, const uint sieve_size, uint maxp)
{
	if (start + sieve_size < ((uint64)maxp) * maxp)
		maxp = isqrt(start + sieve_size) + 1;

	const uchar* pend = bitarray + sieve_size / WHEEL30;
	for (uint p = SmallPrime[FIRST_INDEX].Prime, j = FIRST_INDEX; p < maxp; p = SmallPrime[++j].Prime) {
		uint offset = SmallPrime[j].Si;
#ifdef SM0
		offset += crossSmall0(bitarray, pend, p, offset, SmallPrime[j].Multiple) + 30 * p - sieve_size;
#else
		if ((int)offset < 0) {
			offset = 0 - offset;
			for (int i = 0; i < 24; i += 3) {
				const uchar mask = WheelInit30[offset % WHEEL30].MaskBit;
				if (mask == BIT0)
					break;

				bitarray[offset / WHEEL30] |= mask;
				offset += (SmallPrime[j].Multiple >> i) % 8 * p;
			}
		}
		offset += crossSmall3(bitarray + offset / WHEEL30, pend, p) - sieve_size;
#endif
		SmallPrime[j].Si = offset;
	}
}

//sieve prime multiples in [start, start + sieve_size) with medium algorithm
static void eratSieveMedium1(uchar bitarray[], const uint64 start, const uint sieve_size, const uint spi, uint maxp)
{
	if (start + sieve_size < (uint64)maxp * maxp)
		maxp = isqrt(start + sieve_size) + 1;

	const uint sieve_byte = sieve_size / WHEEL30 + (sieve_size % WHEEL30 ? 1 : 0);
	SievePrime* pSprime = MediumPrime + spi;

#ifndef W210
	const uint minsp = MIN(maxp, (sieve_byte * 80 / 100)) / 30 << SI_BIT;
	pSprime = crossMediumW30(bitarray, sieve_byte, minsp, pSprime);
#endif

	const uint maxsp = (maxp / WHEEL << SI_BIT) + WHEEL_INIT[maxp % WHEEL].WheelIndex;
	bitarray += sieve_byte;

	while (pSprime[1].Sp < maxsp)
		pSprime += crossMedium2(bitarray, sieve_byte, pSprime);

	if (pSprime[0].Sp < maxsp)
		crossMedium1(bitarray, sieve_byte, pSprime + 0);
}

//sieve prime multiples in [start, start + sieve_size) with medium algorithm
static void eratSieveMedium2(uchar bitarray[], const uint64 start, const uint sieve_size, const uint spi, uint maxp)
{
	if (start + sieve_size < (uint64)maxp * maxp)
		maxp = isqrt(start + sieve_size) + 1;

	const uint sieve_byte = sieve_size / WHEEL30 + (sieve_size % WHEEL30 ? 1 : 0);
	SievePrime* pSprime = MediumPrime + spi;

	const uint maxsp = (maxp / WHEEL << SI_BIT) + WHEEL_INIT[maxp % WHEEL].WheelIndex;
	bitarray += sieve_byte;

	while (pSprime[0].Sp < maxsp) {
#if 0
		pSprime += crossMedium1(bitarray, sieve_byte, pSprime);
#else
		uint si = pSprime->Si;
		int offset = (si >> SI_BIT) - sieve_byte;
		if (offset >= 0) {
			pSprime->Si -= sieve_byte << SI_BIT;
			pSprime ++;
			continue;
		}

		const uint wi = pSprime->Sp >> SI_BIT;
		WheelElement* wd = WHEEL_DATA[pSprime->Sp % (1 << SI_BIT)];
		WheelElement we; we.WheelIndex = si % (1 << SI_BIT);

		MEDIUM_SET();
		if (offset < 0)
		{
			MEDIUM_SET();
		}

		pSprime->Si = (offset << SI_BIT) | we.WheelIndex;
		pSprime ++;
#endif
	}
}

//This implementation uses a sieve array with WHEEL210 numbers per byte and
//a modulo wheel that skips multiples of 2, 3, 5 and 7.
static void eratSieveBig(uchar bitarray[], const uint sieve_size)
{
	uint loops = (size_t)Bucket[0].Sprime % MEM_WHEEL / sizeof(SievePrime);
	if (loops % 2)
		crossBig1(bitarray, sieve_size, Bucket[0].Sprime - 1), loops --;
	else if (loops == 0)
		loops = WHEEL_SIZE;

	for (Stock* pStock = Bucket[0].Head; pStock != NULL; loops = WHEEL_SIZE) {
		//push into free block list
		Stock* pnext = pStock->Next; pStock->Next = StockHead; StockHead = pStock;
		SievePrime* pSprime = StockHead->Sprime;
		while (loops) {
			crossBig2(bitarray, sieve_size, pSprime);
			pSprime += 2;
			loops -= 2;
		}
		BucketInfo.CurStock ++;
		pStock = pnext;
	}

	BucketInfo.MaxBucket --;
	memmove(Bucket, Bucket + 1, BucketInfo.LoopSize * sizeof(Bucket[0]));
}

static int segmentProcessed(uchar bitarray[], const uint64 start, uint bytes, PrimeCall* pcall)
{
	//pack last qword
	if (bytes % sizeof(uint64)) {
		*(uint64*)(bitarray + bytes) = ~0;
		bytes += sizeof(uint64) - bytes % sizeof(uint64);
	}

	if (pcall && pcall->func) {
		stype mask = ~(*(stype*)bitarray); //*(stype*) for big edian ?
		const int size = (bytes - 1) / sizeof(mask);
		for (int bi = 0; bi <= size; ) {
			if (mask == 0) {
				mask = ~((stype*)bitarray)[++bi];
				continue;
			}

			const uint64 p = start + bi * WHEEL30 * sizeof(mask) + PRIME_OFFSET(mask); mask &= mask - 1;
			pcall->primes += 1;
			pcall->func(pcall->data, p);//TODO
		}
	}

	return countBit0sArray((uint64*)bitarray, bytes);
}

static void eratSieveL1(uchar bitarray[], const uint64 start, const uint sieve_size)
{
	const uint l1_maxp = Threshold.L1Maxp;
	for (uint offset = 0, l1size = Threshold.L1Size * WHEEL30; offset < sieve_size; offset += l1size) {
		if (l1size + offset > sieve_size)
			l1size = sieve_size - offset;
		eratSieveSmall(bitarray + offset / WHEEL30, start + offset, l1size, l1_maxp);
	}
}

//sieve prime multiples in [start, start + sieve_size)
static int segmentedSieve(uchar bitarray[], const uint64 start, const uint sieve_size)
{
	const uint sqrtp = isqrt(start + sieve_size);
	const uint l1_maxp = Threshold.L1Maxp;
	const uint medium2 = MIN(sqrtp, Threshold.Medium) + 1;
	const bool bmsieve = sqrtp >= l1_maxp;

	preSieve(bitarray, start, sieve_size);

	//copy last L1 seg to the first seg
	const uint copy_from = Config.SieveSize;
	for (uint i = 0; i < l1_maxp; i += sizeof(uint64)) {
		uint64* pc = (uint64*)(bitarray + i + copy_from);
		*(uint64*)(bitarray + i) |= *pc; *pc = 0;
	}

	for (uint offset = 0, l2size = Threshold.L2Size * WHEEL30; offset < sieve_size; offset += l2size) {
		if (l2size + offset > sieve_size)
			l2size = sieve_size - offset;
		eratSieveL1(bitarray + offset / WHEEL30, start + offset, l2size);
		if (bmsieve)
			eratSieveMedium1(bitarray + offset / WHEEL30, start + offset, l2size, Threshold.L1Index, Threshold.L2Maxp);
	}

#ifdef FM
	const uint medium1 = MIN(sqrtp, Threshold.L3Maxp) + 1;
	if (medium1 > Threshold.L2Maxp)
		eratSieveMedium1(bitarray, start, sieve_size, Threshold.L2Index, medium1);
	if (medium2 > Threshold.L3Maxp)
		eratSieveMedium2(bitarray, start, sieve_size, Threshold.L3Index + 1, medium2);
#else
	if (medium2 > Threshold.L2Maxp)
		eratSieveMedium1(bitarray, start, sieve_size, Threshold.L2Index, medium2);
#endif

	if (start >= Threshold.BucketStart)
		eratSieveBig(bitarray, Config.SieveSize);

	return sieve_size / WHEEL30 + (7 + WheelInit30[sieve_size % WHEEL30].WheelIndex) / 8;
}

static int segmentedSieve2(uchar bitarray[], const uint start, const uint sieve_size, bool bcopy)
{
	const uint sqrtp = isqrt(start + sieve_size);
	const uint l1_maxp = Threshold.L1Maxp;
	preSieve(bitarray, start, sieve_size);

#if FSL1
	const uint copy_from = L2_DCACHE_SIZE << 10;
	if (sieve_size <= copy_from * WHEEL30 && bcopy) {
		for (uint i = 0; i < l1_maxp; i += sizeof(uint64)) {
			uint64* pc = (uint64*)(bitarray + i + copy_from);
			*(uint64*)(bitarray + i) |= *pc; *pc = 0;
		}
	}
	eratSieveL1(bitarray, start, sieve_size);
#else
	for (uint offset = 0, l1size = Threshold.L1Size * WHEEL30; offset < sieve_size; offset += l1size) {
		if (l1size + offset > sieve_size)
			l1size = sieve_size - offset;

		//corss out sieve size = L1_SIZE
		const uint64 newstart = start + offset;
		const uint lsqrtp = isqrt(newstart + l1size) + 1;
		uchar* pstart = bitarray + offset / WHEEL30;
		const uchar* pend = bitarray + (offset + l1size) / WHEEL30;
		const uint maxp = MIN(lsqrtp, l1_maxp);

		for (uint p = SmallPrime[FIRST_INDEX].Prime, j = FIRST_INDEX; p < maxp; p = SmallPrime[++j].Prime) {
			uint offset1 = SmallPrime[j].Si;
	#ifndef SM0
			if ((int)offset1 < 0)
				offset1 = 0 - offset1;
	#endif
			SmallPrime[j].Si = offset1 + crossSmall0(pstart, pend, p, offset1, SmallPrime[j].Multiple) + 30 * p - l1size;
		}
	}
#endif

	//	assert(l1_maxp > start / l1_maxp);
	for (uint j = Threshold.L1Index, p = l1_maxp; p <= sqrtp; p = SmallPrime[++j].Prime) {
		uint offset = SmallPrime[j].Si;
		if ((int)offset <= 0) {
			offset = p - start % p;
			if (p * p > start)
				offset = p * p - start;
		}

		const WheelFirst wf = WheelFirst30[offset % WHEEL30][WheelInit30[p % WHEEL30].WheelIndex];
		const uint multiples = (NEXT_MULTIPLE >> (wf.Multiple * 3)) | (NEXT_MULTIPLE << (24 - wf.Multiple * 3));
		SmallPrime[j].Si = crossSmall1(bitarray, bitarray + sieve_size / WHEEL30, p, offset + wf.Correct * p, multiples);
	}

	uint bytes = sieve_size / WHEEL30 + (7 + WheelInit30[sieve_size % WHEEL30].WheelIndex) / 8;
	if (bytes % sizeof(uint64)) {
		*(uint64*)(bitarray + bytes) = ~0;
		bytes += sizeof(uint64) - bytes % sizeof(uint64);
	}

	return bytes;
}

static void setL1Index()
{
	if (Threshold.L1Maxp > 65521)
		Threshold.L1Maxp = 65521; //the bigest prime % WHEEL210 = 1 and < 2^16
	for (uint p = 0, j = 300; ; p = SmallPrime[++j].Prime) {
		if (p >= Threshold.L1Maxp && p % WHEEL210 == 1) {
			Threshold.L1Index = j;
			Threshold.L1Maxp = p;
			break;
		}
	}
}

void setCacheSize(int level, uint cache)
{
	//	cache = 1 << ilog(cache, 2);
	if (level == 1 && cache >= 16 && cache <= (Threshold.L2Size >> 10)) {
		Threshold.L1Size = cache << 10;
		Threshold.L1Maxp = Threshold.L1Size / Config.L1Segs;
		Threshold.L2Size = Threshold.L2Size / Threshold.L1Size * Threshold.L1Size;
		Threshold.L2Maxp = Threshold.L2Size / Config.L2Segs;
		setL1Index();
	} else if (level == 2 && cache >= (Threshold.L1Size >> 10) && cache <= MAX_SEGMENT) {
		Threshold.L2Size = (cache << 10) / Threshold.L1Size * Threshold.L1Size;
		Threshold.L2Maxp = Threshold.L2Size / Config.L2Segs;
	}
}

void setCacheSegs(uint level, uint segs)
{
	if (segs == 0 || segs > 10)
		return;

	if (level == 1) {
		Config.L1Segs = segs;
		Threshold.L1Maxp = Threshold.L1Size / segs;
		setL1Index();
	} else if (level == 2) {
		Config.L2Segs = segs;
		Threshold.L2Maxp = Threshold.L2Size /segs;
	} else if (level == 3 && segs <= 6 && ERAT_BIG > 2) {
		Config.Msegs = segs;
	}
}

//sieve_size : 32k - 4M
uint setSieveSize(uint sieve_size)
{
	const uint l1_size = Threshold.L1Size >> 10, l2_Size = Threshold.L2Size >> 10;
	if (sieve_size <= MAX_SEGMENT && sieve_size > l2_Size) {
		sieve_size = (sieve_size / l2_Size) * l2_Size << 10;// 1 << (ilog(sieve_size, 2) + 10);
	} else if (sieve_size <= l2_Size && sieve_size >= l1_size) {
		sieve_size = (sieve_size / l1_size) * l1_size << 10;
	} else if (sieve_size <= (MAX_SEGMENT >> 10) && sieve_size > 0) {
		sieve_size = sieve_size << 20;
	} else {
		sieve_size = L2_DCACHE_SIZE << 10;
	}

	return Config.SieveSize = sieve_size;
}

static int checkSmall(const uint64 start, const uint64 end, PrimeCall* pcall)
{
	int primes = 0;
	for (int i = 0, p = 2; p < 7; p = SmallPrime[++i].Prime) {
		if (start <= p && p <= end) {
			primes ++;
			if (pcall && pcall->func) {
				pcall->primes += 1;
				pcall->func(pcall->data, p);
			}
		}
	}

	return primes;
}

static int64 pi(uchar* bitarray, uint64 start, uint64 end, PrimeCall* pcall)
{
	const int64 ts = getTime();
	uint64 primes = checkSmall(start, end, pcall);

	uint align210 = (uint)(start % WHEEL210);
	start -= align210;

	if (++end == 0) end --; //overflow if end = 2^64-1

	//pi(n) ~= n/log(n), complexty ~= n*log(log(n)), replaced by n*log(n) more accurate
#ifndef _DEBUG
	double lge = end * log((double)end), lgs = start * log(start + 10.0);
	double pie = PI(end, 1.2), pis = PI(start, 1.2);
#endif

	for (uint si = 0, sieve_size = Config.SieveSize * WHEEL30; start < end; start += sieve_size) {
		if (sieve_size > end - start)
			sieve_size = (uint)(end - start);

		const uint bytes = segmentedSieve(bitarray, start, sieve_size);
		if (align210 > 0) { //why brother ?
			memset(bitarray, ~0, align210 / WHEEL30);
			bitarray[align210 / WHEEL30] |= (1 << WheelInit30[align210 % WHEEL30].WheelIndex) - 1;
			align210 = 0;
		}

		primes += segmentProcessed(bitarray, start, bytes, pcall);

#ifndef _DEBUG
		if ((si ++ & Config.Progress) == 15) {
			const double cur = start + sieve_size;
			double lgc = cur * log(cur), pic = PI(cur, 1.2);
			double tratio = (lgc - lgs) / (lge - lgs) * 100;
			double pratio = (pic - pis) / (pie - pis);
			double timeuse = (getTime() - ts) / (10 * tratio);
			const uint64 picount = (int64)((int64)primes / pratio);
			if (timeuse < 3600)
				printf(">> %.2f%%, time ~= %.2f %s, primes ~= %llu\r", tratio, timeuse, "sec", picount);
			else
				printf(">> %.2f%%, time ~= %.3f %s, primes ~= %llu\r", tratio, timeuse / 60, "min", picount);
			fflush(stdout);
		}
#endif
	}

	return primes;
}

static void convertSci(uint64 n, char buff[40])
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
		if (r % powr == 0 && logr > 4)
			sprintf(buff, "%de%d+%de%d", (int)(n / pown), logn, (int)(r / powr), logr);
		else if (r % powr == 0)
			sprintf(buff, "%de%d+%d", (int)(n / pown), logn, (int)r);
	}
}

static void printResult(const uint64 start, const uint64 end, uint64 primes)
{
	char buff[64] = {0};
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
	uint medium = sieve_size / Config.Msegs + 1;

	uint64 offset = (uint64)medium * medium;
	if (offset > start && sqrtp > medium) {
		offset += sieve_size - offset % sieve_size + start % sieve_size;
		while (offset % WHEEL210 != start % WHEEL210) offset += sieve_size;
		medium = isqrt(offset - offset % WHEEL210 + sieve_size) + 1;
	} else
		offset = start;

	if (medium > sqrtp)
		medium = sqrtp;
	medium += WHEEL210 - medium % WHEEL210;

	Threshold.BucketStart = offset - offset % WHEEL210;
	Threshold.Medium = medium;

	return medium;
}

static uint adjustSieveSize(uint sieve_size, const uint sqrtp)
{
	if (sieve_size < Threshold.L2Size && sqrtp > 10000000)
		sieve_size = setSieveSize(SIEVE_SIZE);
	else if (sieve_size >= Config.Msegs << 21)
		sieve_size = setSieveSize(SIEVE_SIZE);

	if ((sieve_size & (sieve_size - 1)) != 0 && sqrtp > sieve_size * WHEEL30 / Config.Msegs)
		sieve_size = setSieveSize(1 << ilog(sieve_size >> 10, 2));
	return sieve_size * WHEEL30;
}

static void allocMedium(const uint medium)
{
	if (MediumPrime != NULL && (MediumPrime[0].Sp < medium || MediumPrime[0].Si > 2e6)) {
		free(MediumPrime); MediumPrime = NULL;
	}
	if (MediumPrime == NULL) {
		const uint pix = (uint)(PI(medium, 1.2));
		MediumPrime = (SievePrime*) malloc(sizeof(SievePrime) * (pix + PI_65536));
		MediumPrime[0].Sp = medium;
		MediumPrime[0].Si = pix + PI_65536;
	}
}

uint64 doSieve(const uint64 start, const uint64 end, PrimeCall* pcall)
{
	const int64 ts = getTime();
	const uint sqrtp = isqrt(end);
	const uint sieve_size = adjustSieveSize(Config.SieveSize, sqrtp);

	const uint medium = setBucketStart(start, sqrtp, sieve_size);
	allocMedium(medium);

	const uint l1_maxp = Threshold.L1Maxp;
	const uint max_cache = MAX(sieve_size / WHEEL30, L2_DCACHE_SIZE * 1024) + l1_maxp + sizeof(uint64);
	uchar* bitarray = (uchar*) malloc(max_cache);

	//init medium sieve
	const bool bmsieve = sqrtp >= l1_maxp;
	if (bmsieve) {
		setWheelSmall(l1_maxp, isqrt(medium) + 10);
		setWheelMedium(bitarray, sieve_size, medium, start - start % WHEEL210);
	}

	//init big sieve
	uint64 buckets = Threshold.BucketStart;
	if (sqrtp > medium) {
		setBucketInfo(sieve_size, sqrtp, end - buckets);
		setWheelSmall(medium, isqrt(sqrtp) + 1);
		memset(bitarray + (L2_DCACHE_SIZE << 10), 0, l1_maxp + sizeof(uint64));
		setWheelBig(bitarray, medium, sqrtp, buckets, end - buckets);
		if (BucketInfo.CurStock < BucketInfo.LoopSize)
			allocWheelBlock(1);
	} else
		Threshold.BucketStart = end + 1;

	//init small sieve
	{
		setWheelSmall(start, l1_maxp);
		memset(bitarray + sieve_size / WHEEL30, 0, max_cache - sieve_size / WHEEL30);
	}

	const int64 ti = getTime();
	const int64 primes = pi(bitarray, start, end, pcall);

	if (BucketInfo.StockSize > MAX_STOCK / 5) {
		for (uint i = 0; i < BucketInfo.PoolSize; i ++) free(WheelPool[i]);
		memset(&BucketInfo, 0, sizeof(BucketInfo));
	}

	if (Config.Flag & PRINT_RET) {
		const int64 ta = getTime();
		printResult(start, end, primes);
		if (Config.Flag & PRINT_TIME)
			printf(" (%.2f + %.2f = %.2f sec %u kb L%u%u%u)",
					(ti - ts) / 1000.0, (ta - ti) / 1000.0, (ta - ts) / 1000.0, Config.SieveSize >> 10, Config.L1Segs, Config.L2Segs, Config.Msegs);
		putchar('\n');
	}

#ifndef B_R
	assert(BucketInfo.StockSize == BucketInfo.CurStock);
#endif

	free(bitarray);
	return primes;
}

#if X86_64 || X86
static void cpuidInfo(int cpuinfo[4], int id)
{
#if _MSC_VER >= 1400
	__cpuid(cpuinfo, id);
#elif _MSC_VER
	__asm
	{
		mov eax, id
		cpuid
		mov edi, cpuinfo
		mov dword ptr [edi + 0], eax
		mov dword ptr [edi + 4], ebx
		mov dword ptr [edi + 8], ecx
		mov dword ptr [edi +12], edx
	}
#elif __GNUC__ || __TINYC__
	__asm
	(
		"cpuid\n"
		:"=a"(cpuinfo[0]),"=b"(cpuinfo[1]),"=c"(cpuinfo[2]),"=d"(cpuinfo[3])
		:"a"(id)
	);
#elif __clang_version__

#endif
}

static int getCpuInfo()
{
	char cpuName[257] = {-1};
	int  (*pTmp)[4] = (int(*)[4])cpuName;
	cpuidInfo(*pTmp ++, 0x80000002);
	cpuidInfo(*pTmp ++, 0x80000003);
	cpuidInfo(*pTmp ++, 0x80000004);

	for (int i = 0; cpuName[i]; i ++) {
		if (cpuName[i] != ' ' || cpuName[i + 1] != ' ')
			putchar(cpuName[i]);
	}

	int cpuinfo1[4]; cpuidInfo(cpuinfo1, 0x80000005); //work for amd cpu
	int cpuinfo2[4]; cpuidInfo(cpuinfo2, 0x80000006);
//	int cpuinfo3[4]; cpuidInfo(cpuinfo3, 0x2);

	if (cpuinfo1[2] >> 24 >= 16) {
		Threshold.L1Size = (cpuinfo1[2] >> 24) << 10;
		Threshold.L1Maxp = Threshold.L1Size / Config.L1Segs;
	} else if (cpuName[0] == 'I') { //intel cpu is too complexty to get l1size
		cpuinfo1[2] = (32 << 24);
		Threshold.L1Size = (cpuinfo1[2] >> 24) << 10;
		Threshold.L1Maxp = Threshold.L1Size / Config.L1Segs;
	}

	//cat /sys/devices/system/cpu/cpu0/cache/index0/size
	if (cpuinfo2[2] >> 16 >= 64) {
		Threshold.L2Size = (cpuinfo2[2] >> 16) << 10;
		Threshold.L2Maxp = Threshold.L2Size / Config.L2Segs;
	}

	printf("  Cpu Cache L1Dsize = %dk, L2Size = %dk,l3Size = %dk\n", cpuinfo1[2] >> 24, cpuinfo2[2] >> 16, cpuinfo2[3] >> 12);
//	printf(" %x %x %x %x", cpuinfo3[3], cpuinfo3[2], cpuinfo3[1], cpuinfo3[0]);
	return Threshold.L2Size;
}
#endif

void initPrime(int sieve_size)
{
	if (SmallPrime[20].Prime == 0) {
#if X86_64 || X86
#ifndef __clang_version__
		getCpuInfo();
#endif
#endif
		eratoSimple();
		initBitTable();
		initWheel30();
		initWheel210();
		setL1Index();
		setSieveSize(sieve_size);
	}
}

static void fixRangeTest(uint64 lowerBound, const int64 range, uint64 Ret)
{
	const int llog10 = ilog(lowerBound, 10), rlog10 = ilog(range, 10);
	const uint64 maxrange = ipow(2, 34);
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

//		setCacheSize(1, 32 * (rand() % 2 + 1));
		setCacheSize(2, 256 * (rand() % 2 + 1));
		setSieveSize(L2_DCACHE_SIZE << (rand() % 4));
		setCacheSegs(1, rand() % 5 + 1), setCacheSegs(2, rand() % 5 + 2), setCacheSegs(3, rand() % 5 + 2);

		primes += doSieve(lowerBound, end - 1, NULL);
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

static void startTest(int flag)
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
		155428406, // pi(10^12, 2^32)
		143482916, // pi(10^13, 2^32)
		133235063, // pi(10^14, 2^32)
		124350420, // pi(10^15, 2^32)
		116578809, // pi(10^16, 2^32)
		109726486, // pi(10^17, 2^32)
		103626726, // pi(10^18, 2^32)
		98169972,  // pi(10^19, 2^32)
		0,
	};

	int64 ts = getTime();
	Config.Flag &= ~PRINT_RET;
	uint64 primes = 0;
	Config.Progress = 0;
	for (int i = 1; i <= 10; i ++) {
		primes = doSieve(0, ipow(10, i), NULL);
		printf("pi(10^%2d) = %llu\n", i, primes);
	}

	srand((uint)time(0));
	for (int j = 12; primeCounts[j]; j ++) {
		uint64 start = ipow(10, j), end = start + ipow(2, 30);
		setCacheSegs(1, rand() % 6 + 2), setCacheSegs(2, rand() % 6 + 2); setCacheSegs(3, rand() % 6 + 2);
		primes = doSieve(start, end, NULL);
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
		{ipow(01, 12), pow12, pow9 * 37 + 607912018},
		{ipow(10, 17), pow11, 2554712095ul},
		{ipow(10, 12), pow11, 3612791400ul},
		{ipow(10, 13), pow11, 3340141707ul},
		{-1 - pow11,   pow11, 2254197466ul},
		{ipow(10, 14), pow12, pow9 * 31 + 16203073},
		{ipow(10, 15), pow12, pow9 * 28 + 952450479},
		{ipow(10, 16), pow12, pow9 * 27 + 143405794},
		{ipow(10, 18), pow12, pow9 * 24 + 127637783},
		{ipow(10, 19), pow12, pow9 * 22 + 857444126},
		{ipow(10, 19), pow12, pow9 * 22 + 857444126},
		{-1 - pow12,   pow12, pow9 * 22 + 542106206},
//		{-1 - pow12*10, pow12*10, pow9 * 225 + 420940155},
	};

	for (uint k = 0; k < sizeof(rangeData) / sizeof(rangeData[0]); k ++)
		fixRangeTest(rangeData[k][0], rangeData[k][1], rangeData[k][2]);

	Config.Flag |= PRINT_RET;
}

static void printInfo()
{
	const char* sepator =
		"------------------------------------------------------------------------------------------------------------";
	puts(sepator);
	puts("Fast implementation of the segmented sieve of Eratosthenes 2^64\n"
	"Copyright (C) by 2010-2018 Huang Yuanbing 22738078@qq.com/bailuzhou@163.com\n"
	"Compile: g++ -DSIEVE_SIZE=2048 -DFSL1 -DFDIV -march=native -funroll-loops -O3 -s -pipe PrimeNumber.cpp -o prime\n");

	char buff[500];
	char* info = buff;
#ifdef __clang_version__
	info += sprintf(info, "Compiled by clang %s", __clang_version__);
#elif _MSC_VER
	info += sprintf(info, "Compiled by vc++ %d", _MSC_VER);
#elif __GNUC__
	info += sprintf(info, "Compiled by gcc %s", __VERSION__);
#elif __TINYC__
	info += sprintf(info, "Compiled by tcc %s", __TINYC__);
#endif

#if __cplusplus
	info += sprintf(info, " c++ %d", (int)__cplusplus);
#endif

#if X86_64
	info += sprintf(info, " x86-64");
#elif X86
	info += sprintf(info, " x86");
#endif

	info += sprintf(info, " %s %s in %s\n", __TIME__, __DATE__, __FILE__);
	info += sprintf(info, "[MARCO] : MEM_WHEEL = %dM, WHEEL_SIZE = %dk, SIEVE_SIZE = %dk, WHEEL = %d\n",
			MEM_BLOCK * WHEEL_SIZE >> 17 , WHEEL_SIZE >> 7, SIEVE_SIZE, WHEEL);
	info += sprintf(info, "[CACHE] : L1Size = %u, L2Size = %u, SieveSize = %u, Bucket = %u, Block = %u\n",
			Threshold.L1Size >> 10, Threshold.L2Size >> 10, Config.SieveSize >> 10, BucketInfo.LoopSize, BucketInfo.StockSize);
	info += sprintf(info, "[ARGS ] : L1Segs/L2Segs/Mseg = (%u,%u,%u)\n",
			Config.L1Segs, Config.L2Segs, Config.Msegs);
	info += sprintf(info, "[ARGS ] : L1Maxp/L2Maxp/L3Maxp/Medium/Large/SieveSize = (%u,%u,%u,%u,%u,%u)",
			Threshold.L1Maxp, Threshold.L2Maxp, Threshold.L3Maxp, Threshold.Medium, isqrt(Threshold.BucketStart + 1), Config.SieveSize);
	*info = 0;
	puts(buff);
	puts(sepator);
}

static int pi2(uint start, uint end)
{
	uint primes = checkSmall(start, end, NULL);
	const uint sstart = start;
	uint align210 = start % WHEEL210;
	start -= align210;

	if (++end == 0) end --; //overflow if end = 2^64-1
	setWheelSmall(start, isqrt(end) + 1);
	const uint sieve_size = (L2_DCACHE_SIZE * 1024);
	const uint max_cache = sieve_size + sizeof(uint64) + Threshold.L1Maxp;
	uchar* bitarray = (uchar*) malloc(max_cache);
	memset(bitarray + sieve_size, 0, max_cache - sieve_size);

	int64 ts = getTime();
	for (uint l2size = sieve_size * WHEEL30, beg = start; beg < end; beg += l2size) {
		if (l2size > end - beg)
			l2size = end - beg;
		const int bytes = segmentedSieve2(bitarray, beg, l2size, primes > 4);
		if (align210 > 0) { //why brother ?
			memset(bitarray, ~0, align210 / WHEEL30);
			bitarray[align210 / WHEEL30] |= (1 << WheelInit30[align210 % WHEEL30].WheelIndex) - 1;
			align210 = 0;
		}
		primes += segmentProcessed(bitarray, beg, bytes, NULL);
	}

	if (Config.Flag & PRINT_RET) {
		const int64 ta = getTime();
		printResult(sstart, end - 1, primes);
		printf(" (%.2f sec, %u kb)\n", (ta - ts) / 1000.0, sieve_size >> 10);
	} else {
		puts("\r");
	}

	return primes;
}

//get the first digit number index
static int parseCmd(char params[][60])
{
	int cmdi = -1, cdata = 0;

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
//			case 'H': puts(Benchmark); if (n) puts(Help);   break;
			case 'S': setSieveSize(cdata);                  break;
			case 'C': setCacheSize(cdata % 10, cdata / 10); break;
			case 'L': setCacheSegs(cdata % 10, cdata / 10); break;
			case 'M': Config.Progress = (1 << cdata) - 1;   break;
			case 'D': Config.Flag ^= 1 << (n - 'A');        break;
			case 'I': printInfo();                          break;
			default : cmdi = i;                             break;
		}
	}

	return cmdi;
}

//split ccmd string to params array
static int splitCmd(const char* ccmd, char params[][60])
{
	int nwords = 0;

	for (int i = 0; i < 256; i ++) {
		while (isspace(*ccmd) || ',' == *ccmd)
			ccmd ++;

		char* pc = params[i];
		unsigned char c = *ccmd;
		if (c == 0 || c == ';')
			break;

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

static void dumpPrime(void* data, const uint64 prime)
{
	printf("%llu %llu\n", *((uint64*)data), prime);
}

static bool executeCmd(const char* cmd)
{
	while (cmd) {
		char params[10][60] = {0};
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
		if (end == 0)
			end = start, start = 0;
		else if (end < start)
			end += start;
		if (end < start)
			end = -1;

		if (cmdc == 'B') {
			if (isdigit(params[cmdi + 1][0])) {
				int powi = atoi(params[cmdi + 1]);
				uint64 range = powi > 12 ? ipow(2, powi) : ipow(10, powi);
				for (int j = 11; j < 20 && powi > 0; j ++) {
					uint64 start2 = ipow(10, j);
					doSieve(start2, start2 + range, NULL);
				}
			}
			if (isdigit(params[cmdi + 2][0])) {
				int powi = atoi(params[cmdi + 2]);
				uint64 range = powi > 12 ? ipow(2, powi) : ipow(10, powi);
				for (int i = 32; i < 64 && powi > 0; i ++) {
					uint64 start2 = ipow(2, i);
					doSieve(start2, start2 + range, NULL);
				}
			} else
				startTest(0);
		} else if (cmdc == 'P') {
#if PC
			PrimeCall pcall = {0, dumpPrime, 0};
			pcall.data = &pcall.primes;
			doSieve(start, end, &pcall);
#endif
		} else if (cmdi >= 0) {
			doSieve(start, end, NULL);
			if (end >> 32 == 0) pi2(start, end);
		}

		if (pcmd)
			cmd = pcmd + 1;
		else
			break;
	}

	return true;
}

int main(int argc, char* argv[])
{
	srand((uint)time(0));

	initPrime(SIEVE_SIZE);
	if (argc == 2)
		executeCmd(argv[1]);

#if RT
	else if (argc > 2)
	{
	Config.Progress = 0;
	Config.Flag ^= PRINT_RET;
//	Config.Flag &= ~PRINT_TIME;
#if 1
	for (int i = 30; i > 0; i --)
	for (int j = 1; j < 10; j ++) {
		setCacheSize(1, 32 * (rand() % 1 + 1)); setCacheSize(2, 256 * (rand() % 2 + 1));
		uint64 beg = ((uint64)(rand() * rand()) << i) + (uint64)rand() * rand() * rand();
		int sieve_size = setSieveSize(256 << (rand() % 5));
		uint64 range = (rand() % 36) * (uint64)sieve_size * WHEEL30 + rand();
		uint64 rm1 = doSieve(beg, beg + range, NULL);
		if (setSieveSize(L2_DCACHE_SIZE << (rand() % 5 + 0)) == sieve_size)
			setSieveSize(sieve_size*2);

		Config.Flag ^= PRINT_RET;
		setCacheSegs(3, rand() % 5 + 2); setCacheSegs(2, rand() % 6 + 2), setCacheSegs(1, rand() % 5 + 1);
		uint64 rm2 = doSieve(beg, beg + range, NULL);
		if (rm1 != rm2) { printInfo(); printf("%llu %llu\n", beg, range); system("pause"); }
		if ((i + j) % 10 == 0) printf("\r %2d progress ~= %d\n", i, j);
		Config.Flag ^= PRINT_RET;
	}
#endif

//	Config.Flag &= ~PRINT_TIME;
//	Config.Flag ^= PRINT_RET;

	for (int j = 1; j <= 10000; j ++) {
		char cmd[100] = {0};
#if 1
		setCacheSize(1, 32 * (rand() % 1 + 1)); setCacheSize(2, L2_DCACHE_SIZE * (rand() % 2 + 1));
		setCacheSegs(1, rand() % 6 + 1);
		uint64 beg = (uint64)(rand() * rand()) * (uint64)(rand() * rand()) % ipow(10, 12);
		setSieveSize(L2_DCACHE_SIZE * rand() % 8 + L2_DCACHE_SIZE);
		beg -= beg % 2;
		uint64 range = (rand() % 64 + 2) * Config.SieveSize * WHEEL30;
		uint64 rm1 = doSieve(beg, beg + range - 1, NULL);
//		Config.Flag ^= PRINT_RET;

		setCacheSegs(1, rand() % 6 + 2);
		setCacheSegs(2, rand() % 8 + 3);
		setSieveSize(L2_DCACHE_SIZE << ((rand() % 5) + 0));

		uint64 rm2 = doSieve(beg, beg + range, NULL);
		sprintf(cmd, "%llu != %llu [%d] --- %llu %llu s%u L%u1 L%u2 L%u3", rm1, rm2,
			j, beg, range, Config.SieveSize >> 10, Config.L1Segs, Config.L2Segs, Config.Msegs);
		if (rm1 != rm2)
			puts(cmd);
#endif
		for (int i = 2; i <= 6; i ++) {
			setSieveSize(L2_DCACHE_SIZE << ((rand() % 5) + 0));
			uint sieve_size = Config.SieveSize * WHEEL30;
			uint64 medium = sieve_size / i + 1;
			medium -= medium % WHEEL210;
			uint64 start = medium * medium - sieve_size * (rand() % 8 + 1) + rand();
			uint64 end = start + sieve_size * (rand() % 32 + 1);

			setCacheSegs(3, i); setCacheSegs(2, rand() % 6 + 2), setCacheSegs(1, rand() % 6 + 1);
			const uint64 r1 = doSieve(start, end, NULL);
			char cmd2[100] = {0};
			sprintf(cmd2, "  %llu %llu %llu s%u L%u1 L%u2 L%u3",
					r1, start, end, (sieve_size / WHEEL30) >> 10, Config.L1Segs, Config.L2Segs, Config.Msegs);

			setSieveSize(L2_DCACHE_SIZE << ((rand() % 5) + 0));
			setCacheSegs(3, rand() % 5 + 2); setCacheSegs(2, rand() % 6 + 2), setCacheSegs(1, rand() % 5 + 1);
			const uint64 r2 = doSieve(start, end, NULL);
			if (r1 != r2) {
				puts(cmd);
				puts(cmd2);
				printf("   %llu != %llu, %llu %llu s%u L%u1 L%u2 L%u3\n",
						r1, r2, start, end, Config.SieveSize >> 10, Config.L1Segs, Config.L2Segs, Config.Msegs);
			}
		}
		if (j % 10 == 0) printf("\rprogress ~= %d %d%%\n", j, 100*j/1000), fflush(stdout);
	}}
#endif

#ifdef B_R
	executeCmd("e16 e10;");
#else
	if (Threshold.L2Size == 512 << 10 || Config.SieveSize == 4096 << 10)
		executeCmd("L41 L42;");
	executeCmd("1e12 1e10; e14 e10 ;  e10+0");
	executeCmd("10^12 1e9; e16 e9; e18 e9*1; i 0-e9 0-1");
#endif

	while (true) {
		char ccmd[257];
		printf(">> ");
		if (!fgets(ccmd, 100, stdin) || !executeCmd(ccmd))
			break;
	}

	return 0;
}

//TODO
//1.multi-threading
//2.improve bucket algorithm for big range ex[1e14, 1e16]
//3.redesign by c++
// cl /O2 /Oi /GS /GL /D SIEVE_SIZE=4096 /D L2_DCACHE_SIZE=512 /D X86_64 PrimeNumber98.c -o p98_vs
