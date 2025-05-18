/*Fast Single-thread segmented sieve of Eratosthenes prime number n < 2^64 ***/
////http://ntheory.org/sieves/benchmarks.html
static const char* Benchmark =
"g++ -DFDIV -march=native -funroll-loops -O2 -s -pipe PrimeNumber.cpp -o prime\n"
"Windows 10 x64               i3-350M,i5-3470,i7-7500u,i7-6700,r7-1700,r7-5800H,i7-9700\n"
"Pi(0,    1e10) = 455052511    3.10   1.84    1.55     1.32    1.65    1.03     1.19\n"
"Pi(1e11, 1e10) = 394050419    4.35   2.50    2.02     1.78    2.00    1.31     1.50\n"
"Pi(1e12, 1e10) = 361840208    5.30   3.00    2.40     2.04    2.25    1.52     1.77\n"
"Pi(1e13, 1e10) = 334067230    6.52   3.50    2.85     2.39    2.67    1.71     2.06\n"
"Pi(1e14, 1e10) = 310208140    7.90   4.20    3.50     2.87    3.20    2.06     2.55\n"
"Pi(1e15, 1e10) = 289531946    9.90   5.10    4.22     3.49    3.91    2.50     3.14\n"
"Pi(1e16, 1e10) = 271425366    11.7   6.10    4.81     4.12    4.73    3.11     3.75\n"
"Pi(1e17, 1e10) = 255481287    13.9   7.09    5.55     4.84    5.63    3.72     4.45\n"
"Pi(1e18, 1e10) = 241272176    17.2   8.58    6.70     5.82    6.88    4.53     5.31\n"
"Pi(1e19, 1e10) = 228568014    24.0   11.6    9.57     8.00    9.50    6.14     7.45\n"
"Pi(0-1e9,10^9) = 22537866     8.15   4.28    3.92     3.02    3.64    2.28     2.05\n"
"Pi(1e18, 10^6) = 24280        0.65   0.46    0.34     0.48    0.60    0.27     0.35\n"
"Pi(1e18, 10^8) = 2414886      1.30   0.81    0.70     0.60    0.70    0.50     0.65\n"
"Pi(1e18, 10^9) = 24217085     3.50   1.80    1.51     1.26    1.50    0.95     1.25\n"
"Pi(0,    1e12) = 37607912018  500    270     224      200     220\n"
"Pi(1e14, 1e12) = 31016203073  790    420     354      295     320\n"
"Pi(1e16, 1e12) = 27143405794  1160   600     512      420     485\n"
"Pi(1e18, 1e12) = 24127637783  1500   760     622      520     640\n"
"Pi(1e19, 1e12) = 22857444126  1700   830     702      600     665\n";

static const char* Help = "\
	[B: Benchmark (0 - 12, 0 - 40)]\n\
	[D: D[T, R] dump time and result]\n\
	[M: Progress of calculating (0 - 20)]\n\
	[C: Cpu L1/L2 data cache size (L1:16-64, L2:256-4096)k]\n\
	[S: Set the sieve size (32 - 4096)k]\n\
	[L: Set sieve cache segs L(2-12)1, L(2-8)2 L(2-12)3, L(1-32)4]\n\
	[I: Info of programming]\n\
	[P: Print prime in [start, end]]\n\n\
Example:\n\
	1e16 10^10 s2048 c321 c2562";

#include <ctype.h>
#include <math.h>
#include <memory.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#if __x86_64__ || __amd64__ || _M_X64 || __amd64 || __x86_64
# define X86_64       1
#elif __i386__ || _M_IX86 || _X86_ || __i386
# define X86          1
#endif

#if (X86_64 || X86) && (__GNUC__ || _MSC_VER || __clang__ || __TINYC__ || __INTEL_COMPILER)
# define BIT_SCANF    1
#if X86 || _MSC_VER == 0
# define ASM_X86      1
#endif
#endif

// likely/unlikely
#if (__GNUC__ >= 4 || __clang__)
#    define ELIKELY(condition)   __builtin_expect(condition, 1)
#    define EUNLIKELY(condition) __builtin_expect(condition, 0)
#else
#    define ELIKELY(condition)   condition
#    define EUNLIKELY(condition) condition
#endif

#if _MSC_VER
typedef unsigned __int64 uint64;
typedef __int64 int64;
#else
typedef unsigned long long uint64;
typedef long long int64;
#endif

typedef unsigned int uint;
typedef unsigned char uchar;
typedef unsigned short ushort;

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

#ifndef WMH11
#define WHEEL210  210
#else
#define WHEEL210  2310
#endif

#ifndef M210 //big sieve wheel
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
typedef struct PrimeCall     PrimeCall;
#endif

#define MIN(a, b)         (a < b ? a : b)
#define MAX(a, b)         (a > b ? a : b)

enum ECONST
{
	ERAT_SMALL = 6, //4 - 16
	ERAT_MEDIUM = 4, //2 - 6
	ERAT_BIG = 8, //2 - 8
#ifndef CHAR_BIT
# define CHAR_BIT   8
#endif

	WHEEL30 = 30,
	PRIME_PRODUCT = 210 * 11 * 13 * 17 * 19,
	FIRST_INDEX = PRIME_PRODUCT / 9699690 + 7,
	NEXT_MULTIPLE = 0x5A28A6, //magic number from 0x799b799b. 4 bit shrift into 3 bit
	PI_65536 = 6542 + 1, //pi(2^16) + 1
	MAX_SEGMENT = 8 << 10,  //4Mkb
#ifndef UINT_MAX
	UINT_MAX = 0 - 1u,
#endif

#if WHEEL210 < 2310
	PWS = 48,
	SP_BIT = 6, //PWS < 2^SP_BIT < 210 [6-7]
	SI_BIT = 8, //SI_BIT >= SP_BIT [8 - 10]
	MAX_WHEEL_GAP = 11 - 1,   //max prime wheel 210 gap
#else
	PWS = 480,
	SP_BIT = 10, //PWS < 2^SP_BIT < 2310 [9-10]
	SI_BIT = 10, //SI_BIT >= SP_BIT [8 - 10]S
	MAX_WHEEL_GAP = 15 - 1,   //max prime wheel 210 gap
#endif

#if (L2_DCACHE_SIZE < 256 || L2_DCACHE_SIZE > 4096)
	L2_DCACHE_SIZE = 256,
#endif
#ifndef SIEVE_SIZE
	SIEVE_SIZE = 4096
#endif
};

enum EBUCKET
{
	UINT_PIMAX = 203280221, //= pi(2^32)
	MAX_BUCKET = 5465 * 2,  //= 0xffffffff / (256 * (1 << 10) * 3) + 4, 5465 (sqrtp * (8 + 2) / sieve_size + 2);
	WHEEL_SIZE = 1 << 12,   //= 2048 11: 16k, [10 - 13]
	MEM_WHEEL = WHEEL_SIZE * sizeof(int) * 2, //=32768
	MEM_BLOCK = (1 << 19) / WHEEL_SIZE, //19:4 MB
	MAX_STOCK = UINT_PIMAX / WHEEL_SIZE + MAX_BUCKET + MEM_BLOCK, //=55221
	MAX_POOL = UINT_PIMAX / (MEM_BLOCK * WHEEL_SIZE) + 100 //=487
};

enum EFLAG
{
	PRINT_TIME = 1 << ('T' - 'A'),
	SLOW_TEST = 1 << ('A' - 'A'),
	SAVE_DATA = 1 << ('F' - 'A'),
	PRINT_RET = 1 << ('R' - 'A')
};

enum EBITMASK
{
	BIT0 = 1 << 0, BIT1 = 1 << 1,
	BIT2 = 1 << 2, BIT3 = 1 << 3,
	BIT4 = 1 << 4, BIT5 = 1 << 5,
	BIT6 = 1 << 6, BIT7 = 1 << 7
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
	uint Bsegs;
};

struct WheelElement
{
#if WHEEL210 < 2310
	uchar WheelIndex;
	uchar Multiple;
	uchar Correct;
	uchar MaskBit;
#else
	ushort WheelIndex;
	ushort Multiple;
	ushort Correct;
	uchar MaskBit;
#endif
};

/**
struct WheelFirst
{
#if WHEEL210 < 2310
	uchar WheelIndex;
	uchar Multiple;
	uchar Correct;
#else
	ushort WheelIndex;
	ushort Multiple;
	ushort Correct;
#endif
};**/

struct WheelInit
{
#if WHEEL210 < 2310
	uchar WheelIndex;
	uchar MaskBit;
#else
	ushort WheelIndex;
	uchar MaskBit;
#endif
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
void initPrime();
//uint64 doSieve(const uint64 start, const uint64 end, PrimeCall* pcall);
uint setSieveSize(uint sieve_size);
void setCacheSegs(uint level, uint cachesegs);
void setCacheSize(uint level, uint cachecpu);
#endif

//common cache
typedef WheelElement WheelFirst;
static WheelFirst WheelFirst210[WHEEL210][PWS];
static WheelFirst WheelFirst30[WHEEL30][8];

static WheelInit WheelInit210[WHEEL210];
static WheelElement WheelData210[PWS][PWS];

static WheelInit WheelInit30[WHEEL30];
static uchar Pattern30[64];
static WheelElement WheelData30[8][8];
//presieved with prime <= 19
static uchar PreSieved[PRIME_PRODUCT / WHEEL30];

#if BIT_SCANF == 0
static uchar Lsb[1 << 16];
#endif

#if _MSC_VER >= 1400
# include <intrin.h>
#endif

#if POPCNT == 0
//number of bits 1 binary representation table in Range[0-2^16)
static uchar WordNumBit1[1 << 16];
#if (__GNUC__ || __clang__) && (X86 || X86_64)
# include <popcntintrin.h>
#endif
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
*******************************/

static BucketInfo_ BucketInfo;
static Stock* StockHead;
static Stock  StockCache[MAX_STOCK];

//big/bucket wheel pool
static Bucket_ Bucket[MAX_BUCKET];
static SievePrime* WheelPool[MAX_POOL];

//medium prime
typedef struct SievePrime MediumSieve;
static MediumSieve* MediumPrime;

//small cache
static SievePrime SmallPrime[PI_65536];

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
	1 << 3,
};

////////////////basic func //////////////////////////////

//the simple sieve of Eratosthenes p < 2^16
static int eratoSimple()
{
	int primes = 1;
	const uint maxp = (1 << 16) + 1;
	uchar bitarray[maxp >> 4] = { 0 };

	for (uint p = 3; p < maxp; p += 2) {
		if (0 == (bitarray[p >> 4] & (1 << (p / 2 & 7)))) {
			SmallPrime[primes++].Sp = p;
			for (uint j = p * p; j <= maxp; j += p * 2)
				bitarray[j >> 4] |= 1 << (j / 2 & 7);
		}
	}

	SmallPrime[primes++].Sp = maxp;
	return primes;
}

//The first presieved template, cross off the first 8th prime multiples
static void initBitTable()
{
	uint i = 0;
#if POPCNT == 0
	const uint nbitsize = sizeof(WordNumBit1) / sizeof(WordNumBit1[0]);
	for (i = 1; i < nbitsize; i++)
		WordNumBit1[i] = WordNumBit1[i >> 1] + (i & 1);
#endif

	const uchar pattern[] = { 1, 7, 11, 13, 17, 19, 23, 29 };
	for (i = 0; i < sizeof(Pattern30) / sizeof(Pattern30[0]); i++)
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

	for (int j = 0, wj = 0; j < WHEEL30; j++) {
		WheelInit30[j].MaskBit = 0;
		WheelInit30[j].WheelIndex = wj;
		if (j == Pattern30[wj])
			WheelInit30[j].MaskBit = 1 << (wj++);
	}

	for (i = 1; PRIME_PRODUCT % pattern[i] == 0; i++) {
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
	for (int i = 0; i < WHEEL30; i++) {
		for (int pi = 0; pi < psize; pi++) {
			int multiples = 1 - (i % 2);
			int next = i + multiples * Pattern30[pi];

			while (WheelInit30[next % WHEEL30].MaskBit == 0)
				next = i + Pattern30[pi] * (multiples += 2);

			int wi = WheelInit30[next % WHEEL30].WheelIndex;
			WheelFirst* wf = &WheelFirst30[i][pi];
			wf->WheelIndex = wi;
			wf->Multiple = (nextMultiple[wi] >> (pi * 4)) & 15;
			wf->Correct = multiples;
		}
	}

	for (int wi = 0; wi < psize; wi++) {
		for (int pi = 0; pi < psize; pi++) {
			int multiples = 2;
			int next = Pattern30[wi] + Pattern30[pi] * multiples;
			while (WheelInit30[next % WHEEL30].MaskBit == 0)
				next = Pattern30[wi] + Pattern30[pi] * (multiples += 2);

			WheelElement* we30 = &WheelData30[pi][wi];
			we30->WheelIndex = WheelInit30[next % WHEEL30].WheelIndex;
			we30->Multiple = multiples;
			we30->Correct = next / WHEEL30 - Pattern30[wi] / WHEEL30;
			we30->MaskBit = 1 << wi;
		}
	}
}

static void initWheel210()
{
	int wi = 0, i = 0;
	const int psize = PWS;
	uint wpattern[PWS];

#if WHEEL210 < 2310
	const ushort pattern[] = { 3, 5, 7, 2311 };
#else
	const ushort pattern[] = { 3, 5, 7, 11, 2311 };
#endif

	for (int j = 1; j < WHEEL210; j += 2) {
		WheelInit210[j - 1].WheelIndex = WheelInit210[j].WheelIndex = wi;
		int i = 0;
		while (i < sizeof(pattern) / sizeof(pattern[0]) && j % pattern[i] != 0)
			i++;
		if (i == sizeof(pattern) / sizeof(pattern[0]))
			wpattern[wi++] = j;
	}

	for (i = 0; i < WHEEL210; i++) {
		for (int pi = 0; pi < psize; pi++) {
			int multiples = 1 - (i % 2);
			uint next = i + wpattern[pi] * multiples;

			WheelFirst* wf = &WheelFirst210[i][pi];
			while (WheelInit210[next % WHEEL210].WheelIndex == WheelInit210[(next + 1) % WHEEL210].WheelIndex)
				next = i + wpattern[pi] * (multiples += 2);
			wf->WheelIndex = WheelInit210[next % WHEEL210].WheelIndex;
			wf->Correct = multiples;
			assert(multiples < MAX_WHEEL_GAP);
		}
	}

	for (int pi = 0; pi < psize; pi++) {
		for (wi = 0; wi < psize; wi++) {
			int multiples = 2;
			uint next = wpattern[wi] + wpattern[pi] * 2;

			while (WheelInit210[next % WHEEL210].WheelIndex == WheelInit210[(next + 1) % WHEEL210].WheelIndex) {
				next += wpattern[pi] * 2;
				multiples += 2;
			}

			assert(multiples <= MAX_WHEEL_GAP);
			WheelElement* we210 = &WheelData210[pi][wi];
			we210->Correct = next / WHEEL30 - wpattern[wi] / WHEEL30;
			we210->MaskBit = WheelInit30[wpattern[wi] % WHEEL30].MaskBit;
			we210->WheelIndex = WheelInit210[next % WHEEL210].WheelIndex;
			we210->Multiple = multiples * (WHEEL210 / WHEEL30);
		}
	}
}

#if __linux__ || __unix__
#include <sys/resource.h>
#elif _WIN32
#include <windows.h>
#endif

static int64 getTime()
{
#ifdef _WIN32
	FILETIME ptime[4] = { 0 };
	GetThreadTimes(GetCurrentThread(), &ptime[0], &ptime[1], &ptime[2], &ptime[3]);
	return (ptime[2].dwLowDateTime + ptime[3].dwLowDateTime) / 10000;
#elif __linux__ || __unix__
	struct rusage rup;
	getrusage(RUSAGE_SELF, &rup);
	uint64 sec = rup.ru_utime.tv_sec + rup.ru_stime.tv_sec;
	uint64 usec = rup.ru_utime.tv_usec + rup.ru_stime.tv_usec;
	return sec * 1000 + usec / 1000;
#elif _WIN32
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
		logn++;
		x /= n;
	}

	return logn;
}

//x^n < 2^64, n < 64
static uint64 mpow(const uint x, uint n)
{
	uint64 pown = 1;
	while (n--)
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
#if _WIN32
	unsigned long index;
#if _WIN64
	_BitScanForward64(&index, n);
#else
	_BitScanForward(&index, n);
#endif
#elif __GNUC__ || __clang__
#if X86_64
	//	uint index = __builtin_ffsll(n) - 1;
	uint index = __builtin_ctzll(n);
#else
	uint index = __builtin_ctzl(n);
#endif
#elif ASM_X86
	stype index;
#if X86_64
#if __GNUC__ || __TINYC__ || __clang__
	__asm__("bsfq %1, %0\n" : "=r" (index) : "rm" (n) : "cc");
#else
	__asm
	{
		bsfq eax, n
		mov index, eax
	}
#endif
#else
#if __GNUC__ || __TINYC__ || __clang__
	__asm__("bsfl %1, %0\n" : "=r" (index) : "rm" (n) : "cc");
#else
	__asm
	{
		bsf eax, n
		mov index, eax
	}
#endif
#endif
#endif

	return (uint)index;
}
#endif

//(n / p) < 2^32 is more efficient
inline static uint fastllDiv(const uint64 n, uint p)
{
#if ASM_X86 == 0
	p = (uint)(n % p);
#elif __GNUC__ || __TINYC__ || __clang__
	const uint loww = (uint)n, higw = (uint)(n >> 32);
	__asm__("divl %%ecx\n" : "=d" (p) : "d"(higw), "a"(loww), "c"(p));
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
		ret = ret * 10 + *str++ - '0';

	if (*str && isdigit(str[1])) {
		if (str[0] == '^') {
			ret = mpow((uint)ret, atoi(str + 1));
		}
		else if (str[0] == 'e' || str[0] == 'E') {
			if (ret == 0)
				ret = 1;
			ret *= mpow(10, atoi(str + 1));
		}
	}

	const char* ps = str;
	if ((ps = strchr(str, '+')) != NULL) {
		ret += atoint64(ps + 1);
	}
	else if ((ps = strchr(str, '-')) != NULL) {
		ret -= atoint64(ps + 1);
	}
	else if ((ps = strchr(str, '*')) != NULL) {
		ret *= atoi(ps + 1);
	}
	else if ((ps = strchr(str, '/')) != NULL) {
		ret /= atoi(ps + 1);
	}

	return ret;
}

//count number of bit 0 in binary representation of array, bytes % 8 = 0
static uint countBit0sArray(uint64 bitarray[], const uint bytes)
{
	uint bit1s = 0, bit2s = 0;
	uint loops = bytes / sizeof(uint64);
	while (loops-- > 0) {
		const uint64 qw = *bitarray++;
#if (__GNUC__ || __clang__ ) && X86_64 || POPCNT
		bit2s += _mm_popcnt_u64(qw);
#elif POPCNT
		bit2s += __builtin_popcountll(qw);
#else
		const uint hig = (uint)(qw >> 32);
		const uint low = (uint)(qw);
		//if (sizeof(uint64) != sizeof(uint)) //uint64 == uint ?
		bit1s += WordNumBit1[(ushort)hig] + WordNumBit1[hig >> 16];
		bit2s += WordNumBit1[(ushort)low] + WordNumBit1[low >> 16];
#endif
	}

	return bytes * CHAR_BIT - bit1s - bit2s;
}

//////////////////////// core code //////////////////////////////////
static void allocWheelBlock(const uint blocks)
{
#if 0
	SievePrime* pSprime = (SievePrime*)_aligned_malloc(MEM_BLOCK * MEM_WHEEL, MEM_WHEEL);
	WheelPool[BucketInfo.PoolSize++] = pSprime;
#else
	SievePrime* pSprime = (SievePrime*)malloc((MEM_BLOCK + 1) * MEM_WHEEL);
	WheelPool[BucketInfo.PoolSize++] = pSprime;
	pSprime = (SievePrime*)((size_t)pSprime + MEM_WHEEL - (size_t)pSprime % MEM_WHEEL);
#endif
	//assert (BucketInfo.PoolSize < sizeof(WheelPool) / sizeof(WheelPool[0]));
	//assert (BucketInfo.StockSize + blocks < sizeof(StockCache) / sizeof(StockCache[0]));

	Stock* pStock = StockCache + BucketInfo.StockSize;
	for (uint i = 0; i < MEM_BLOCK; i++) {
		pStock->Sprime = pSprime + WHEEL_SIZE * i;
		pStock->Next = pStock + 1;
		pStock++;
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
	BucketInfo.MaxBucket = range / sieve_size + 1;//overflow for big range
	BucketInfo.LoopSize = sqrtp / (sieve_size / MAX_WHEEL_GAP) + 2; //add the first, last bucket
	BucketInfo.Log2Size = ilog(sieve_size / WHEEL30, 2);
	BucketInfo.ModuloSize = (1 << BucketInfo.Log2Size) - 1;
	assert(BucketInfo.Log2Size <= (32 - SI_BIT) && SI_BIT >= SP_BIT);
	return 0;
}

static void initSmallWheel(const uint64 start, const uint maxp)
{
	const uint segsize = Threshold.L1Size * WHEEL30;
	uint64 offset = start - start % WHEEL210;

	for (uint p = SmallPrime[FIRST_INDEX].Sp, j = FIRST_INDEX; p < maxp; p = SmallPrime[++j].Sp) {
		const uint64 p2 = (uint64)p * p;
		if (p2 > offset && p2 > offset + segsize) //overflow
			offset += (p2 - offset) / segsize * segsize;

		uint sieve_index = p - (uint)(offset % p);
		if (p2 > offset)
			sieve_index = (uint)(p2 - offset);

		const int pi = WheelInit30[p % WHEEL30].WheelIndex;
		const WheelFirst& wf = WheelFirst30[sieve_index % WHEEL30][pi];
		sieve_index += wf.Correct * p;
#ifndef SM0
		SmallPrime[j].Si = 0 - ((sieve_index << 4) | wf.Multiple);
#else
		SmallPrime[j].Si = (sieve_index << 4) | wf.Multiple;
#endif
	}
}

static int segmentedSieve2(uchar bitarray[], const uint start, const uint sieve_size, bool bcopy);
static void initMediumWheel(uchar* bitarray, const uint sieve_size, const uint medium, const uint64 start)
{
	uint j = Threshold.L1Index, l1_maxp = Threshold.L1Maxp - Threshold.L1Maxp % WHEEL210;
	Threshold.L2Index = 0;

	uint segsize = MIN(sieve_size, Threshold.L2Size * WHEEL30);
	const int bytes = segmentedSieve2(bitarray, l1_maxp, medium - l1_maxp, false);

	stype mask = 0, * sbitarray = (stype*)bitarray;
	uint noffset = l1_maxp - sizeof(mask) * WHEEL30;
	const uint pn = 1 + bytes / sizeof(mask);

	for (uint i = 0; i < pn; ) {
		if (EUNLIKELY(mask == 0)) {
			mask = ~sbitarray[i++];
			noffset += sizeof(mask) * WHEEL30;
			continue;
		}

		const uint p = noffset + PRIME_OFFSET(mask); mask &= mask - 1;
		if (EUNLIKELY(Threshold.L2Index == 0 && p >= Threshold.L2Maxp)) {
			MediumPrime[j + 0].Sp = MediumPrime[j + 1].Sp = UINT_MAX;
			Threshold.L2Maxp = p;
			Threshold.L2Index = j += 2;
			segsize = sieve_size;
		}

		uint64 offset = start;
		const uint64 p2 = (uint64)p * p;
		if (EUNLIKELY(p2 > offset && p2 >= offset + segsize)) //overflow
			offset += (p2 - offset) / segsize * segsize;

		uint sieve_index = p - (uint)(offset % p);
		if (EUNLIKELY(p2 > offset))
			sieve_index = (uint)(p2 - offset);

		if (Threshold.L2Index == 0) {
			const int pi = WheelInit30[p % WHEEL30].WheelIndex;
			const WheelFirst& wf = WheelFirst30[sieve_index % WHEEL30][pi];
			sieve_index += wf.Correct * p;
#ifdef SL2
			MediumPrime[j + 0].Sp = p;
#ifndef SM0
			MediumPrime[j++].Si = 0 - ((sieve_index << 4) + wf.Multiple);
#else
			MediumPrime[j++].Si = (sieve_index << 4) + wf.Multiple;
#endif
#else
			MediumPrime[j + 0].Sp = (p / WHEEL30 << SI_BIT) + pi;
			MediumPrime[j++].Si = (sieve_index / WHEEL30 << SI_BIT) + wf.WheelIndex;
#ifdef SM0
			MediumPrime[j - 1].Si = -MediumPrime[j - 1].Si;
#endif
#endif
		}
		else {
			const int pi = WHEEL_INIT[p % WHEEL].WheelIndex;
			const WheelFirst& wf = WHEEL_FIRST[(sieve_index + (uint)(offset % WHEEL)) % WHEEL][pi];
			sieve_index += wf.Correct * p;
			MediumPrime[j + 0].Sp = (p / WHEEL << SI_BIT) + pi;
			MediumPrime[j++].Si = (sieve_index / WHEEL30 << SI_BIT) + wf.WheelIndex;
		}
	}

	MediumPrime[j + 0].Sp = MediumPrime[j + 1].Sp = UINT_MAX;
}

inline static void pushBucket(const uint offset, const uint sp, const ushort wi)
{
	const uint next_bucket = offset >> BucketInfo.Log2Size;
#ifndef B_R
	if (next_bucket >= BucketInfo.MaxBucket)
		return;
#endif

	SievePrime* pSprime = Bucket[next_bucket].Sprime++;
	if (EUNLIKELY((size_t)(pSprime) % MEM_WHEEL == 0)) {
		Bucket_* pbucket = Bucket + next_bucket;
#ifndef P_B
		if (EUNLIKELY(BucketInfo.CurStock-- == 0)) allocWheelBlock(1);
#else
		BucketInfo.CurStock--;
#endif
		//pop from free stock list
		Stock* psfree = StockHead; StockHead = psfree->Next;
		psfree->Next = pbucket->Head;
		pbucket->Head = psfree;
		pbucket->Sprime = (pSprime = psfree->Sprime) + 1;
	}

	//	*(uint64*)pSprime = ((uint64)sp << 32) | (offset & BucketInfo.ModuloSize) << SI_BIT | wi;
	pSprime->Sp = sp;
	pSprime->Si = (offset & BucketInfo.ModuloSize) << SI_BIT | wi;
}

static void initBigWheel(uchar* bitarray, uint medium, uint sqrtp, const uint64 start, const uint64 range)
{
	uint nextp = 0; uint64 remp = 0;
	if (ELIKELY(sqrtp < UINT_MAX)) sqrtp++; //watch overflow if sqrtp = 2^32 - 1

	const uint irange = (uint)((range >> 32) > WHEEL30 ? UINT_MAX : range / WHEEL30);
	for (uint l2_size = L2_DCACHE_SIZE * WHEEL30 << 10; medium < sqrtp; medium += l2_size) {
		if (l2_size > sqrtp - medium)
			l2_size = sqrtp - medium;

		const int bytes = segmentedSieve2(bitarray, medium, l2_size, nextp > 0);
		stype mask = 0, * sbitarray = (stype*)bitarray; //little endian
		uint offset = medium - sizeof(mask) * WHEEL30;
		const uint pn = 1 + bytes / sizeof(mask);

		for (uint j = 0; j < pn; ) {
			if (EUNLIKELY(mask == 0)) {
				mask = ~sbitarray[j++];
				offset += sizeof(mask) * WHEEL30;
				continue;
			}

			const uint p = offset + PRIME_OFFSET(mask); mask &= mask - 1;
#ifdef FDIV
			if (EUNLIKELY(p > nextp)) { //ugly & difficult to understand but efficient
				remp = start / (nextp = p + (uint64)p * p / (uint)(start >> 32));
				if (EUNLIKELY(p > nextp)) //overflow
					remp = start >> 32, nextp = UINT_MAX;
			}
			uint sieve_index = p - fastllDiv(start - remp * p, p);
#else
			uint sieve_index = p - (uint)(start % p);
#endif

#ifndef P_C
			if (EUNLIKELY(sieve_index > range)) continue;
#endif
			const uint sp = (p / WHEEL210 << SP_BIT) + WheelInit210[p % WHEEL210].WheelIndex;
			const uint modulo_210 = sieve_index % WHEEL210;
			const WheelFirst* wf = &WheelFirst210[modulo_210][sp % (1 << SP_BIT)];

			if (sizeof(int*) == sizeof(uint64))
				sieve_index = (sieve_index + wf->Correct * (uint64)p) / WHEEL30;
			else {
				sieve_index = (sieve_index / WHEEL210 + (sp >> SP_BIT) * wf->Correct) * (WHEEL210 / WHEEL30);
				sieve_index += (wf->Correct * (p % WHEEL210) + modulo_210) / WHEEL30;
			}

#ifndef B_R
			if (sieve_index > irange)
				continue;
#endif

#ifndef P_B
			pushBucket(sieve_index, sp, wf->WheelIndex);
#else
			const uint next_bucket = sieve_index >> BucketInfo.Log2Size;
			SievePrime* pSprime = Bucket[next_bucket].Sprime++;
			if (EUNLIKELY((size_t)(pSprime) % MEM_WHEEL == 0)) {
				Bucket_* pbucket = Bucket + next_bucket;
				if (EUNLIKELY(BucketInfo.CurStock-- == 0)) allocWheelBlock(1);
				//pop from free stock list
				Stock* psfree = StockHead; StockHead = psfree->Next;
				psfree->Next = pbucket->Head;
				pbucket->Head = psfree;
				pbucket->Sprime = (pSprime = psfree->Sprime) + 1;
			}
			pSprime->Sp = sp;
			pSprime->Si = (sieve_index & BucketInfo.ModuloSize) << SI_BIT | wf->WheelIndex;
#endif
		}
	}
}

static int crossSmall0(uchar bitarray[], const uchar* pend, const uint p, uint offset, uint multiples)
{
	uchar* ppbeg[8], wi;
	for (int i = 0; i < 8; i++) {
		wi = WheelInit30[offset % WHEEL30].WheelIndex;
		ppbeg[wi] = bitarray + offset / WHEEL30;
		offset += (multiples % 8) * p; multiples /= 8;
	}

#define OR_ADD(n) *ps##n |= BIT##n, ps##n += p
	uchar* ps0 = ppbeg[0], * ps1 = ppbeg[1], * ps2 = ppbeg[2], * ps3 = ppbeg[3];
	uchar* ps4 = ppbeg[4], * ps5 = ppbeg[5], * ps6 = ppbeg[6], * ps7 = ppbeg[7];

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
	for (int i = 0; i < 4; i++) {
		uchar* ps1 = bitarray + offset / WHEEL30;
		const uchar masks1 = WheelInit30[offset % WHEEL30].MaskBit;
		offset += (multiples % 8) * p; multiples /= 8;

		uchar* ps2 = bitarray + offset / WHEEL30;
		const uchar masks2 = WheelInit30[offset % WHEEL30].MaskBit;
		offset += (multiples % 8) * p; multiples /= 8;

		while (ps2 <= pend) {
			*ps2 |= masks2, * ps1 |= masks1;
			ps1 += p; ps2 += p;
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
		case 0: POSET(6,0,1, 10,0,2, 12,0,3,  16,0,4,  18,0,5,  22,0,6,  28,0, 7) break;
		case 1: POSET(4,0,7, 6,1,3,  10,2,2,  16,3,6,  18,4,1,  24,5,5,  28,6, 4) break;
		case 2: POSET(2,0,6, 6,2,1,  8,2,7,   12,4,3,  18,6,5,  20,7,2,  26,9, 4) break;
		case 3: POSET(4,1,6, 6,2,5,  10,4,2,  12,5,1,  16,6,7,  22,9,4,  24,10,3) break;
		case 4: POSET(6,3,3, 8,4,4,  14,7,7,  18,10,1, 20,11,2, 24,13,5, 26,14,6) break;
		case 5: POSET(4,2,4, 10,6,2, 12,7,5,  18,11,3, 22,13,7, 24,15,1, 28,17,6) break;
		case 6: POSET(2,1,4, 6,4,5,  12,9,1,  14,10,6, 20,15,2, 24,18,3, 26,19,7) break;
		case 7: POSET(2,1,7, 8,7,6,  12,11,5, 14,13,4, 18,17,3, 20,19,2, 24,23,1) break;
	}

	return (int)(ps - pbeg) * WHEEL30;
}

#define MEDIUM_SET(n) \
	we##n = wd##n[(int)we##n.WheelIndex]; \
	bitarray[(int)offset##n] |= we##n.MaskBit; \
	offset##n += we##n.Correct + we##n.Multiple * wi##n

#define MEDIUM_SET2() \
	we1 = wd1[(int)we1.WheelIndex]; \
	we2 = wd2[(int)we2.WheelIndex]; \
	bitarray[(int)offset1] |= we1.MaskBit; \
	bitarray[(int)offset2] |= we2.MaskBit; \
	offset1 += we1.Correct + we1.Multiple * wi1; \
	offset2 += we2.Correct + we2.Multiple * wi2; \

//cross out one medium prime from array
static void crossMedium1(uchar bitarray[], const uint sieve_byte, MediumSieve* pSprime)
{
	const uint si = pSprime->Si;
	int offset = (si >> SI_BIT) - sieve_byte;

	const uint wi = pSprime->Sp >> SI_BIT;
	//WheelElement* wd = WheelData30[pSprime->Sp % (1 << SI_BIT)];
	WheelElement wd[8]; /*= WheelData30[sp % (1 << SI_BIT)];*/ memcpy(wd, WheelData30[pSprime->Sp % (1 << SI_BIT)], sizeof(wd));
	WheelElement we; we.WheelIndex = si % (1 << SI_BIT);

	while (offset < 0) {
		MEDIUM_SET();
	}

	pSprime->Si = (offset << SI_BIT) | we.WheelIndex;
}

//cross out 2 medium prime from array
static void crossMedium2(uchar bitarray[], uint sieve_byte, uint maxsp, MediumSieve* pSprime)
{
	//stack stack local
	WheelElement wheelData30[8][8]; memcpy(wheelData30, WheelData30, sizeof(wheelData30));

	while (pSprime[1].Sp < maxsp) {
		const uint si1 = pSprime[0].Si, sp1 = pSprime[0].Sp;
		const uint si2 = pSprime[1].Si, sp2 = pSprime[1].Sp;

		uint wi1 = sp1 >> SI_BIT, offset1 = (si1 >> SI_BIT) - sieve_byte;
		uint wi2 = sp2 >> SI_BIT, offset2 = (si2 >> SI_BIT) - sieve_byte;

#if 1
		WheelElement* wd2 = wheelData30[sp2 % (1 << SI_BIT)];// memcpy(wd2, WheelData30[sp2 % (1 << SI_BIT)], sizeof(wd2));
		WheelElement* wd1 = wheelData30[sp1 % (1 << SI_BIT)];// memcpy(wd1, WheelData30[sp1 % (1 << SI_BIT)], sizeof(wd1));
#else
		WheelElement wd1[8]; WheelElement wd2[8];
		memcpy(wd2, wheelData30[sp2 % (1 << SI_BIT)], sizeof(wd2));
		memcpy(wd1, wheelData30[sp1 % (1 << SI_BIT)], sizeof(wd1));
#endif

		WheelElement we1; we1.WheelIndex = si1 % (1 << SI_BIT);
		WheelElement we2; we2.WheelIndex = si2 % (1 << SI_BIT);

		while ((int)(offset1 & offset2) < 0) {
			MEDIUM_SET2();
		}

		while ((int)offset2 < 0) {
			MEDIUM_SET(2);
		}

		while ((int)offset1 < 0) {
			MEDIUM_SET(1);
		}

		pSprime[0].Si = offset1 << SI_BIT | we1.WheelIndex;
		pSprime[1].Si = offset2 << SI_BIT | we2.WheelIndex;
		pSprime += 2;
	}

	if (pSprime[0].Sp < maxsp)
		crossMedium1(bitarray, sieve_byte, pSprime);
}

static void crossMedium4(uchar bitarray[], uint sieve_byte, uint maxsp, MediumSieve* pSprime)
{
	WheelElement wheelData[PWS][PWS]; memcpy(wheelData, WheelData210, sizeof(WHEEL_DATA));

	while (pSprime[1].Sp < maxsp) {
		uint& si1 = pSprime[0].Si, sp1 = pSprime[0].Sp;
		uint& si2 = pSprime[1].Si, sp2 = pSprime[1].Sp;
		pSprime += 2;

		uint wi1 = sp1 >> SI_BIT, offset1 = (si1 >> SI_BIT) - sieve_byte;
		uint wi2 = sp2 >> SI_BIT, offset2 = (si2 >> SI_BIT) - sieve_byte;

		//use stack local
#if 1
		WheelElement* wd2 = wheelData[sp2 % (1 << SI_BIT)];// memcpy(wd2, WheelData30[sp2 % (1 << SI_BIT)], sizeof(wd2));
		WheelElement* wd1 = wheelData[sp1 % (1 << SI_BIT)];// memcpy(wd1, WheelData30[sp1 % (1 << SI_BIT)], sizeof(wd1));
#else
		WheelElement wd1[48]; WheelElement wd2[48];
		memcpy(wd2, wheelData[sp2 % (1 << SI_BIT)], sizeof(wd2));
		memcpy(wd1, wheelData[sp1 % (1 << SI_BIT)], sizeof(wd1));
#endif

		WheelElement we1; we1.WheelIndex = si1 % (1 << SI_BIT);
		WheelElement we2; we2.WheelIndex = si2 % (1 << SI_BIT);

		while ((int)(offset1 & offset2) < 0) {
			MEDIUM_SET2();
		}

		while ((int)offset2 < 0) {
			MEDIUM_SET(2);
		}

		while ((int)offset1 < 0) {
			MEDIUM_SET(1);
		}

		si1 = offset1 << SI_BIT | we1.WheelIndex;
		si2 = offset2 << SI_BIT | we2.WheelIndex;
	}

	uint& si = pSprime->Si;
	int offset = (si >> SI_BIT) - sieve_byte;

	const uint wi = pSprime->Sp >> SI_BIT;
	WheelElement* wd = WHEEL_DATA[pSprime->Sp % (1 << SI_BIT)];
	WheelElement we; we.WheelIndex = si % (1 << SI_BIT);

	while (offset < 0) {
		MEDIUM_SET();
	}
	si = (offset << SI_BIT) | we.WheelIndex;
}

//more efficient for medium if prime < sieve_size / 30
static MediumSieve* crossMediumW30(uchar bitarray[], const uint sieve_byte, const uint minsp, MediumSieve* pSprime)
{
	//use stack local
	uchar pattern30[sizeof(Pattern30)]; memcpy(pattern30, Pattern30, sizeof(pattern30));
	WheelElement wheelData30[8][8]; memcpy(wheelData30, WheelData30, sizeof(wheelData30));

	const uchar* pend = bitarray + sieve_byte;
	for (uint sp = pSprime->Sp; sp < minsp; sp = (++pSprime)->Sp) {
		const uint wi = sp >> SI_BIT, pi = sp % (1 << SI_BIT);
		uchar* ps = bitarray + (pSprime->Si >> SI_BIT);
		const uint p = wi * WHEEL30 + pattern30[pi];

		WheelElement* wd = wheelData30[pi];
		WheelElement* we = wd + pSprime->Si % (1 << SI_BIT);
		const uchar* pmin = pend + p;

		int nwi = 0, lw = -1;
#ifdef W2
		for (int i = 0; i < 4; i += 1) {
#else
		for (int i = 0; i < 8; i += 1) {
#endif
			uchar* ps1 = ps;
			const uchar masks1 = we->MaskBit;
			ps += we->Correct + we->Multiple * wi;
			uchar cw1 = we->WheelIndex; we = wd + we->WheelIndex;

#ifdef W2
			uchar* ps2 = ps;
			const uchar masks2 = we->MaskBit;
			ps += we->Correct + we->Multiple * wi;
			uchar cw2 = we->WheelIndex; we = wd + we->WheelIndex;

			while (ps2 < pend) {
				*ps1 |= masks1; *ps2 |= masks2;
				ps1 += p; ps2 += p;
			}

			if (ps1 < pend) {
				*ps1 |= masks1; ps1 += p;
			}
#else
			uchar* ps2 = ps1 + p;
			while (ps2 < pend) {
				*ps1 |= masks1; ps1 += p * 2;
				*ps2 |= masks1; ps2 += p * 2;
			}
			if (ps1 < pend) {
				*ps1 |= masks1; ps1 += p;
			}
#endif

			if (ps1 < pmin) {
				pmin = ps1, nwi = lw;
			}
#ifndef W2
			lw = cw1;
#else
			if (ps2 < pmin) {
				pmin = ps2, nwi = cw1;
			}
			lw = cw2;
#endif
		}

		if (nwi < 0)
			nwi = lw;
		pSprime->Si = ((uint)(pmin - pend) << SI_BIT) | nwi;
	}

	return pSprime;
}

//cross out big prime's multiples from bucket
static inline void crossBig1(uchar bitarray[], const uint sieve_size, const SievePrime* pSprime)
{
	uint offset = pSprime->Si, sp = pSprime->Sp;
	const WheelElement* wd = WheelData210[sp % (1 << SP_BIT)];
	const WheelElement* we = wd + offset % (1 << SI_BIT);
	bitarray[offset >>= SI_BIT] |= we->MaskBit;
	offset += we->Correct + we->Multiple * (sp >> SP_BIT);

#ifdef FB
	if (EUNLIKELY(offset < sieve_size)) {
		we = wd + we->WheelIndex;
		bitarray[offset] |= we->MaskBit;
		offset += we->Correct + we->Multiple * (sp >> SP_BIT);
	}
#endif

	pushBucket(offset, sp, we->WheelIndex);
}

//cross out 2 big prime's multiples from bucket, 15% improvement
static void crossBig2(uchar bitarray[], const uint sieve_size, const SievePrime* pSprime)//, const WheelElement wheelData210[PWS][PWS])
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

#ifdef FB
	if (EUNLIKELY(offset1 < sieve_size)) {
		we1 = wd1 + we1->WheelIndex;
		bitarray[offset1] |= we1->MaskBit;
		offset1 += we1->Correct + we1->Multiple * (sp1 >> SP_BIT);
	}
	if (EUNLIKELY(offset2 < sieve_size)) {
		we2 = wd2 + we2->WheelIndex;
		bitarray[offset2] |= we2->MaskBit;
		offset2 += we2->Correct + we2->Multiple * (sp2 >> SP_BIT);
	}
#endif

	pushBucket(offset1, sp1, we1->WheelIndex);
	pushBucket(offset2, sp2, we2->WheelIndex);
}

static void preSieve(uchar bitarray[], const uint64 start, const uint sieve_size)
{
	const uint offset = (uint)(start % PRIME_PRODUCT) / WHEEL30;
	const uint bits = sieve_size / WHEEL30 * 8 + WheelInit30[sieve_size % WHEEL30].WheelIndex;
	uint bytes = (bits + 7) / CHAR_BIT, remains = sizeof(PreSieved) - offset;
	if (remains > bytes) {
		memcpy(bitarray, PreSieved + offset, bytes);
	}
	else {
		memcpy(bitarray, PreSieved + offset, remains);
		for (; remains + sizeof(PreSieved) < bytes; remains += sizeof(PreSieved))
			memcpy(bitarray + remains, PreSieved, sizeof(PreSieved));
		memcpy(bitarray + remains, PreSieved, bytes - remains);
	}

	//wheel pattern < WHEEL30 is prime except 1
	if (EUNLIKELY(start == 0))
		bitarray[0] = BIT0;

	//pack the last byte with bit 1
	if (bits % CHAR_BIT != 0)
		bitarray[bits / CHAR_BIT] |= ~((1 << (bits % CHAR_BIT)) - 1);
}

static int packWheel30(uchar bitarray[], uint offset, uint p)
{
	uint multiples = offset % 16 * 3; offset >>= 4;
	multiples = (NEXT_MULTIPLE >> multiples) | (NEXT_MULTIPLE << (24 - multiples));
	for (int i = 0; i < 24; i += 3) {
		const uchar mask = WheelInit30[offset % WHEEL30].MaskBit;
		if (mask == BIT0)
			break;
		bitarray[offset / WHEEL30] |= mask;
		offset += (multiples >> i) % 8 * p;
	}
	return offset;
}

//sieve prime multiples in [start, start + sieve_size) with small algorithm
static void eratSieveLeve1(uchar bitarray[], const uint64 start, const uint sieve_size, uint maxp)
{
	if (start + sieve_size < ((uint64)maxp) * maxp)
		maxp = isqrt(start + sieve_size) + 1;

	const uchar* pend = bitarray + sieve_size / WHEEL30;
	for (uint p = SmallPrime[FIRST_INDEX].Sp, j = FIRST_INDEX; p < maxp; p = SmallPrime[++j].Sp) {
		uint& offset = SmallPrime[j].Si;
#ifdef SM0
		uint multiples = offset % 16 * 3; uint mask = multiples; offset >>= 4;
		multiples = (NEXT_MULTIPLE >> multiples) | (NEXT_MULTIPLE << (24 - multiples));
		offset += crossSmall0(bitarray, pend, p, offset, multiples) + 30 * p - sieve_size;
		offset = offset * 16 + mask;
#else
		if (EUNLIKELY((int)offset < 0))
			offset = packWheel30(bitarray, 0 - offset, p);
		offset += crossSmall3(bitarray + offset / WHEEL30, pend, p) - sieve_size;
#endif
	}
}

//sieve prime multiples in [start, start + sieve_size) with small algorithm
static void eratSieveLeve2(uchar bitarray[], const uint64 start, const uint sieve_size, uint maxp)
{
	if (start + sieve_size < ((uint64)maxp) * maxp)
		maxp = isqrt(start + sieve_size) + 1;

	const uchar* pend = bitarray + sieve_size / WHEEL30;
	for (uint p = MediumPrime[Threshold.L1Index].Sp, j = Threshold.L1Index; p < maxp; p = MediumPrime[++j].Sp) {
		uint& offset = MediumPrime[j].Si;
#ifdef SM0
		uint multiples = offset % 16 * 3, mask = multiples; offset >>= 4;
		multiples = (NEXT_MULTIPLE >> multiples) | (NEXT_MULTIPLE << (24 - multiples));
		offset += crossSmall0(bitarray, pend, p, offset, multiples) + 30 * p - sieve_size;
		offset = offset * 16 + mask;
#else
		if ((int)offset < 0)
			offset = packWheel30(bitarray, 0 - offset, p);
		offset += crossSmall3(bitarray + offset / WHEEL30, pend, p) - sieve_size;
#endif
	}
}

//sieve prime multiples in [start, start + sieve_size) with medium algorithm
static void eratSieveMedium1(uchar bitarray[], const uint64 start, const uint sieve_size, const uint spi, uint maxp)
{
	if (start + sieve_size < (uint64)maxp * maxp)
		maxp = isqrt(start + sieve_size) + 1;

	MediumSieve* pSprime = MediumPrime + spi;
	const uint sieve_byte = sieve_size / WHEEL30 + sieve_size % WHEEL30;
	const uint maxsp = (maxp / WHEEL30 << SI_BIT) + WheelInit30[maxp % WHEEL30].WheelIndex;

	//	if (maxp <= Threshold.L2Maxp) {
	const uint minsb = sieve_byte / 32 << SI_BIT;
	const uint minsp = MIN(maxsp, minsb);
	pSprime = crossMediumW30(bitarray, sieve_byte, minsp, pSprime);
	//	}
	crossMedium2(bitarray + sieve_byte, sieve_byte, maxsp, pSprime);
}

//sieve prime multiples in [start, start + sieve_size) with medium algorithm
static void eratSieveMedium2(uchar bitarray[], const uint64 start, const uint sieve_size, const uint spi, uint maxp)
{
	if (start + sieve_size < (uint64)maxp * maxp)
		maxp = isqrt(start + sieve_size) + 1;

	const uint sieve_byte = sieve_size / WHEEL30 + (sieve_size % WHEEL30 ? 1 : 0);
	const uint maxsp = (maxp / WHEEL << SI_BIT) + WHEEL_INIT[maxp % WHEEL].WheelIndex;
	crossMedium4(bitarray + sieve_byte, sieve_byte, maxsp, MediumPrime + spi);
}

//This implementation uses a sieve array with WHEEL210 numbers per byte and
//a modulo wheel that skips multiples of 2, 3, 5 and 7.
static void eratSieveBig(uchar bitarray[], const uint sieve_size)
{
	//use stack local
	//WheelElement wheelData210[PWS][PWS];
	//memcpy(wheelData210, WheelData210, sizeof(wheelData210));

	uint loops = (size_t)Bucket[0].Sprime % MEM_WHEEL / sizeof(SievePrime);
	if (EUNLIKELY(loops == 0))
		loops = WHEEL_SIZE;

	for (Stock* pStock = Bucket[0].Head; pStock != NULL; loops = WHEEL_SIZE) {
		//push into free block list
		Stock* pnext = pStock->Next; pStock->Next = StockHead; StockHead = pStock;
		pStock = pnext;

		SievePrime* pSprime = StockHead->Sprime;
		while (loops >= 2) {
#ifndef B2
			crossBig2(bitarray, sieve_size, pSprime);
#else
			crossBig1(bitarray, sieve_size, pSprime + 0);
			crossBig1(bitarray, sieve_size, pSprime + 1);
#endif
			loops -= 2;
			pSprime += 2;
		}
		if (loops)
			crossBig1(bitarray, sieve_size, pSprime);

		BucketInfo.CurStock++;
	}

	BucketInfo.MaxBucket--;
	memmove(Bucket, Bucket + 1, BucketInfo.LoopSize * sizeof(Bucket[0]));
}

static uint segmentProcessed(uchar bitarray[], const uint64 start, uint bytes, PrimeCall* pcall)
{
	//pack last qword
	if (bytes % sizeof(uint64)) {
		memset(bitarray + bytes, ~0, sizeof(uint64));
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

static void eratSieveSmall(uchar bitarray[], const uint64 start, const uint sieve_size, const uint l1_maxp)
{
	for (uint offset = 0, l1_size = Threshold.L1Size * WHEEL30; offset < sieve_size; offset += l1_size) {
		if (EUNLIKELY(l1_size + offset > sieve_size))
			l1_size = sieve_size - offset;
		eratSieveLeve1(bitarray + offset / WHEEL30, start + offset, l1_size, l1_maxp);
	}
}

//sieve prime multiples in [start, start + sieve_size)
static int segmentedSieve(uchar bitarray[], const uint64 start, const uint sieve_size)
{
	const uint l1_maxp = Threshold.L1Maxp;
	const uint l2_maxp = Threshold.L2Maxp;
	const uint medium2 = Threshold.Medium + 1;

	preSieve(bitarray, start, sieve_size);

	//or last l1_size seg with the first seg for overflow
	const uint copy_from = Config.SieveSize;
#ifndef SL2
	for (uint i = 0; i < l1_maxp; i += sizeof(uint64)) {
#else
	for (uint i = 0; i < l2_maxp; i += sizeof(uint64)) {
#endif
		uint64& c = *(uint64*)(bitarray + i + copy_from);
		*(uint64*)(bitarray + i) |= c; c = 0;
	}

	for (uint offset = 0, l2_size = Threshold.L2Size * WHEEL30; offset < sieve_size; offset += l2_size) {
		if (EUNLIKELY(l2_size + offset > sieve_size))
			l2_size = sieve_size - offset;

		eratSieveSmall(bitarray + offset / WHEEL30, start + offset, l2_size, l1_maxp);
#ifdef SL2
		eratSieveLeve2(bitarray + offset / WHEEL30, start + offset, l2_size, l2_maxp);
#else
		eratSieveMedium1(bitarray + offset / WHEEL30, start + offset, l2_size, Threshold.L1Index, l2_maxp);
#endif
	}

	if (medium2 > l2_maxp)
#if M210
		eratSieveMedium2(bitarray, start, sieve_size, Threshold.L2Index, medium2);
#else
		eratSieveMedium1(bitarray, start, sieve_size, Threshold.L2Index, medium2);
#endif

	if (start >= Threshold.BucketStart) {
		const uint big_size = 1 << BucketInfo.Log2Size;
		for (uint offset = 0; offset < copy_from; offset += big_size)
			eratSieveBig(bitarray + offset, big_size);
	}

	return sieve_size / WHEEL30 + (7 + WheelInit30[sieve_size % WHEEL30].WheelIndex) / 8;
}

static int segmentedSieve2(uchar bitarray[], const uint start, const uint sieve_size, bool bcopy)
{
	const uint sqrtp = isqrt(start + sieve_size);
	const uint l1_maxp = Threshold.L1Maxp;
	preSieve(bitarray, start, sieve_size);

#ifndef FSL1
	const uint copy_from = L2_DCACHE_SIZE << 10;
	if (sieve_size <= copy_from * WHEEL30 && bcopy) {
		for (uint i = 0; i < l1_maxp; i += sizeof(uint64)) {
			uint64* pc = (uint64*)(bitarray + i + copy_from);
			*(uint64*)(bitarray + i) |= *pc; *pc = 0;
		}
	}
	eratSieveSmall(bitarray, start, sieve_size, l1_maxp);
#else
	for (uint offset = 0, l1_size = Threshold.L1Size * WHEEL30; offset < sieve_size; offset += l1_size) {
		if (l1_size + offset > sieve_size)
			l1_size = sieve_size - offset;

		//corss out sieve size = L1_SIZE
		const uint64 newstart = start + offset;
		const uint lsqrtp = isqrt(newstart + l1_size) + 1;
		uchar* pstart = bitarray + offset / WHEEL30;
		const uchar* pend = bitarray + (offset + l1_size) / WHEEL30;
		const uint maxp = MIN(lsqrtp, l1_maxp);

		for (uint p = SmallPrime[FIRST_INDEX].Sp, j = FIRST_INDEX; p < maxp; p = SmallPrime[++j].Sp) {
			uint& offset1 = SmallPrime[j].Si;
#ifndef SM0
			if ((int)offset1 < 0)
				offset1 = 0 - offset1;
#endif
			uint mask = multiples = offset1 % 16 * 3; offset1 >>= 4;
			multiples = (NEXT_MULTIPLE >> multiples) | (NEXT_MULTIPLE << (24 - multiples));
			offset1 += crossSmall0(pstart, pend, p, offset1, multiples) + 30 * p - l1_size;
			offset = offset * 16 + mask;
		}
	}
#endif

	//	assert(l1_maxp > start / l1_maxp);
	for (uint j = Threshold.L1Index, p = l1_maxp; p <= sqrtp; p = SmallPrime[++j].Sp) {
		uint& offset = SmallPrime[j].Si;
		if ((int)offset <= 0) {
			offset = p - start % p;
			if (p * p > start)
				offset = p * p - start;
		}

		const WheelFirst& wf = WheelFirst30[offset % WHEEL30][WheelInit30[p % WHEEL30].WheelIndex];
		const uint multiples = (NEXT_MULTIPLE >> (wf.Multiple * 3)) | (NEXT_MULTIPLE << (24 - wf.Multiple * 3));
		offset = crossSmall1(bitarray, bitarray + sieve_size / WHEEL30, p, offset + wf.Correct * p, multiples);
	}

	uint bytes = sieve_size / WHEEL30 + (7 + WheelInit30[sieve_size % WHEEL30].WheelIndex) / 8;
	if (bytes % sizeof(uint64)) {
		memset(bitarray + bytes, ~0, sizeof(uint64));
		bytes += sizeof(uint64) - bytes % sizeof(uint64);
	}

	return bytes;
}

static void setL1Index()
{
#if WHEEL210 < 2310
	if (Threshold.L1Maxp > 65521) //
		Threshold.L1Maxp = 65521; //the bigest prime % WHEEL210 = 1 and < 2^16
#else
	if (Threshold.L1Maxp > 57751) //57751
		Threshold.L1Maxp = 57751; //the bigest prime % 2310 = 1 and < 2^16
#endif

	for (uint p = 0, j = 300; ; p = SmallPrime[++j].Sp) {
		if (p >= Threshold.L1Maxp && p % 210 == 1) {
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
		Threshold.L1Size = (1 << ilog(cache, 2)) << 10;
		Threshold.L1Maxp = Threshold.L1Size / Config.L1Segs;
		Threshold.L2Size = Threshold.L2Size / Threshold.L1Size * Threshold.L1Size;
		Threshold.L2Maxp = Threshold.L2Size / Config.L2Segs;
		setL1Index();
	}
	else if (level == 2 && cache >= (Threshold.L1Size >> 10) && cache <= MAX_SEGMENT) {
		//Threshold.L2Size = (cache << 10) / Threshold.L1Size * Threshold.L1Size;
		Threshold.L2Size = (1 << ilog(cache, 2)) << 10;
		Threshold.L2Maxp = Threshold.L2Size / Config.L2Segs;
	}
}

void setCacheSegs(uint level, uint segs)
{
	if (segs == 0 || segs > 20)
		return;

	if (level == 1) {
		Config.L1Segs = segs;
		Threshold.L1Maxp = Threshold.L1Size / segs;
		setL1Index();
	}
	else if (level == 2) {
		Config.L2Segs = segs;
		Threshold.L2Maxp = Threshold.L2Size / segs;
	}
	else if (level == 3) {
		Config.Msegs = segs;
	}
	else if (level == 4) {
		Config.Bsegs = 1 << ilog(segs, 2);
	}
}

//sieve_size : 32k - 8m
uint setSieveSize(uint sieve_size)
{
	const uint l1_size = Threshold.L1Size >> 10, l2_Size = Threshold.L2Size >> 10;
	if (sieve_size <= MAX_SEGMENT && sieve_size > l2_Size) {
		sieve_size = (sieve_size / l2_Size) * l2_Size << 10;// 1 << (ilog(sieve_size, 2) + 10);
	}
	else if (sieve_size <= l2_Size && sieve_size >= l1_size) {
		sieve_size = (sieve_size / l1_size) * l1_size << 10;
	}
	else if (sieve_size <= (MAX_SEGMENT >> 10) && sieve_size > 0) {
		sieve_size = sieve_size << 20;
	}
	else {
		sieve_size = SIEVE_SIZE << 10;
	}

	return Config.SieveSize = sieve_size;
}

static int checkSmall(const uint64 start, const uint64 end, PrimeCall * pcall)
{
	int primes = 0;
	for (int i = 0, p = 2; p < 7; p = SmallPrime[++i].Sp) {
		if (start <= p && p <= end) {
			primes++;
			if (pcall && pcall->func) {
				pcall->primes += 1;
				pcall->func(pcall->data, p);
			}
		}
	}

	return primes;
}

static uint64 pi(uchar * bitarray, uint64 start, uint64 end, PrimeCall * pcall)
{
	const int64 ts = getTime();
	uint64 primes = checkSmall(start, end, pcall);

	uint align210 = (uint)(start % WHEEL210);
	start -= align210;

	if (EUNLIKELY(++end == 0)) end--; //overflow if end = 2^64-1

	//pi(n) ~= n/log(n), complexty ~= n*log(log(n)), replaced by n*log(n) more accurate
#ifndef PTIME
	const double logr = 1.2;
	double lge = end * log((double)end), lgs = start * log(start + 10.0);
	double pie = PI(end, logr), pis = PI(start, logr);
#endif

	for (uint si = 0, sieve_size = Config.SieveSize * WHEEL30; start < end; start += sieve_size) {
		if (EUNLIKELY(sieve_size > end - start))
			sieve_size = (uint)(end - start);

		const uint bytes = segmentedSieve(bitarray, start, sieve_size);
		if (EUNLIKELY(align210 > 0)) { //why brother ?
			memset(bitarray, ~0, align210 / WHEEL30);
			bitarray[align210 / WHEEL30] |= (1 << WheelInit30[align210 % WHEEL30].WheelIndex) - 1;
			align210 = 0;
		}

		primes += segmentProcessed(bitarray, start, bytes, pcall);

#ifndef PTIME
		if (EUNLIKELY((si++ & Config.Progress) == 15)) {
			const double cur = start + sieve_size;
			double lgc = cur * log(cur), pic = PI(cur, logr);
			double tratio = (lgc - lgs) / (lge - lgs) * 100;
			double pratio = (pic - pis) / (pie - pis);
			double timeuse = (getTime() - ts) / (10 * tratio);
			const uint64 picount = (int64)((int64)primes / pratio);
			if (ELIKELY(timeuse < 3600))
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
	const uint64 pown = mpow(10, logn);
	if (n % pown == 0)
		sprintf(buff, "%de%d", (int)(n / pown), logn);
	else if (n % (pown / 10) == 0)
		sprintf(buff, "%de%d", (int)(n / (pown / 10)), logn - 1);
	else if ((n & (n - 1)) == 0)
		sprintf(buff, "2^%d", ilog(n, 2));
	else if (n > 1000000000l) {
		uint64 r = n - (int)(n / pown) * pown;
		const int logr = ilog(r, 10);
		const uint64 powr = mpow(10, logr);
		if (r % powr == 0 && logr > 4)
			sprintf(buff, "%de%d+%de%d", (int)(n / pown), logn, (int)(r / powr), logr);
		else if (r % powr == 0)
			sprintf(buff, "%de%d+%d", (int)(n / pown), logn, (int)r);
	}
}

static void printResult(const uint64 start, const uint64 end, uint64 primes)
{
	char buff[128] = { 0 };
	char begbuff[40] = { 0 }, endbuff[40] = { 0 }, rangebuff[40] = { 0 };
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
	}
	else if (begbuff[0] && endbuff[0]) {
		if (rangebuff[0] && range <= start / 10)
			sprintf(buff, "%s,+%s", begbuff, rangebuff);
		else
			sprintf(buff, "%s,%s", begbuff, endbuff);
	}
	else if (rangebuff[0]) {
		if (begbuff[0])
			sprintf(buff, "%s,+%s", begbuff, rangebuff);
		else if (endbuff[0])
			sprintf(buff, "%s-%s,%s", endbuff, rangebuff, endbuff);
		else
			sprintf(buff, "%llu,+%s", start, rangebuff);
	}
	else {
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
	}
	else
		offset = start;

	if (medium > sqrtp)
		medium = sqrtp;
	medium += WHEEL210 - medium % WHEEL210;

	Threshold.BucketStart = offset - offset % WHEEL210;
	Threshold.Medium = medium;

	return medium;
}

static uint adjustConfig(uint sieve_size, const uint sqrtp)
{
	uint medium = (sieve_size * WHEEL30) / Config.Msegs + 1;
	if (sieve_size < Threshold.L2Size && sqrtp > medium)
		sieve_size = setSieveSize(SIEVE_SIZE);
	else if (sieve_size >= Config.Msegs << 21)
		sieve_size = setSieveSize(SIEVE_SIZE);

	if ((sieve_size & (sieve_size - 1)) != 0 && sqrtp > (sieve_size * WHEEL30) / Config.Msegs)
		sieve_size = setSieveSize(1 << ilog(sieve_size >> 10, 2));

	medium = (sieve_size * WHEEL30) / Config.Msegs + 1;
#ifndef FB
	if (sqrtp > medium) {
		const int l2segs = sieve_size / (Threshold.L2Size);
		assert(l2segs > 0 && sieve_size >= (128 << 10));
		if (Config.Bsegs * 2 < Config.Msegs)
			Config.Bsegs = 1 << ilog(Config.Msegs, 2);
		if (Config.Bsegs > l2segs) {
			Config.Bsegs = l2segs;
			Config.Msegs = 2;
		}
	}
#endif

	return sieve_size;
}

static void allocMedium(const uint medium)
{
	if (MediumPrime != NULL && MediumPrime[0].Sp < medium) {
		free(MediumPrime); MediumPrime = NULL;
	}
	if (MediumPrime == NULL) {
		const uint pix = (uint)(PI(medium, 1.2));
		MediumPrime = (SievePrime*)malloc(sizeof(SievePrime) * (pix + PI_65536));
		MediumPrime[0].Sp = medium;
		MediumPrime[0].Si = pix + PI_65536;
		MediumPrime[Threshold.L1Index + 0].Sp = MediumPrime[Threshold.L1Index + 1].Sp = UINT_MAX;
		MediumPrime[pix + PI_65536 - 2].Sp = MediumPrime[pix + PI_65536 - 1].Sp = UINT_MAX;
	}
}

uint64 doSieve(const uint64 start, const uint64 end, PrimeCall * pcall)
{
	const int64 ts = getTime();
	const uint sqrtp = isqrt(end);
	const uint sieve_byte = adjustConfig(Config.SieveSize, sqrtp);

	const uint medium = setBucketStart(start, sqrtp, sieve_byte * WHEEL30);
	allocMedium(medium);

	const uint l1_maxp = Threshold.L1Maxp;
	const uint max_cache = sieve_byte + Threshold.L2Size + Threshold.L1Maxp + sizeof(uint64);
	uchar* bitarray = (uchar*)malloc(max_cache); //TODO align64
	assert((size_t)bitarray % sizeof(stype) == 0);

	//init medium sieve
	const bool bmsieve = sqrtp >= l1_maxp;
	if (bmsieve) {
		initSmallWheel(l1_maxp, isqrt(medium) + 10);
		initMediumWheel(bitarray, sieve_byte * WHEEL30, medium, start - start % WHEEL210);
	}

	//init big sieve
	uint64 buckets = Threshold.BucketStart;
	if (sqrtp > medium) {
		setBucketInfo(sieve_byte / Config.Bsegs * WHEEL30, sqrtp, end - buckets);
		initSmallWheel(medium, isqrt(sqrtp) + 1);
		memset(bitarray + (L2_DCACHE_SIZE << 10), 0, l1_maxp + sizeof(uint64));
		initBigWheel(bitarray, medium, sqrtp, buckets, end - buckets);
		if (BucketInfo.CurStock < BucketInfo.LoopSize)
			allocWheelBlock(1);
	}
	else
		Threshold.BucketStart = end + 1;

	//init small sieve
	initSmallWheel(start, l1_maxp);
	memset(bitarray + sieve_byte, 0, max_cache - sieve_byte);

	const int64 ti = getTime();
	const uint64 primes = pi(bitarray, start, end, pcall);

	if (BucketInfo.StockSize > MAX_STOCK / 5) {
		for (uint i = 0; i < BucketInfo.PoolSize; i++) free(WheelPool[i]);
		memset(&BucketInfo, 0, sizeof(BucketInfo));
	}

	if (Config.Flag & PRINT_RET) {
		const int64 ta = getTime();
		printResult(start, end, primes);
		if (Config.Flag & PRINT_TIME)
			printf(" (%.2f + %.2f = %.2f sec %u kb L%u%u%u %d)",
				(ti - ts) / 1000.0, (ta - ti) / 1000.0, (ta - ts) / 1000.0, Config.SieveSize >> 10,
				Config.L1Segs, Config.L2Segs, Config.Msegs, Config.Bsegs);
		putchar('\n');
	}

#ifndef B_R
	assert(BucketInfo.StockSize == BucketInfo.CurStock);
#endif

	free(bitarray);
	return primes;
}

#if X86_64 || X86
static void cpuidInfo(int regs[4], int id, int ext)
{
#if _MSC_VER >= 1600 //2010
	__cpuidex(regs, id, ext);
#elif __GNUC__ || __TINYC__ || __clang__
	__asm__(
		"cpuid\n"
		: "=a"(regs[0]), "=b"(regs[1]), "=c"(regs[2]), "=d"(regs[3])
		: "a"(id), "c"(ext)
	);
#elif ASM_X86
	__asm
	{
		mov eax, id
		mov ecx, ext
		cpuid
		mov edi, regs
		mov dword ptr[edi + 0], eax
		mov dword ptr[edi + 4], ebx
		mov dword ptr[edi + 8], ecx
		mov dword ptr[edi + 12], edx
	}
#endif
}

static int getIntelCache(int* l3size, int* l2size, int* l1size)
{
	int regs[4];

	int i = 0;
	do {
		cpuidInfo(regs, 4, i++);
		int cache_type = regs[0] & 0x0F;
		if (cache_type == 0) {
			break;
		}

		int cache_level = ((regs[0] >> 5) & 0x7);
		int lsize = (regs[1] & 0xFFF) + 1;
		int ways = ((regs[1] >> 22) & 0x3FF) + 1;
		int partitions = ((regs[1] >> 12) & 0x3FF) + 1;
		int sets = regs[2] + 1;
		int cache_size = (ways * partitions * lsize * sets) >> 10;

		if (cache_type == 3 && cache_level == 3)
			*l3size = cache_size;
		if (cache_type == 3 && cache_level == 2)
			*l2size = cache_size;
		else if (cache_type == 1 && cache_level == 1) //data cache
			*l1size = cache_size;
		printf("\n cache_level %d, cache_type = %d cache_size = %d kb", cache_level, cache_type, cache_size);
	} while (i < 16);

	return 0;
}

static int getCpuInfo()
{
	int regs[4] = { 0 };
	char vendor[0x40] = { 0 };
	int (*pTmp)[4] = (int(*)[4])vendor;
	cpuidInfo(*pTmp++, 0x80000002, 0);
	cpuidInfo(*pTmp++, 0x80000003, 0);
	cpuidInfo(*pTmp++, 0x80000004, 0);

	for (int i = 0; vendor[i]; i++) {
		if (vendor[i] != ' ' || vendor[i + 1] != ' ')
			putchar(vendor[i]);
	}

	cpuidInfo(regs, 1, 0); int cpuCores = (regs[1] >> 16) & 0xff; // EBX[23:16]

	int l2cache[4]; cpuidInfo(l2cache, 0x80000006, 0);
	int l3Size = 0, l1Size = 64, l2Size = ((uint)l2cache[2]) >> 16;

	if (strstr(vendor, "AMD")) {
		l3Size = ((uint)l2cache[3] >> 18) * 512; //edx
		cpuidInfo(regs, 0x80000005, 0); l1Size = (uint)regs[2] >> 24;
		cpuidInfo(regs, 0x80000008, 0), cpuCores = (regs[2] & 0xFF) + 1;
		setSieveSize(l3Size / 4);
	}
	else if (strstr(vendor, "Intel")) {
		getIntelCache(&l3Size, &l2Size, &l1Size);
		cpuidInfo(regs, 0, 0);
		if (regs[0] >= 0xB) { //max_leaf
			uint ebx[2] = { 0 }, words[4] = { 0 }, trys = 0;
			cpuidInfo((int*)words, 0xB, 0);
			while (words[0] != 0 || words[1] != 0) {
				//printf("\ntrys = %d, %d %d\n", trys, words[1] & 0xFFFF, words[1]);
				if (trys < 2) ebx[trys] = words[1];
				cpuidInfo((int*)words, 0xB, ++trys);
			}
			const int ht = (ebx[0] & 0xFFFF), tc = (ebx[1] & 0xFFFF);
			cpuCores = tc;
		}
		else {
			cpuidInfo(regs, 4, 0);
			if ((regs[0] & 0x1f) != 0)
				cpuCores = ((regs[0] >> 26) & 0x3f) + 1; // EAX[31:26] + 1
			else
				cpuCores /= 2;
		}
		if (cpuCores == 0)
			cpuCores = 2;
		setSieveSize((l3Size >> 10) / 2);
	}

	setCacheSize(1, l1Size);
	setCacheSize(2, l2Size);

	printf(" L1Dsize/L2Size/L3Size = %d/%d/%d kb, cores = %u\n\n", l1Size, l2Size, l3Size, cpuCores);
	return 0;
}
#endif

void initPrime()
{
	srand((uint)time(0));
	if (SmallPrime[1].Sp == 0) {
		eratoSimple();
#if (X86_64 || X86) && !(__clang__ && __llvm__ == 0)
		getCpuInfo();
#endif
		initBitTable();
		initWheel30();
		initWheel210();
		setL1Index();
	}
}

static void fixRangeTest(uint64 lowerBound, const int64 range, uint64 Ret)
{
	const int llog10 = ilog(lowerBound, 10), rlog10 = ilog(range, 10);
	const uint64 maxrange = mpow(10, 10);
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
		setCacheSegs(1, rand() % 5 + 1), setCacheSegs(2, rand() % 5 + 2), setCacheSegs(3, rand() % 8 + 2);
		setCacheSegs(4, rand() % 8 + 1);

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
		98169972, // pi(10^19, 2^32)
		0,
	};

	int64 ts = getTime();
	Config.Flag &= ~PRINT_RET;
	uint64 primes = 0;
	Config.Progress = 0;
	for (int i = 1; i <= 10; i++) {
		primes = doSieve(0, mpow(10, i), NULL);
		printf("pi(10^%2d) = %llu\n", i, primes);
	}

	for (int j = 12; primeCounts[j]; j++) {
		uint64 start = mpow(10, j), end = start + mpow(2, 30);
		setCacheSegs(1, rand() % 6 + 2), setCacheSegs(2, rand() % 6 + 2); setCacheSegs(3, rand() % 6 + 2);
		primes = doSieve(start, end, NULL);
		if (primes == primeCounts[j])
			printf("pi(10^%d, 10^%d+2^32) = %llu                 \n", j, j, primes);
	}

	Config.Progress = 0;
	printf("Time elapsed %.f sec\n\n", (getTime() - ts) / 1000.0);
	puts("All Big tests passed SUCCESSFULLY!\nStart Rand Test");

	const uint64 pow11 = mpow(10, 11);
	const uint64 pow12 = pow11 * 10, pow9 = pow11 / 100;

	const uint64 rangeData[][3] =
	{
		{mpow(01, 11), pow11, 4118054813ul},
		{mpow(01, 12), pow12, pow9 * 37 + 607912018},
		{mpow(10, 17), pow11, 2554712095ul},
		{mpow(10, 12), pow11, 3612791400ul},
		{mpow(10, 13), pow11, 3340141707ul},
		{mpow(10, 14), pow12, pow9 * 31 + 16203073},
		{mpow(10, 15), pow12, pow9 * 28 + 952450479},
		{mpow(10, 16), pow12, pow9 * 27 + 143405794},
		{mpow(10, 18), pow12, pow9 * 24 + 127637783},
		{mpow(10, 19), pow12, pow9 * 22 + 857444126},
		{mpow(10, 19), pow12, pow9 * 22 + 857444126},
		{-1 - pow11,   pow11, 2254197466ul},
		{-1 - pow12,   pow12, pow9 * 22 + 542106206},
		//		{-1 - pow12*10, pow12*10, pow9 * 225 + 420940155},
	};

	for (uint k = 0; k < sizeof(rangeData) / sizeof(rangeData[0]); k++)
		fixRangeTest(rangeData[k][0], rangeData[k][1], rangeData[k][2]);

	Config.Flag |= PRINT_RET;
}

static void printInfo()
{
	const char* sepator =
		"------------------------------------------------------------------------------------------------------------";
	puts(sepator);
	puts("Fast implementation of the segmented sieve of Eratosthenes n < 2^64\n"
		"Copyright (C) by 2010-2025 Huang Yuanbing bailuzhou at 163.com\n"
		"Compile: g++ -march=native -funroll-loops -O3 -pipe PrimeNumber.cpp -o prime\n");

	char buff[500] = { 0 };
	char* info = buff;
#ifdef __clang__
	info += sprintf(info, "clang %s", __clang_version__); //vc/gcc/llvm
#if __llvm__
	info += sprintf(info, " on llvm/");
#endif
#endif

#if _MSC_VER
	info += sprintf(info, "Compiled by vc++ %d", _MSC_VER);
#elif __GNUC__
	info += sprintf(info, "Compiled by gcc %d.%d.%d", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#elif __INTEL_COMPILER
	info += sprintf(info, "Compiled by intel c++ %d", __INTEL_COMPILER);
#elif __TINYC__
	info += sprintf(info, "Compiled by tcc %d", __TINYC__);
#endif

#if __cplusplus
	info += sprintf(info, " __cplusplus = %d", (int)__cplusplus);
#endif

#if X86_64
	info += sprintf(info, " x86-64");
#elif X86
	info += sprintf(info, " x86");
#elif __arm64__ || __aarch64__
	info += sprintf(info, " arm64");
#elif __arm__
	info += sprintf(info, " arm");
#else
	info += sprintf(info, " unknow");
#endif

	info += sprintf(info, " %s %s in %s\n", __TIME__, __DATE__, __FILE__);
	info += sprintf(info, "[MARCO] : MEM_WHEEL = %d mb, WHEEL_SIZE = %d kb, SIEVE_SIZE = %d kb, WHEEL/WH210 = %d/%d\n",
		MEM_BLOCK * WHEEL_SIZE >> 17, WHEEL_SIZE >> 7, SIEVE_SIZE, WHEEL, WHEEL210);
	info += sprintf(info, "[CACHE] : L1Size = %u, L2Size = %u, SieveSize = %u, Bucket = %u, Block = %u\n",
		Threshold.L1Size >> 10, Threshold.L2Size >> 10, Config.SieveSize >> 10, BucketInfo.LoopSize, BucketInfo.StockSize);
	info += sprintf(info, "[ARGS ] : L1Segs/L2Segs/Mseg/Bsegs = (%u,%u,%u,%u)\n",
		Config.L1Segs, Config.L2Segs, Config.Msegs, Config.Bsegs);
	info += sprintf(info, "[ARGS ] : L1Maxp/L2Maxp/Medium/Large/SieveSize = (%u,%u,%u,%u,%u)",
		Threshold.L1Maxp, Threshold.L2Maxp, Threshold.Medium, isqrt(Threshold.BucketStart + 1), Config.SieveSize);
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

	if (++end == 0) end--; //overflow if end = 2^64-1
	initSmallWheel(start, isqrt(end) + 1);
	const uint sieve_size = (L2_DCACHE_SIZE * 1024);
	const uint max_cache = sieve_size + sizeof(uint64) + Threshold.L1Maxp;
	uchar* bitarray = (uchar*)malloc(max_cache); //TODO align64
	memset(bitarray + sieve_size, 0, max_cache - sieve_size);

	int64 ts = getTime();
	for (uint l2_size = sieve_size * WHEEL30, beg = start; beg < end; beg += l2_size) {
		if (l2_size > end - beg)
			l2_size = end - beg;
		const int bytes = segmentedSieve2(bitarray, beg, l2_size, primes > 4);
		if (align210 > 0) { //why brother ?
			memset(bitarray, ~0, align210 / WHEEL30);
			bitarray[align210 / WHEEL30] |= (1 << WheelInit30[align210 % WHEEL30].WheelIndex) - 1;
			align210 = 0;
		}
		primes += segmentProcessed(bitarray, beg, bytes, NULL);
	}

	if (Config.Flag & PRINT_RET) {
		printResult(sstart, end - 1, primes);
		printf(" (%.2f sec, %u kb)\n", (getTime() - ts) / 1000.0, sieve_size >> 10);
	}
	else {
		putchar('\r');
	}

	return primes;
}

//get the first digit number index
static int parseCmd(char params[][60])
{
	int cmdi = -1, cdata = 0;

	for (int i = 0; params[i][0]; i++) {
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
		case 'H': puts(Help); if (n) puts(Benchmark);   break;
		case 'S': setSieveSize(cdata);                  break;
		case 'C': setCacheSize(cdata % 10, cdata / 10); break;
		case 'L': setCacheSegs(cdata % 10, cdata / 10); break;
		case 'M': Config.Progress = (1 << cdata) - 1;   break;
		case 'D': Config.Flag ^= 1 << (n - 'A');        break;
		case 'I': printInfo();                          break;
		default: cmdi = i;                             break;
		}
	}

	return cmdi;
}

//split ccmd string to params array
static int splitCmd(const char* ccmd, char params[][60])
{
	int nwords = 0;

	for (int i = 0; i < 256; i++) {
		while (isspace(*ccmd) || ',' == *ccmd)
			ccmd++;

		char* pc = params[i];
		unsigned char c = *ccmd;
		if (c == 0 || c == ';')
			break;

		bool isvalid = false;
		while (isalnum(c) || c == '^' || c == '/' || c == '+' || c == '-' || c == '*') {
			*pc++ = c;
			c = *++ccmd;
			isvalid = true;
		}
		if (isvalid)
			nwords++;
		else
			ccmd++;
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
		char params[10][60] = { 0 };
		char* pcmd = (char*)strchr(cmd, ';');
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
			puts(Benchmark);
			if (isdigit(params[cmdi + 2][0])) {
				int powi = atoi(params[cmdi + 2]);
				uint64 range = powi > 12 ? mpow(2, powi) : mpow(10, powi);
				for (int i = 32; i < 64 && powi > 0; i++) {
					uint64 start2 = mpow(2, i);
					doSieve(start2, start2 + range, NULL);
				}
			}
			if (isdigit(params[cmdi + 1][0])) {
				int64 range = atoint64(params[cmdi + 1]);
				if (range < 64)
					range = range > 12 ? mpow(2, range) : mpow(10, range);
				for (int j = 11; j < 20; j++) {
					uint64 start2 = mpow(10, j);
					doSieve(start2, start2 + range, NULL);
				}
			}
			else
				startTest(0);
		}
		else if (cmdc == 'P') {
			PrimeCall pcall = { NULL, dumpPrime, 0 };
			pcall.data = &pcall.primes;
			doSieve(start, end, &pcall);
		}
		else if (cmdi >= 0) {
			doSieve(start, end, NULL);
#if FCHECK
			if (end >> 32 == 0 || sizeof(uint64) == sizeof(uint))
				pi2(start, end);
#endif
		}

		if (pcmd)
			cmd = pcmd + 1;
		else
			break;
	}

	return true;
}

static void randTest()
{
	Config.Progress = 0;
	Config.Flag ^= PRINT_RET;
	//	Config.Flag &= ~PRINT_TIME;
#if 1
	for (int i = 30; i > 0; i--)
		for (int j = 1; j < 10; j++) {
			setCacheSize(1, 32 * (rand() % 1 + 1)); setCacheSize(2, 256 * (rand() % 2 + 1));
			uint64 beg = ((uint64)(rand() * rand()) << i) + (uint64)rand() * rand() * rand();
			uint sieve_size = setSieveSize(256 << (rand() % 5));
			uint64 range = (rand() % 36) * (uint64)sieve_size * WHEEL30 + rand();
			uint64 rm1 = doSieve(beg, beg + range, NULL);
			if (setSieveSize(L2_DCACHE_SIZE << (rand() % 5 + 0)) == sieve_size)
				setSieveSize(sieve_size * 2);

			Config.Flag ^= PRINT_RET;
			setCacheSegs(3, rand() % 12 + 2); setCacheSegs(2, rand() % 6 + 2), setCacheSegs(1, rand() % 5 + 1);
			setCacheSegs(4, rand() % 8 + 1);
			uint64 rm2 = doSieve(beg, beg + range, NULL);
			if (rm1 != rm2) { printInfo(); printf("%llu %llu\n", beg, range); system("pause"); }
			if ((i + j) % 10 == 0) printf("\r %2d progress ~= %d\n", i, j);
			Config.Flag ^= PRINT_RET;
		}
#endif

	//	Config.Flag &= ~PRINT_TIME;
	//	Config.Flag ^= PRINT_RET;

	for (int j = 1; j <= 10000; j++) {
#if 1
		char cmd[100] = { 0 };
		setCacheSize(1, 32 * (rand() % 1 + 1)); setCacheSize(2, L2_DCACHE_SIZE * (rand() % 2 + 1));

		setCacheSegs(1, rand() % 6 + 1);
		setCacheSegs(2, rand() % 6 + 1);
		setCacheSegs(3, rand() % 12 + 2);
		uint64 beg = (uint64)(rand() * rand()) * (uint64)(rand() * rand()) % mpow(10, 12);
		setSieveSize(L2_DCACHE_SIZE * rand() % 8 + L2_DCACHE_SIZE);
		beg -= beg % 2;
		uint64 range = (rand() % 64 + 2) * Config.SieveSize * WHEEL30;
		const uint64 rm1 = doSieve(beg, beg + range - 1, NULL);
		sprintf(cmd, "%llu %llu s%u L%u1 L%u2 L%u3 c%d2",
			beg, range, Config.SieveSize >> 10, Config.L1Segs, Config.L2Segs, Config.Msegs, Threshold.L2Size >> 10);
		//		Config.Flag ^= PRINT_RET;

		setCacheSegs(1, rand() % 8 + 1); setCacheSegs(2, rand() % 6 + 2); setCacheSegs(4, rand() % 8 + 1);
		setSieveSize(L2_DCACHE_SIZE << ((rand() % 5) + 0));

		const uint64 rm2 = doSieve(beg, beg + range, NULL);
		if (rm1 != rm2) {
			printf("%llu != %llu --- %llu %llu s%u L%u1 L%u2 L%u3 c%d2 | ", rm1, rm2,
				beg, range, Config.SieveSize >> 10, Config.L1Segs, Config.L2Segs, Config.Msegs, Threshold.L2Size >> 10);
			puts(cmd);
		}
#endif
		if (j % 10 == 0) printf("\rprogress ~= %d %d%%\n", j, 100 * j / 10000), fflush(stdout);
	}

#if 1
	for (int k = 1; k <= 1000; k++) {
		for (int i = 6; i <= 12; i++) {
			setSieveSize(L2_DCACHE_SIZE << ((rand() % 5) + 0));
			uint sieve_size = Config.SieveSize * WHEEL30;
			uint64 medium = sieve_size / i + 1;
			medium -= medium % WHEEL210;
			uint64 start = medium * medium - sieve_size * (rand() % 8 + 1) + rand();
			uint64 end = start + sieve_size * (rand() % 32 + 1);

			setCacheSegs(3, i); setCacheSegs(2, rand() % 6 + 2), setCacheSegs(1, rand() % 6 + 1);
			const uint64 r1 = doSieve(start, end, NULL);
			char cmd2[100] = { 0 };
			sprintf(cmd2, "  %llu %llu s%u L%u1 L%u2 L%u3",
				start, end, Config.SieveSize >> 10, Config.L1Segs, Config.L2Segs, Config.Msegs);

			setSieveSize(L2_DCACHE_SIZE << ((rand() % 5) + 0));
			setCacheSegs(3, rand() % 10 + 2); setCacheSegs(2, rand() % 6 + 2), setCacheSegs(1, rand() % 5 + 1);
			const uint64 r2 = doSieve(start, end, NULL);
			if (r1 != r2) {
				printf("   %llu != %llu, %llu %llu s%u L%u1 L%u2 L%u3\n  %s\n",
					r1, r2, start, end, Config.SieveSize >> 10, Config.L1Segs, Config.L2Segs, Config.Msegs, cmd2);
			}
		}
		if (k % 10 == 0) printf("\rmedium progress ~= %d %d%%\n", k, 100 * k / 1000), fflush(stdout);
	}
#endif
}

int main(int argc, char* argv[])
{
	initPrime(); //5ms
	if (argc == 2) {
		executeCmd(argv[1]);
		return 0;
	}
#ifdef RT
	else if (argc > 2)
		randTest();
#endif

#ifdef B_R
	executeCmd("e16 e10");
#elif GCOV
	executeCmd("e18 e9");
#else
	//	if (Threshold.L2Size >= 512 << 10 || Config.SieveSize >= 4096 << 10)
	//		executeCmd("L41 L42 2^31;");
	executeCmd("e10;1e12 1e10; e14 e10;");
	executeCmd("10^12 1e9; e16 e9; e18 e9*1; i");
	executeCmd("e14 e9; e16 e10 s4");
#endif

	while (true) {
		char ccmd[257];
		printf(">> ");
		if (!fgets(ccmd, 255, stdin) || !executeCmd(ccmd))
			break;
	}

	return 0;
}

//TODO:
//1.multi-threading
//2.redesign by c++17
// cl /O2 /Oi /GS /GL  /D SIEVE_SIZE=4096 /D L2_DCACHE_SIZE=512 /D X86_64 PrimeNumber.cpp -o p98_vs
