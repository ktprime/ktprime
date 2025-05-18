/***
the most fast segmented sieving of twin prime before 2019
doc:
	http://sweet.ua.pt/tos/software/prime_sieve.html
	http://primesieve.org/
***/

const char* Benchmark =
"Mingw: g++ 5.1.0 bailuzhou@163\n"
":g++ -DSIEVE_SIZE=2048 -march=native -funroll-loops -O3 -s -pipe\n"
"Windows  10 x64  on x86_64   i3-350M, i5-3470, i7-7500u, i7-6700, R7-1700\n"
"pi2(1e11, 1e10) = 20498568    3.75     1.90     1.80      1.56     1.84\n"
"pi2(1e12, 1e10) = 17278660    4.30     2.23     2.10      1.80     2.11\n"
"pi2(1e13, 1e10) = 14735239    5.10     2.60     2.26      2.06     2.41\n"
"pi2(1e14, 1e10) = 12706059    6.20     3.22     3.00      2.47     2.83\n"
"pi2(1e15, 1e10) = 11069670    7.50     3.81     3.43      3.02     3.44\n"
"pi2(1e16, 1e10) = 9728024     8.70     4.45     4.22      3.47     3.95\n"
"pi2(1e17, 1e10) = 8614943     10.3     5.25     4.71      4.05     4.70\n"
"pi2(1e18, 1e10) = 7687050     13.4     6.45     5.55      4.97     5.86\n"
"pi2(1e19, 1e10) = 6895846     20.2     9.41     8.20      7.32     8.83\n"
"pi2(0-1e9,2^64) = 670362      7.02     3.90     3.72      3.02     3.56\n"
"pi2(1e18, 1e6)  = 794         0.75     0.46     0.37      0.50     0.58\n"
"pi2(1e18, 1e8)  = 77036       1.30     0.78     0.70      0.60     0.70\n"
"pi2(1e18, 1e9)  = 769103      3.20     1.52     1.40      1.22     1.36\n"
"pi2(1e14, 1e12) = 1270127074  655      324      244       245      280\n"
"pi2(1e16, 1e12) = 972773783   920      430      362       330      390\n"
"pi2(1e18, 1e12) = 768599834   1120     550      492       411      500\n"
"pi2(1e19, 1e12) = 689816098   1200     600      522       480      540\n";

static const char* Help = "\
	[B: Benchmark (0 - 12, 0 - 40)]\n\
	[D: D[T, R] dump time and result]\n\
	[M: Progress of calculating (0 - 20)]\n\
	[C: Cpu L1/L2 data cache size (L1:16-64, L2:256-4096)k]\n\
	[S: Set the sieve size (32 - 4096)k]\n\
	[L: Set sieve cache segs L(2-12)1, L(2-8)2 L(2-6)3, L(1-32)4]\n\
	[I: Info of programming]\n\
	[Z: Compile self and run]\n\
	[P: Print prime in [start, end]]\n\n\
Example:\n\
	1e16 1e16+10^10 s1024\n\
	p 1e12+100 100";

#include <ctype.h>
#include <math.h>
#include <memory.h>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#if __x86_64__ || __amd64__ || _M_X64 || __amd64 || __x86_64
	# define X86_64       1
#elif __i386__ || _M_IX86 || _X86_ || __i386
	# define X86          1
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

#if X86_64 && _MSC_VER
	# define ASM_X86      0
	# define BIT_SCANF    1
#elif X86_64 || X86
	#ifndef __clang__
	# define ASM_X86      1
	# define BIT_SCANF    1
	#endif
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

#if 0
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
	ERAT_SMALL    = 6, //4 - 16
	ERAT_MEDIUM   = 6, //2 - 6
	ERAT_BIG      = 6, //2 - 6

#ifndef CHAR_BIT
# define CHAR_BIT   8
#endif

	WHEEL30       = 30,
	WHEEL210      = 210,
	PRIME_PRODUCT = 210 * 11 * 13 * 17 * 19,
	FIRST_INDEX   = PRIME_PRODUCT / 9699690 + 7,
	NEXT_MULTIPLE = 0x799b799b,
	MAX_WHEEL_GAP = 28,
	PI_65536      = 6542 + 1, //pi(2^16) + 1
	MAX_SEGMENT   = 8 << 10,  //4Mkb
#ifndef UINT_MAX
	UINT_MAX      = 0-1u,
#endif

	SP_BIT        = 6, //48 < 2^SP_BIT < WHEEL210
	SI_BIT        = 8, //[8 - 10]

#if (L2_DCACHE_SIZE < 256 || L2_DCACHE_SIZE > 4096)
	L2_DCACHE_SIZE= 256,
#endif
#ifndef SIEVE_SIZE
	SIEVE_SIZE    = 4096
#endif
};

enum ECMD
{
	COUNT_PRIME = 0,
	COPY_BITS,
	SAVE_PRIME,
	SAVE_BYTE,
	PCALL_BACK
};

enum EBUCKET
{
	UINT_PIMAX = 203280221, //= pi(2^32)
	MAX_BUCKET = 5464 * 3, //= 28*2^32 / 2^18 * 30 + 4
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
	PRINT_RET  = 1 << ('R' - 'A')
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
	uint L1Segs;

	uint L2Maxp;
	uint L2Index;
	uint L2Segs;

	uint Medium;

	uint64 BucketStart; //min bucket start
};

struct Config_
{
	uint SieveSize;
	uint Progress;
	uint Flag;
	uint Msegs;
	uint Bsegs;
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
	uchar PrimeIndex;
	uchar Reserved; //remove slow ?
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

//api
struct PrimeCall
{
	int Oper;
	uchar* Data;
	uint64 Primes;
};

typedef void (*sieve_call)(uint64, uint64);
static void printInfo();
void initPrime(int sieve_size);
uint64 doSieve(const uint64 start, const uint64 end, PrimeCall* cmd);
uint setSieveSize(uint sieve_size);
void setCacheSegs(uint level, uint cachesegs);
void setCacheSize(uint level, uint cachecpu);

typedef WheelElement WheelFirst;
static WheelFirst WheelFirst30[WHEEL30][8];
static WheelFirst WheelFirst210[WHEEL210][48];

static WheelInit WheelInit30[WHEEL30];
static WheelInit WheelInit210[WHEEL210];

static WheelElement WheelData30[8][8];
static WheelElement WheelData210[48][32];
//presieved with prime <= 19
static uchar PreSieved[PRIME_PRODUCT / WHEEL30];

#if _MSC_VER >= 1400
# include <intrin.h>
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
static Stock StockCache [MAX_STOCK];

//big/bucket wheel pool
static Bucket_ Bucket [MAX_BUCKET];
static SievePrime* WheelPool [MAX_POOL];

//medium prime
static SievePrime* MediumPrime;
static SievePrime SmallPrime[PI_65536];

//small prime
static uint Prime[PI_65536];

static struct Threshold_ Threshold =
{
	32 << 10, L2_DCACHE_SIZE << 10,
	(32 << 10) / ERAT_SMALL, 0, ERAT_SMALL,
	(L2_DCACHE_SIZE << 10) / ERAT_MEDIUM, 0, ERAT_MEDIUM,
	SIEVE_SIZE * (WHEEL30 << 10) / ERAT_BIG, ERAT_BIG,
};

static struct Config_ Config =
{
	SIEVE_SIZE << 10,
	(1 << 6) - 1,
	PRINT_RET | PRINT_TIME,
	ERAT_BIG,
	8,
};

#ifndef PGAP
#define PGAP        2
#endif

#if BIT_SCANF == 0
static uchar Lsb[1 << 16];
#endif

//number of bits 1 binary representation table in Range[0-2^16)
static uchar WordNumBit1[1 << 16];

static uchar TwBits[6] =
{
#if PGAP == 2
	1, 3, 6, 9, 11, 14 //PGAP = 2
#else
	0, 2, 4, 8, 10, 12 //PGAP = 4
#endif
};

static uchar Pattern30[64] =
{
	1, 7, 11, 13, 17, 19, 23, 29,
	31, 37, 41, 43, 47, 49, 53, 59, 61
};

static int Pattern210[WHEEL210];

//the simple sieve of Eratosthenes implementated by bit packing
static int eratoSimple()
{
	int primes = 1;
	const uint maxp = (1 << 16) + 1;
	uchar bitarray[1 << 12] = {0};

	for (uint p = 3; p < maxp; p += 2) {
		if (0 == (bitarray[p >> 4] & (1 << (p / 2 & 7)))) {
			Prime[primes ++] = p;
			for (uint j = p * (p / 2) + p / 2; j <= maxp / 2; j += p)
				bitarray[j >> 3] |= 1 << (j & 7);
		}
	}
	Prime[primes ++] = maxp;
	return primes;
}

//The first presieved template, cross off the first 8th prime multiples
static void initPreSieved()
{
	for (int i = 3; PRIME_PRODUCT % Prime[i] == 0; i ++) {
		int p = Prime[i];
		for (int offset = p; offset < sizeof(PreSieved) * WHEEL30; offset += p * 2) {
			PreSieved[offset / WHEEL30] |= WheelInit30[offset % WHEEL30].MaskBit;
		}
	}
}

static void initBitTable( )
{
	int i = 0, n = 0;
	for (i = 2; Pattern30[i]; i++) {
		if (Pattern30[i] - Pattern30[i - 1] == PGAP) {
			TwBits[n++] = i - 2;
		}
	}

	for (i = 0; i < sizeof(WordNumBit1) / sizeof (WordNumBit1[0]); i ++) {
		for (int j = 0; j < sizeof(TwBits) / sizeof(TwBits[0]); j++)
			WordNumBit1[i] += (i >> TwBits[j]) % 4 == 0;
	}

	for (i = 8; i < sizeof(Pattern30) / sizeof (Pattern30[0]); i ++) {
		Pattern30[i] = Pattern30[i - 8] + WHEEL30;
	}

#if BIT_SCANF == 0
	int nbitsize2 = sizeof(Lsb) / sizeof(Lsb[0]);
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
		WheelInit30[j].MaskBit = 0;
		WheelInit30[j].WheelIndex = wj;
		WheelInit30[j].PrimeIndex = wj;
		if (j == Pattern30[wj]) {
			WheelInit30[j].MaskBit = 1 << (wj ++);
		}
	}

	for (int i = 0; i < WHEEL30; i ++) {
		for (int pi = 0; pi < 8; pi ++) {
			int multiples = 0, offset = i;
			if (i % 2 == 0) {
				multiples = 1;
				offset += Pattern30[pi];
			}
			while (WheelInit30[offset % WHEEL30].MaskBit == 0) {
				offset += Pattern30[pi] * 2;
				multiples += 2;
			}
			int wi = WheelInit30[offset % WHEEL30].WheelIndex;
			WheelElement& wf = WheelFirst30[i][pi];
			wf.Multiple = (nextMultiple[wi] >> (pi * 4)) & 15;
			wf.WheelIndex = wi;
			wf.Correct = multiples;
			wf.MaskBit = 1 << wi;
		}
	}

	for (int wi = 0; wi < 8; wi ++) {
		for (int pi = 0; pi < 8; pi ++) {
			int multiples = 2;
			int next = Pattern30[wi] + Pattern30[pi] * 2;
			while (WheelInit30[next % WHEEL30].MaskBit == 0) {
				next += Pattern30[pi] * 2;
				multiples += 2;
			}

			WheelElement& we30 = WheelData30[pi][wi];
			we30.Multiple = multiples * (WHEEL30 / WHEEL30);
			we30.WheelIndex = WheelInit30[next % WHEEL30].WheelIndex;
			we30.Correct = next / WHEEL30 - Pattern30[wi] / WHEEL30;
			we30.MaskBit = WheelInit30[Pattern30[wi]].MaskBit;
		}
	}
}

static void initWheel210()
{
	int wi = 0, i = 0;
	const int psize = 48, wsize = 30;
	int wpattern[WHEEL210] = {0};

	for (i = 0; i < WHEEL210; i ++) {
		const uchar mask = WheelInit30[i % WHEEL30].MaskBit;
		WheelInit210[i].MaskBit = mask;
		WheelInit210[i].WheelIndex = -1;
		WheelInit210[i].PrimeIndex = wi;
		if (mask && i % (WHEEL210 / WHEEL30))
			Pattern210[wi ++] = i;
	}

	wi = 0, i = 1;
	if (PGAP == 2) {
		WheelInit210[i].WheelIndex = wi;
		wpattern[wi ++] = i;
	}

	for (; i < WHEEL210; i += 2) {
		const uchar mask = WheelInit210[i % WHEEL30].MaskBit;
		if (mask && i % (WHEEL210 / WHEEL30) &&
			(i + PGAP) % (WHEEL210 / WHEEL30) &&
			WheelInit210[(i + PGAP) % WHEEL30].MaskBit) {
			WheelInit210[i].WheelIndex = wi;

			wpattern[wi ++] = i;
			if (i + PGAP < WHEEL210) {
				WheelInit210[i + PGAP].WheelIndex = wi;
				wpattern[wi ++] = (i + PGAP) % WHEEL210;
			}
		}
	}

	for (wi = 0; wi < wsize; wi ++) {
		for (int pi = 0; pi < psize; pi ++) {
			int multiples = 2;
			int next = wpattern[wi] + Pattern210[pi] * 2;
			while (WheelInit210[next % WHEEL210].WheelIndex < 0) {
				next += Pattern210[pi] * 2;
				multiples += 2;
			}

			WheelElement& we210 = WheelData210[pi][wi];
			we210.Correct = next / WHEEL30 - wpattern[wi] / WHEEL30;
			we210.MaskBit = WheelInit210[wpattern[wi]].MaskBit;
			we210.WheelIndex = WheelInit210[next % WHEEL210].WheelIndex;
			we210.Multiple = multiples * (WHEEL210 / WHEEL30);
		}
	}

	for (i = 0; i < WHEEL210; i ++) {
		for (int pi = 0; pi < psize; pi ++) {
			int multiples = 0, next = i;
			if (i % 2 == 0) {
				multiples = 1;
				next += Pattern210[pi];
			}

			while (WheelInit210[next % WHEEL210].WheelIndex < 0) {
				next += Pattern210[pi] * 2;
				multiples += 2;
			}

			WheelFirst210[i][pi].WheelIndex = WheelInit210[next % WHEEL210].WheelIndex;
			WheelFirst210[i][pi].Correct = multiples;
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
	FILETIME ptime[4] = {0};
	GetThreadTimes(GetCurrentThread(), &ptime[0], &ptime[1], &ptime[2], &ptime[3]);
	return (/*((uint64)(ptime[2].dwHighDateTime + ptime[3].dwHighDateTime) << 32) +*/ ptime[2].dwLowDateTime + ptime[3].dwLowDateTime) / 10000;
	//return clock();
#elif __linux__ || __unix__
	struct rusage rup;
	getrusage(RUSAGE_SELF, &rup);
	long sec  = rup.ru_utime.tv_sec  + rup.ru_stime.tv_sec;
	long usec = rup.ru_utime.tv_usec + rup.ru_stime.tv_usec;
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
#if _MSC_VER > 1400
	unsigned long index;
	#if X86_64
	_BitScanForward64(&index, n);
	#else
	_BitScanForward(&index, n);
	#endif
#elif __GNUC__
	#if X86_64
//	uint index = __builtin_ffsll(n) - 1;
	uint index = __builtin_ctzll(n);
	#else
	uint index = __builtin_ctzl(n);
	#endif
#elif ASM_X86
	stype index;
	#if X86_64
	#if __GNUC__ || __TINYC__
	__asm__ ("bsfq %1, %0\n" : "=r" (index) : "rm" (n) : "cc");
	#else
	__asm
	{
		bsfq eax, n
		mov index, eax
	}
	#endif
	#else
	#if __GNUC__ || __TINYC__
	__asm__ ("bsfl %1, %0\n" : "=r" (index) : "rm" (n) : "cc");
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
#elif __GNUC__ || __TINYC__
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

inline static int countBitsTable(const uint64 n)
{
	const uint hig = (uint)(n >> 32), low = (uint)n;
	int sum =  WordNumBit1[(ushort)hig] + WordNumBit1[(ushort)low];
	sum += WordNumBit1[hig >> 16] + WordNumBit1[low >> 16];
	return sum;
}

static int crossOffWheelFactor2(uchar* ps, const uchar* pend, const uint p)
{

	#define poSet(s0,o0,b0,s1,o1,b1,s2,o2,b2,s3,o3,b3,s4,o4,b4,s5,o5,b5)\
	while (ps <= pend) {\
		ps[o * s0 + o0] |= BIT##b0,  ps[o * s1 + o1] |= BIT##b1;\
		ps[o * s2 + o2] |= BIT##b2, ps[o * s3 + o3] |= BIT##b3;\
		ps[o * s4 + o4] |= BIT##b4, ps[o * s5 + o5] |= BIT##b5;\
		ps += p;\
	}

	const uchar* pbeg = ps;
	const uint o = p / WHEEL30;
	switch (WheelInit30[p % WHEEL30].WheelIndex)
	{
#if PGAP == 2 //remove 7->1, 23->6
		case 3 : poSet(0,0,0, 6,2,5, 10,4,2, 16,6,7, 22,9,4, 24,10,3);	break;
		case 7 : poSet(0,0,0, 2,1,7, 12,11,5, 14,13,4, 18,17,3,20,19,2);break;
		case 2 : poSet(0,0,0, 8,2,7, 12,4,3, 18,6,5, 20,7,2, 26,9,4);	break;
		case 1 : poSet(0,0,0, 4,0,7, 6,1,3, 10,2,2,  24,5,5, 28,6,4);	break;
		case 6 : poSet(0,0,0, 2,1,4, 6,4,5, 20,15,2, 24,18,3, 26,19,7);	break;
		case 4 : poSet(0,0,0, 6,3,3, 8,4,4, 14,7,7, 20,11,2, 24,13,5);	break;
		case 0 : poSet(0,0,0, 10,0,2, 12,0,3,  16,0,4, 18,0,5, 28,0,7);	break;
		case 5 : poSet(0,0,0, 4,2,4, 10,6,2, 12,7,5, 18,11,3, 22,13,7);	break;
#elif PGAP == 4  //remove 1->0, 29->7
		case 3 : poSet(4,1,6, 6,2,5, 10,4,2, 12,5,1, 22,9,4, 24,10,3);	 break;
		case 7 : poSet(8,7,6, 12,11,5,14,13,4, 18,17,3, 20,19,2,24,23,1);break;
		case 2 : poSet(2,0,6, 6,2,1, 12,4,3, 18,6,5, 20,7,2, 26,9,4);	 break;
		case 1 : poSet(6,1,3, 10,2,2, 16,3,6, 18,4,1, 24,5,5, 28,6,4);	 break;
		case 6 : poSet(2,1,4, 6,4,5, 12,9,1, 14,10,6, 20,15,2, 24,18,3); break;
		case 4 : poSet(6,3,3, 8,4,4, 18,10,1, 20,11,2, 24,13,5, 26,14,6);break;
		case 0 : poSet(6,0,1, 10,0,2, 12,0,3, 16,0,4, 18,0,5, 22,0,6);	 break;
		case 5 : poSet(4,2,4, 10,6,2, 12,7,5, 18,11,3, 24,15,1, 28,17,6);break;
#endif
	}

	return (int)(ps - pbeg) * WHEEL30;
}

static int crossOffWheelFactor(uchar* ps, const uchar* pend, const uint p)
{
	#define poSet2(s1,o1,b1,s2,o2,b2,s3,o3,b3,s4,o4,b4,s5,o5,b5,s6,o6,b6,s7,o7,b7)\
	while (ps <= pend) {\
		ps[o * 00 + 00] |= 1 << 00, ps[o * s1 + o1] |= BIT##b1;\
		ps[o * s2 + o2] |= BIT##b2, ps[o * s3 + o3] |= BIT##b3;\
		ps[o * s4 + o4] |= BIT##b4, ps[o * s5 + o5] |= BIT##b5;\
		ps[o * s6 + o6] |= BIT##b6, ps[o * s7 + o7] |= BIT##b7;\
		ps += p;\
	}

	const uchar* pbeg = ps;
	const uint o = p / WHEEL30;
	switch (WheelInit30[p % WHEEL30].WheelIndex) //37216505
	{
		case 3: poSet2(4,1,6, 6,2,5,  10,4,2,  12,5,1,  16,6,7,  22,9,4,  24,10,3) break;
		case 7: poSet2(2,1,7, 8,7,6,  12,11,5, 14,13,4, 18,17,3, 20,19,2, 24,23,1) break;
		case 2: poSet2(2,0,6, 6,2,1,  8,2,7,   12,4,3,  18,6,5,  20,7,2,  26,9, 4) break;
		case 1: poSet2(4,0,7, 6,1,3,  10,2,2,  16,3,6,  18,4,1,  24,5,5,  28,6, 4) break;
		case 6: poSet2(2,1,4, 6,4,5,  12,9,1,  14,10,6, 20,15,2, 24,18,3, 26,19,7) break;
		case 4: poSet2(6,3,3, 8,4,4,  14,7,7,  18,10,1, 20,11,2, 24,13,5, 26,14,6) break;
		case 0: poSet2(6,0,1, 10,0,2, 12,0,3,  16,0,4,  18,0,5,  22,0,6,  28,0, 7) break;
		case 5: poSet2(4,2,4, 10,6,2, 12,7,5,  18,11,3, 22,13,7, 24,15,1, 28,17,6) break;
	}

	return (int)(ps - pbeg) * WHEEL30;
}

//count number of bit 0 in binary representation of array
//bit 0 mean it's a prime position
static int countBit0sArray(const uint64 bitarray[], const int bitsize)
{
	int bit0s = 0, loops = bitsize / 64;
	while (loops -- >= 0) {
		bit0s += countBitsTable((*bitarray >> 1) | (bitarray[1] << 63));
		bitarray ++;
	}

	return bit0s;
}

static void dumpPrime(uint64 index, uint64 prime)
{
	printf("(%llu, %llu)\n", prime, prime + PGAP);
}

static int doCall(const ushort bitarray[], uint64 offset, const int size, uint64 sum_prime, sieve_call func)
{
	int primes = 0;
	for (int bi = 0; bi <= size / 2; bi++) {
		ushort mask = (bitarray[bi] >> 1) | (bitarray[bi + 1] << 15);
		for (int i = 0; i < sizeof(TwBits); i++) {
			const uint bits = TwBits[i];
			if ((mask >> bits) % 4 == 0)
				func(++primes + sum_prime, offset + Pattern30[(1 + bits) % 8] + (1 + bits) / 8 * WHEEL30);
		}
		offset += WHEEL30 * sizeof(mask);
	}

	return primes;
}

static int savePrimeByte(const stype bitarray[], const int bytes, uchar* prime)
{
	int primes = 0, lastp = (char)prime[0];
	stype mask = ~bitarray[0];
	int size = (bytes - 1) / sizeof(mask);

	for (int bi = 0; bi <= size; ) {
		if (mask == 0) {
			mask = ~bitarray[++bi];
			continue;
		}

		const int curp = bi * WHEEL30 * sizeof(mask) + PRIME_OFFSET(mask);
		mask &= mask - 1;
		*prime++ = (curp - lastp) + (curp - lastp) / 256; //save adjective prime difference
		primes++;
		lastp = curp;
	}

	//assert(bytes * WHEEL30 - lastp) < 256
	*prime = lastp - bytes * WHEEL30;

	return primes;
}

//core memory alloc
static void allocWheelBlock(const uint blocks)
{
	SievePrime *pSprime = (SievePrime*) malloc((blocks * MEM_BLOCK + 1) * MEM_WHEEL);
	WheelPool[BucketInfo.PoolSize ++] = pSprime;
	//assert (BucketInfo.PoolSize < sizeof(WheelPool) / sizeof(WheelPool[0]));
	//assert (BucketInfo.StockSize + blocks < sizeof(StockCache) / sizeof(StockCache[0]));

	//align by MEM_WHEEL
	pSprime = (SievePrime*)((size_t)pSprime + MEM_WHEEL - (uint)((size_t)pSprime) % MEM_WHEEL);

	for (int j = 0; j < blocks; j++) {
		for (uint i = 0; i < MEM_BLOCK; i ++) {
			Stock* pStock = StockCache + i + BucketInfo.StockSize;
			pStock->Sprime = pSprime + WHEEL_SIZE * i;
			pStock->Next = pStock + 1;
		}
		StockCache[BucketInfo.StockSize + MEM_BLOCK - 1].Next = StockHead;
		StockHead = StockCache + BucketInfo.StockSize;

		BucketInfo.StockSize += MEM_BLOCK;
		BucketInfo.CurStock +=  MEM_BLOCK;
		pSprime += WHEEL_SIZE * MEM_BLOCK;
	}
}

#define	PI(x, r) (x / log((double)x) * (1 + r / log((double)x)))
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
	const uint segsize = Threshold.L1Size * WHEEL30;
	uint64 offset = start - start % WHEEL210;
	uint p = 0, j = Threshold.L1Index;

	for (; p < maxp; p = Prime[++j])
		SmallPrime[j].Si = 0;

	for (p = Prime[FIRST_INDEX], j = FIRST_INDEX; p < maxp; p = Prime[++j]) {
		const uint64 p2 = (uint64)p * p;
		if (p2 > offset && p2 > offset + segsize) //overflow
			offset += (p2 - offset) / segsize * segsize;

		uint sieve_index = p - (uint)(offset % p);
		if (p2 > offset)
			sieve_index = (uint)(p2 - offset);

		const WheelFirst& wf = WheelFirst30[sieve_index % WHEEL30][WheelInit30[p % WHEEL30].WheelIndex];
		const ushort multiples = NEXT_MULTIPLE >> (wf.Multiple * 2);
		sieve_index += wf.Correct * p;

		SmallPrime[j].Sp = multiples;
		SmallPrime[j].Si = 0 - sieve_index;
	}

}

static int segmentedSieve2(uchar bitarray[], uint start, uint sieve_size, bool bcopy);
static void setWheelMedium(uchar* bitarray, const uint sieve_size, const uint medium, const uint64 start)
{
	uint j = Threshold.L1Index, l1_maxp = Threshold.L1Maxp - Threshold.L1Maxp % WHEEL210;
	Threshold.L2Index = 0;
	uint segsize = MIN(sieve_size, Threshold.L2Size * WHEEL30);

	const uint pix = PI(medium, 1.2);
	MediumPrime = (SievePrime*) malloc(sizeof(SievePrime) * (pix + 100));
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
			segsize = sieve_size;
		}

		uint64 offset = start;
		const uint64 p2 = (uint64)p * p;
		if (p2 > offset && p2 >= offset + segsize) //overflow
			offset += (p2 - offset) / segsize * segsize;

		uint sieve_index = p - (uint)(offset % p);
		if (p2 > offset)
			sieve_index = (uint)(p2 - offset);

		const int pi = WHEEL_INIT[p % WHEEL].PrimeIndex;
		const WheelFirst& wf = WHEEL_FIRST[(sieve_index + uint(offset % WHEEL)) % WHEEL][pi];
		sieve_index += wf.Correct * p;
		//assert(sieve_index / WHEEL30 < (-1u >> SP_BIT));
		MediumPrime[j + 0].Sp = (p / WHEEL << SI_BIT) + pi;
		MediumPrime[j ++ ].Si = (sieve_index / WHEEL30 << SP_BIT) | wf.WheelIndex;
	}

	MediumPrime[j + 0].Sp = MediumPrime[j + 1].Sp = MediumPrime[j + 2].Sp = UINT_MAX;
}

static void pushBucket(const uint offset, const uint sp, const uchar wi)
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
	for (uint l2_size = L2_DCACHE_SIZE * WHEEL30 << 10; medium < sqrtp; medium += l2_size) {
		if (l2_size > sqrtp - medium)
			l2_size = sqrtp - medium;

		const int bytes = segmentedSieve2(bitarray, medium, l2_size, nextp > 0);
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
			uint sieve_index = p - fastllDiv(start - remp * p, p);
#else
			uint sieve_index = p - (uint)(start % p);
#endif

#if 1
			if (sieve_index > range) continue;
#endif
			const uint sp = (p / WHEEL210 << SP_BIT) + WheelInit210[p % WHEEL210].PrimeIndex;
			const uint modulo_210 = sieve_index % WHEEL210;
			const WheelFirst& wf = WheelFirst210[modulo_210][sp % (1 << SP_BIT)];
#if X86_64
			sieve_index = (sieve_index + wf.Correct * (uint64)p) / WHEEL30;
#else
			sieve_index = (sieve_index / WHEEL210 + (sp >> SP_BIT) * wf.Correct) * (WHEEL210 / WHEEL30);
			sieve_index += (wf.Correct * (p % WHEEL210) + modulo_210) / WHEEL30;
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
			pSprime->Si = (sieve_index & BucketInfo.ModuloSize) << SI_BIT | wf.WheelIndex;
		}
	}
}

static int crossSmall0(uchar bitarray[], const uchar* pend, const uint p, uint offset, ushort multiples)
{
	uchar* ppbeg[8], wi;
	for (int i = 0; i < 8; i ++) {
		wi = WheelInit30[offset % WHEEL30].WheelIndex;
		ppbeg[wi] = bitarray + offset / WHEEL30;
		offset += (multiples % 4) * 2 * p; multiples /= 4;
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

static int crossSmall1(uchar bitarray[], const uchar* pend, const uint p, uint offset, ushort multiples)
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

static int crossSmall3(uchar bitarray[], const uchar* pend, const uint p, uint offset, ushort multiples)
{
	for (int i = 0; i < 8; i ++) {
		const uchar mask = WheelInit30[offset % WHEEL30].MaskBit;
		if (mask == BIT0)
			break;

		bitarray[offset / WHEEL30] |= mask;
		offset += (multiples % 4) * 2 * p; multiples /= 4;
	}

	return crossOffWheelFactor2(bitarray + offset / WHEEL30, pend, p) + offset;
}

//sieve big prime from bucket
static void crossBig1(uchar bitarray[], const uint sieve_size, const SievePrime* pSprime)
{
	uint offset = pSprime->Si, sp = pSprime->Sp;
	const WheelElement* wd = WheelData210[sp % (1 << SP_BIT)];
	const WheelElement* we = wd + offset % (1 << SI_BIT);
	bitarray[offset >>= SI_BIT] |= we->MaskBit;
	offset += we->Correct + we->Multiple * (sp >> SP_BIT);

#if ERAT_BIG > 12
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

#if ERAT_BIG > 12
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
	if (bits % CHAR_BIT != 0)
		bitarray[bits / CHAR_BIT] |= ~((1 << (bits & 7)) - 1);
}

static void eratSieveL0(uchar bitarray[], const uint64 start, const uint sieve_size, uint maxp)
{
	const uchar* pend = bitarray + sieve_size / WHEEL30;
	if (start + sieve_size < ((uint64)maxp) * maxp)
		maxp = isqrt(start + sieve_size) + 1;

#ifndef FSL1
	for (uint p = Prime[FIRST_INDEX], j = FIRST_INDEX; p < maxp; p = Prime[++j]) {
		uint& offset = SmallPrime[j].Si;
		if ((int)offset < 0)
			offset = 0 - offset;
		offset += crossSmall0(bitarray, pend, p, offset, SmallPrime[j].Sp) + 30 * p - sieve_size;
	}
#else
	for (uint p = Prime[FIRST_INDEX], j = FIRST_INDEX; p < maxp; p = Prime[++j]) {
		uint& offset = SmallPrime[j].Si;
		if ((int)offset < 0) {
			offset = 0 - offset;
			for (int i = 0; i < 16; i += 2) {
				const uchar mask = WheelInit30[offset % WHEEL30].MaskBit;
				if (mask == BIT0)
					break;

				bitarray[offset / WHEEL30] |= mask;
				offset += (SmallPrime[j].Sp >> i) % 4 * p * 2;
			}
		}
		offset += crossOffWheelFactor(bitarray + offset / WHEEL30, pend, p) - sieve_size;
	}
#endif
}

static void eratSieveL1(uchar bitarray[], const uint64 start, const uint sieve_size, uint maxp)
{
	const uchar* pend = bitarray + sieve_size / WHEEL30;
	if (start + sieve_size < ((uint64)maxp) * maxp)
		maxp = isqrt(start + sieve_size) + 1;

	for (uint p = Prime[FIRST_INDEX], j = FIRST_INDEX; p < maxp; p = Prime[++j]) {
		uint& offset = SmallPrime[j].Si;
		if ((int)offset > 0)
			offset += crossOffWheelFactor2(bitarray + offset / WHEEL30, pend, p) - sieve_size;
		else
			offset = crossSmall3(bitarray, pend, p, 0 - offset, SmallPrime[j].Sp) - sieve_size;
	}
}

static void eratSieveSmall(uchar bitarray[], const uint64 start, const uint sieve_size)
{
	for (uint offset = 0, l1_size = Threshold.L1Size * WHEEL30; offset < sieve_size; offset += l1_size) {
		if (l1_size + offset > sieve_size)
			l1_size = sieve_size - offset;
		eratSieveL1(bitarray + offset / WHEEL30, start + offset, l1_size, Threshold.L1Maxp);
	}
}

#define SAFE_SET(n) \
	we##n = wd##n[(int)we##n.WheelIndex]; \
	bitarray[(int)offset##n] |= we##n.MaskBit; \
	offset##n += we##n.Correct + (we##n.Multiple) * wi##n

//sieve 1 medium prime from array
static void crossMedium1(uchar bitarray[], const uint sieve_byte, SievePrime* pSprime)
{
	uint& wheel = pSprime->Si;
	int offset = (wheel >> SP_BIT) - sieve_byte;
	if (offset >= 0) {
		wheel -= sieve_byte << SP_BIT;
		return;
	}

	const uint wi = pSprime->Sp >> SI_BIT;
	WheelElement* wd = WHEEL_DATA[pSprime->Sp % (1 << SI_BIT)];
	WheelElement we; we.WheelIndex = wheel % (1 << SP_BIT);

	do {
		we = wd[we.WheelIndex];
		bitarray[(int)offset] |= we.MaskBit;
		offset += we.Correct + (we.Multiple) * wi;
	} while (offset < 0);

	wheel = (offset << SP_BIT) | we.WheelIndex;
}

//sieve 3 medium prime from array
inline static int crossMedium3(uchar bitarray[], const uint sieve_byte, SievePrime* pSprime)
{
	uint& si0 = pSprime[0].Si, sp0 = pSprime[0].Sp;
	uint wi0 = sp0 >> SI_BIT, offset0 = (si0 >> SP_BIT) - sieve_byte;

	uint& si1 = pSprime[1].Si, sp1 = pSprime[1].Sp;
	uint wi1 = sp1 >> SI_BIT, offset1 = (si1 >> SP_BIT) - sieve_byte;

	uint& si2 = pSprime[2].Si, sp2 = pSprime[2].Sp;
	uint wi2 = sp2 >> SI_BIT, offset2 = (si2 >> SP_BIT) - sieve_byte;

	WheelElement* wd0 = WHEEL_DATA[sp0 % (1 << SI_BIT)];
	WheelElement* wd1 = WHEEL_DATA[sp1 % (1 << SI_BIT)];
	WheelElement* wd2 = WHEEL_DATA[sp2 % (1 << SI_BIT)];

	WheelElement we0, we1, we2;
	we0.WheelIndex = si0 % (1 << SP_BIT);
	we1.WheelIndex = si1 % (1 << SP_BIT);
	we2.WheelIndex = si2 % (1 << SP_BIT);

	while ((int)offset0 < 0) {
		SAFE_SET(0);
		if ((int)offset1 >= 0) break;
		SAFE_SET(1);
		if ((int)offset2 >= 2) break;
		SAFE_SET(2);
	}

	while ((int)offset0 < 0) { SAFE_SET(0); }
	while ((int)offset1 < 0) { SAFE_SET(1); }
	while ((int)offset2 < 0) { SAFE_SET(2); }

	si0 = offset0 << SP_BIT | we0.WheelIndex;
	si1 = offset1 << SP_BIT | we1.WheelIndex;
	si2 = offset2 << SP_BIT | we2.WheelIndex;

	return 3;
}

//sieve prime multiples in [start, start + sieve_size) with medium algorithm
static void eratSieveMedium(uchar bitarray[], const uint64 start, const uint sieve_size, const uint wmi, uint maxp)
{
	if (start + sieve_size < (uint64)maxp * maxp)
		maxp = isqrt(start + sieve_size) + 1;

	const uint sieve_byte = sieve_size / WHEEL30 + (sieve_size % WHEEL30 ? 1 : 0);
	SievePrime* pSprime = MediumPrime + wmi;

	const uint maxsp = (maxp / WHEEL << SI_BIT) + WHEEL_INIT[maxp % WHEEL].PrimeIndex;
	bitarray += sieve_byte;

	while (pSprime[2].Sp < maxsp) {
		pSprime += crossMedium3(bitarray, sieve_byte, pSprime);
	}
	if (pSprime[0].Sp < maxsp)
		crossMedium1(bitarray, sieve_byte, pSprime ++);
	if (pSprime[0].Sp < maxsp)
		crossMedium1(bitarray, sieve_byte, pSprime);
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
			pSprime += crossBig2(bitarray, sieve_size, pSprime);
			loops -= 2;
		}
		BucketInfo.CurStock ++;
		pStock = pnext;
	}

	BucketInfo.MaxBucket --;
	memmove(Bucket, Bucket + 1, BucketInfo.LoopSize * sizeof(Bucket[0]));
}

static int segmentProcessed(uchar bitarray[], const uint64 start, const uint bytes, PrimeCall* pcall)
{
	int primes = 0;
	uint64& lastqw = *(uint64*)(bitarray + bytes);
	const uint64 tmp = lastqw;
	lastqw = ~0;

	const int oper = pcall ? pcall->Oper : COUNT_PRIME;
	if (COUNT_PRIME == oper) {
		primes = countBit0sArray((uint64*)bitarray, bytes * sizeof(uint64));
	} else if (PCALL_BACK == oper) {
		primes = doCall((ushort*)bitarray, start, bytes, pcall->Primes, (sieve_call)pcall->Data);
	} else if (SAVE_BYTE == oper) {
		primes = savePrimeByte((stype*)bitarray, bytes, pcall->Data + pcall->Primes);
	}

	lastqw = tmp;
	if (pcall)
		pcall->Primes += primes;
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
	const uint copy_size = Threshold.L1Maxp, copy_from = Config.SieveSize;
	for (uint i = 0; i < copy_size; i += sizeof(uint)) {
		uint& c = *(uint*)(bitarray + i + copy_from);
		*(uint*)(bitarray + i) |= c, c = 0;
	}

	for (uint offset = 0, l2_size = Threshold.L2Size * WHEEL30; offset < sieve_size; offset += l2_size) {
		if (l2_size + offset > sieve_size)
			l2_size = sieve_size - offset;

		uchar* buffer = bitarray + offset / WHEEL30;
		eratSieveSmall(buffer, start + offset, l2_size);
		if (bmsieve)
			eratSieveMedium(buffer, start + offset, l2_size, Threshold.L1Index, Threshold.L2Maxp);
	}

	if (medium > Threshold.L2Maxp)
		eratSieveMedium(bitarray, start, sieve_size, Threshold.L2Index, medium);

	if (start >= Threshold.BucketStart) {
		const uint big_size = 1 << BucketInfo.Log2Size;
		for (uint offset = 0; offset < copy_from; offset += big_size)
			eratSieveBig(bitarray + offset, big_size);
	}

	return sieve_size / WHEEL30 + (7 + WheelInit30[sieve_size % WHEEL30].WheelIndex) / 8;
}

static int segmentedSieve2(uchar bitarray[], const uint start, const uint sieve_size, bool bcopy = true)
{
	const uint sqrtp = isqrt(start + sieve_size);
	const uint l1_maxp = Threshold.L1Maxp;
	preSieve(bitarray, start, sieve_size);

#if FSL1
	const uint copy_from = L2_DCACHE_SIZE << 10;
	if (sieve_size <= copy_from * WHEEL30 && bcopy) {
		for (uint i = 0; i < l1_maxp; i += sizeof(uint)) {
			uint& c = *(uint*)(bitarray + i + copy_from);
			*(uint*)(bitarray + i) |= c, c = 0;
		}
	}
#endif

	for (uint offset = 0, l1_size = Threshold.L1Size * WHEEL30; offset < sieve_size; offset += l1_size) {
		if (l1_size + offset > sieve_size)
			l1_size = sieve_size - offset;
		eratSieveL0(bitarray + offset / WHEEL30, start + offset, l1_size, l1_maxp);
	}

	const uchar* pend = bitarray + sieve_size / WHEEL30;
	for (uint j = Threshold.L1Index, p = l1_maxp; p <= sqrtp; p = Prime[++j]) {
		uint offset = p - start % p;
//		if ((int)offset <= 0) {
			offset = p - start % p;
			if (p * p > start)
				offset = p * p - start;
//		}

		const WheelFirst& wf = WheelFirst30[offset % WHEEL30][WheelInit30[p % WHEEL30].WheelIndex];
		const ushort multiples = NEXT_MULTIPLE >> (wf.Multiple * 2);
		offset = crossSmall1(bitarray, pend, p, offset + wf.Correct * p, multiples);
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
	if (Threshold.L1Maxp > 65521)
		Threshold.L1Maxp = 65521; //the bigest prime % 210 = 1 and < 2^16
	for (uint p = 0, j = 300; ; p = Prime[++j]) {
		if (p >= Threshold.L1Maxp && p % WHEEL210 == 1) {
			Threshold.L1Index = j;
			Threshold.L1Maxp = p;
			break;
		}
	}
}

void setCacheSize(int level, uint cache)
{
	cache = 1 << ilog(cache, 2);
	if (level == 1 && cache >= 16 && cache <= (Threshold.L2Size >> 10)) {
		Threshold.L1Size = cache << 10;
		Threshold.L1Maxp = Threshold.L1Size / Threshold.L1Segs;
	} else if (level == 2 && cache >= (Threshold.L1Size >> 10) && cache <= MAX_SEGMENT) {
		Threshold.L2Size = (cache << 10) / Threshold.L1Size * Threshold.L1Size;
		Threshold.L2Maxp = Threshold.L2Size / Threshold.L2Segs;
	}
}

void setCacheSegs(uint level, uint segs)
{
	if (segs == 0 || segs > 10)
		return;

	if (level == 1) {
		Threshold.L1Segs = segs;
		Threshold.L1Maxp = Threshold.L1Size / segs;
		setL1Index();
	} else if (level == 2) {
		Threshold.L2Segs = segs;
		Threshold.L2Maxp = Threshold.L2Size /segs;
	} else if (level == 3) {
		Config.Msegs = segs;
	} else if (level == 4) {
		Config.Bsegs = 1 << ilog(segs, 2);
	}
}

//sieve_size : 32k - 4M
uint setSieveSize(uint sieve_size)
{
	if (sieve_size <= MAX_SEGMENT && sieve_size > L2_DCACHE_SIZE) {
		sieve_size = (sieve_size / L2_DCACHE_SIZE) * L2_DCACHE_SIZE << 10;// 1 << (ilog(sieve_size, 2) + 10);
	} else if (sieve_size <= L2_DCACHE_SIZE && sieve_size >= 32) {
		sieve_size = (sieve_size / 32) * 32 << 10;
	} else if (sieve_size <= (MAX_SEGMENT >> 10) && sieve_size > 0) {
		sieve_size = sieve_size << 20;
	} else {
		sieve_size = SIEVE_SIZE << 10;
	}

	return Config.SieveSize = sieve_size;
}

static int checkSmall(const uint64 start, const uint64 end, PrimeCall* cmd)
{
	const uint smallPrime[] = {3, 5, 7, 0, 0, 0};
	uint primes = 0;
	for (int i = 0; i < 5; i++) {
		uint p = smallPrime[i];
		if (start <= p && (p == smallPrime[i + 1] - PGAP || p == smallPrime[i + 2] - PGAP) && smallPrime[i + 1] <= end) {
			primes ++;
			if (cmd && cmd->Oper == PCALL_BACK) {
				(*(sieve_call)cmd->Data)(primes, p);
				cmd->Primes += 1;
			}
		}
	}

	return primes;
}

static uint64 pi2(uchar* bitarray,uint64 start, uint64 end, PrimeCall* pcall)
{
	const int64 ts = getTime();
	uint64 primes = checkSmall(start, end, pcall);

	uint align210 = (uint)(start % WHEEL210);
	start -= align210;

	if (++end == 0) end --; //overflow if end = 2^64-1

	uint64 last_qword = ~0;
#ifndef _DEBUG
	double lge = end * log(end), lgs = start * log(start + 1);
	double pie = PI(end, 1.2), pis = PI(start, 1.2);
#endif

	for (uint si = 0, sieve_size = Config.SieveSize * WHEEL30; start < end; start += sieve_size) {
		if (sieve_size > end - start)
			sieve_size = (uint)(end - start);

		const uint bytes = segmentedSieve(bitarray + 8, start, sieve_size);
		if (align210 > 0) {
			memset(bitarray + 8, ~0, align210 / WHEEL30);
			bitarray[align210 / WHEEL30 + 8] |= (1 << WheelInit30[align210 % WHEEL30].WheelIndex) - 1;
			align210 = 0;
		}

		//pack adjacent segment
		*(uint64*)(bitarray) = (last_qword << 56) | ~((uint64)1 << 63);
		last_qword = *(bitarray + bytes + 7);
		primes += segmentProcessed(bitarray, start - 8 * WHEEL30, bytes + 8, pcall);

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
	printf("\rpi2(%s) = %llu", buff, primes);
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

static uint64 adjustSieveSize(uint sieve_size, const uint sqrtp)
{
	if (sieve_size < Threshold.L2Size && sqrtp > 10000000)
		sieve_size = setSieveSize(SIEVE_SIZE);
	else if (sieve_size >= Config.Msegs << 21)
		sieve_size = setSieveSize(SIEVE_SIZE);

	if ((sieve_size & (sieve_size - 1)) != 0 && sqrtp > sieve_size * WHEEL30 / Config.Msegs)
		sieve_size = setSieveSize(1 << ilog(sieve_size >> 10, 2));
	return sieve_size * WHEEL30;
}

static uint allocMedium(const uint medium)
{
	if (MediumPrime != NULL && (MediumPrime[0].Sp < medium || medium > 2e7)) {
		free(MediumPrime); MediumPrime = NULL;
	}
	if (MediumPrime == NULL) {
		const uint pix = PI(medium, 1.2);
		MediumPrime = (SievePrime*) malloc(sizeof(SievePrime) * (pix + PI_65536));
		MediumPrime[0].Sp = medium;
		MediumPrime[0].Si = pix + PI_65536;
	}
	return medium;
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
	uchar* bitarray = (uchar*) malloc(max_cache + 1024);

	const bool bmsieve = sqrtp >= l1_maxp;
	if (bmsieve) {
		setWheelSmall(l1_maxp, isqrt(medium) + 10);
		setWheelMedium(bitarray, sieve_size, medium, start - start % WHEEL210); //TODO
	}

	//init bucket sieve
	uint64& buckets = Threshold.BucketStart;
	if (sqrtp > medium) {
		if (Config.Bsegs * 2 < Config.Msegs) Config.Bsegs = 1 << ilog(Config.Msegs, 2);
		setBucketInfo(sieve_size / Config.Bsegs, sqrtp, end - buckets);
		setWheelSmall(medium, isqrt(sqrtp) + 1);
		memset(bitarray + (L2_DCACHE_SIZE << 10), 0, l1_maxp + sizeof(uint64));
		setWheelBig(bitarray, medium, sqrtp, buckets, end - buckets);
		if (BucketInfo.CurStock < BucketInfo.LoopSize)
			allocWheelBlock(1);
	} else
		buckets = end + 1;

	//init small sieve
	{
		setWheelSmall(start, l1_maxp);
		memset(bitarray + sieve_size / WHEEL30 + 8, 0, max_cache - sieve_size / WHEEL30);
	}

	const int64 ti = getTime();
	const uint64 primes = pi2(bitarray, start, end, pcall);

	if (BucketInfo.StockSize > MAX_STOCK / 10) {
		for (uint i = 0; i < BucketInfo.PoolSize; i ++) free(WheelPool[i]);
		memset(&BucketInfo, 0, sizeof(BucketInfo));
	}

	if (Config.Flag & PRINT_RET) {
		const int64 ta = getTime();
		printResult(start, end, primes);
		if (Config.Flag & PRINT_TIME)
			printf(" (%.2f + %.2f = %.2f sec %d kb L%d%d%d %d)",
					(ti - ts) / 1000.0, (ta - ti) / 1000.0, (ta - ts) / 1000.0, Config.SieveSize >> 10, Threshold.L1Segs, Threshold.L2Segs, Config.Msegs, Config.Bsegs);
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
#elif __GNUC__ || __TINYC__
	__asm__ (
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
		mov dword ptr [edi + 0], eax
		mov dword ptr [edi + 4], ebx
		mov dword ptr [edi + 8], ecx
		mov dword ptr [edi +12], edx
	}
#endif
}

static int getIntelCache(int* l3size, int* l1size)
{
	int regs[4];

#if 0
	cpuidInfo(regs, 0, 0); /* Maximum Input Value */
	int max_leaf = regs[0];
	cpuidInfo(regs, 1, 0); /* Additional Information */
	int family = (regs[0] >> 8) & 0xF;
	int model  = (regs[0] >> 4) & 0xF;

	cpuidInfo(regs, 2, 0); /* Cache and TLB Information */

	regs[0] &= 0xFFFFFF00; /* least significant byte of EAX is invalid */
	for (int j = 0; j < 4; j++) {
		if (regs[j] < 0) { /* invalid if most significant bit set */
			regs[j] = 0;
		}
	}

	unsigned char *descriptors = (unsigned char *) regs;
	const int mb = 1024;

	#define SETL3(s) *l3size = s; break
	#define SETL1(s) *l1size = s; break

	int use_leaf_4 = 0;
	for (int j = 1; j < 16; j++) {
		switch (descriptors[j]) {
			case 0x0A: SETL1(8);
			case 0x0C: SETL1(16);
			case 0x0D: SETL1(16);
			case 0x0E: SETL1(24);
			case 0x10: SETL1(16);
//			case 0x11: SETL1(16);
			case 0x15: SETL1(16);
			case 0x2C: SETL1(32);
			case 0x30: SETL1(32);
			case 0x60: SETL1(16);
			case 0x66: SETL1(8);
			case 0x67: SETL1(16);
			case 0x68: SETL1(32);

			case 0x22: SETL3(512 );
			case 0x23: SETL3(1 * mb);
			case 0x25: SETL3(2 * mb);
			case 0x29: SETL3(4 * mb);
			case 0x40: SETL3(0); /* no L3 cache */
			case 0x46: SETL3(4 * mb);
			case 0x47: SETL3(8 * mb);
			case 0x49:
					if (family == 0x0F && model == 0x06) {
						SETL3(4 * mb);
					}
					break;
			case 0x4A: SETL3(6 * mb);
			case 0x4B: SETL3(8 * mb);
			case 0x4C: SETL3(12 * mb);
			case 0x4D: SETL3(16 * mb);
			case 0x88: SETL3(2 * mb);
			case 0x89: SETL3(4 * mb);
			case 0x8A: SETL3(8 * mb)
			case 0x8D: SETL3(3 * mb)
			case 0xD0: SETL3(512 );
			case 0xD1: SETL3(1 * mb);
			case 0xD2: SETL3(2 * mb);
			case 0xD6: SETL3(1 * mb);
			case 0xD7: SETL3(2 * mb);
			case 0xD8: SETL3(4 * mb);
			case 0xDC: SETL3(1 * mb + 512);
			case 0xDD: SETL3(3 * mb);
			case 0xDE: SETL3(6 * mb);
			case 0xE2: SETL3(2 * mb);
			case 0xE3: SETL3(4 * mb);
			case 0xE4: SETL3(8 * mb);
			case 0xEA: SETL3(12 * mb);
			case 0xEB: SETL3(18 * mb);
			case 0xEC: SETL3(24 * mb);
			case 0xFF:
					 use_leaf_4 = 1;
					 break;
			default: break;
		}
	}

	if (*l3size > 0 && *l1size < 64)
		return 0;
#endif

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
		else if (cache_type == 1 && cache_level == 1) //data cache
			*l1size = cache_size;
		printf("\n cache_level %d, cache_type = %d cache_size = %d kb", cache_level, cache_type, cache_size);
	} while (i < 16);

	return 0;
}

static int getCpuInfo()
{
	int regs[4] = {0};
	char vendor[0x40] = {0};
	int (*pTmp)[4] = (int(*)[4])vendor;
	cpuidInfo(*pTmp ++, 0x80000002, 0);
	cpuidInfo(*pTmp ++, 0x80000003, 0);
	cpuidInfo(*pTmp ++, 0x80000004, 0);

	for (int i = 0; vendor[i]; i ++) {
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
	} else if (strstr(vendor, "Intel")) {
		getIntelCache(&l3Size, &l1Size);
		cpuidInfo(regs, 0, 0);
		int max_leaf = regs[0], trys = 0;
		if (max_leaf >= 0xB) {
			int ebx[4] = {0};
			do {
				cpuidInfo((int*)regs, 0xB, trys);
				ebx[trys++] = regs[1];
			} while ((regs[0] != 0 || regs[1] != 0) && trys < 4);
			int ht = (ebx[0] & 0x00FF), tc = (ebx[1] & 0xFFFF);
			cpuCores = tc;
		} else {
			cpuidInfo(regs, 4, 0);
			if ((regs[0] & 0x1f) != 0)
				cpuCores = ((regs[0] >> 26) & 0x3f) + 1; // EAX[31:26] + 1
			else
				cpuCores /= 2;
		}

		setSieveSize((l3Size >> 10) / cpuCores * 1);
	}

	setCacheSize(2, l2Size);
	setCacheSize(1, l1Size);

	printf(" L1Dsize/L2Size/L3Size = %d/%d/%d kb, Cpu cores = %u\n\n", l1Size, l2Size, l3Size, cpuCores);
	return 0;
}
#endif

void initPrime(int sieve_size)
{
	if (Prime[2] == 0) {
#if (X86_64 || X86) && !(__clang__ && __llvm__ == 0)
		getCpuInfo();
#endif
		eratoSimple();
		initBitTable();
		initWheel30();
		initWheel210();
		initPreSieved();
		setL1Index();
		setSieveSize(sieve_size);
	}
}

static void fixRangeTest(uint64 lowerBound, const int64 range, uint64 Ret)
{
	const uint llog10 = ilog(lowerBound, 10), rlog10 = ilog(range, 10);
	const uint64 maxrange = ipow(10, 9);
	printf("Sieving the primes within (10^%u, 10^%u+10^%u) randomly\n", llog10, llog10, rlog10);
	uint64 primes = 0, upperBound = lowerBound + range;
	if (upperBound + 1 == 0)
		printf("Sieving Pi[2^64-10^%d, 2^64-1] with range 10^%d\n", rlog10, ilog(maxrange, 10));
	else
		printf("Sieving pi[10^%d, 10^%d+10^%d] with range 10^%d randomly\n", llog10, llog10, rlog10, ilog(maxrange, 10));

	while (lowerBound < upperBound) {
		uint64 rd = rand() * rand();
		uint64 end = lowerBound + (rd * rand() * rd) % maxrange + ipow(10, 4);
		end = end - end % WHEEL210 + 6;
		if (end > upperBound)
			end = upperBound;
		setSieveSize(rand() % MAX_SEGMENT + 128);
		setCacheSegs(1, rand() % 5 + 2), setCacheSegs(2, rand() % 5 + 2);
		setCacheSegs(3, rand() % 3 + 4);
		primes += doSieve(lowerBound, end - 1, NULL);
		printf("Remaining chunk: %.2f%%\r", (int64)(upperBound - lowerBound) * 100.0 / range);
		lowerBound = end;
	}

	if (upperBound + 1 == 0)
		printf("2^64-10^%d, 2^64] = %llu\n", rlog10, primes);
	else
		printf("Pi[10^%u, 10^%u+10^%u] = %llu\n", llog10, llog10, rlog10, primes);
	assert(primes == Ret);
}

static void startBenchmark()
{
	const uint primeCounts[] =
	{
		2,         // pi2(10^1)
		8,         // pi2(10^2)
		35,        // pi2(10^3)
		205,       // pi2(10^4)
		1224,      // pi2(10^5)
		8169,      // pi2(10^6)
		58980,     // pi2(10^7)
		440312,    // pi2(10^8)
		3424506,   // pi2(10^9)
		27412679,  // pi2(10^10)
		224376048, // pi2(10^11)
		12739574,  // pi2(2^32)
		7423125,   // pi2(10^12, 10^12 + 2^32)
		6330622,   // pi2(10^13, 10^13 + 2^32)
		5458483,   // pi2(10^14, 10^14 + 2^32)
		4753245,   // pi2(10^15, 10^15 + 2^32)
		4179060,   // pi2(10^16, 10^16 + 2^32)
		3699173,   // pi2(10^17, 10^17 + 2^32)
		3301635,   // pi2(10^18, 10^18 + 2^32)
		2961514,   // pi2(10^19, 10^19 + 2^32)
		0,
	};

	int64 ts = getTime();
	Config.Flag &= ~PRINT_RET;
	uint64 primes = 0;
	Config.Progress = 0;
	for (int i = 1; i <= 10; i ++) {
		primes = doSieve(0, ipow(10, i), NULL);
		printf("pi2(10^%02d) = %llu\n", i, primes);
	}

	srand(time(0));
	for (int j = 12; primeCounts[j]; j ++) {
		uint64 start = ipow(10, j), end = start + ipow(2, 32);
		setCacheSegs(1, rand() % 5 + 2), setCacheSegs(2, rand() % 5 + 2), setCacheSegs(3, rand() % 4 + 3);
		primes = doSieve(start, end, NULL);
		printf("pi2(10^%d, 10^%d+2^32) = %llu                 \n", j, j, primes);
		if (primes != primeCounts[j])
			printf("   pi2(10^%d, 10^%d+2^32) = %u                 \n", j, j, primeCounts[j]);
	}

	Config.Progress = 0;
	printf("Time elapsed %.f sec\n\n", (getTime() - ts) / 1000.0);
	puts("All Big tests passed SUCCESSFULLY!");
	const uint64 pow11 = ipow(10, 11);
	const uint64 pow12 = pow11 * 10, pow9 = pow11 / 100;

	const uint64 RangeData[][3] =
	{
		{ipow(10, 11), pow11, 199708605},
		{-1 - pow12,   pow12, 670910567},
		{ipow(10, 15), pow11, 110670758},
		{ipow(10, 17), pow11, 86176910},
		{ipow(10, 19), pow11, 68985092},
		{ipow(10, 14), pow12, 1270127074},
		{ipow(10, 16), pow12, 972773783},
		{ipow(10, 18), pow12, 768599834},
	};

	for (int k = 0; k < sizeof(RangeData) / sizeof(RangeData[0]); k ++) {
		fixRangeTest(RangeData[k][0], RangeData[k][1], RangeData[k][2]);
	}

	Config.Flag |= PRINT_RET;
}


static void printInfo( )
{
	const char* sepator =
		"------------------------------------------------------------------------------------------------------------";
	puts(sepator);
	puts("Fast implementation of the segmented sieve of Eratosthenes 2^64\n"
	"Copyright (C) by 2010-2025 Huang Yuanbing 22738078@qq.com/bailuzhou@163.com\n"
	"Compile: g++ -DSIEVE_SIZE=2048 -DFSL1 -DFDIV -march=native -funroll-loops -O3 -s -pipe PrimeNumber.cpp\n");

	char buff[500];
	char* info = buff;
#ifdef __clang__
	info += sprintf(info, "clang %s", __clang_version__); //vc/gcc/llvm
#if __llvm__
	info += sprintf(info, " on llvm/");
#endif
#endif

#if _MSC_VER
	info += sprintf(info, "Compiled by vc %d", _MSC_VER);
#elif __GNUC__
	info += sprintf(info, "Compiled by gcc %d.%d.%d", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#elif __TINYC__
	info += sprintf(info, "Compiled by tcc %d", __TINYC__);
#elif  __INTEL_COMPILER
	info += sprintf(info, "Compiled by intel %d", __INTEL_COMPILER);
#endif

#if __cplusplus
	info += sprintf(info, " c++ version %d", (int)__cplusplus);
#endif

#if X86_64
	info += sprintf(info, " x86-64");
#elif X86
	info += sprintf(info, " x86");
#elif __arm64__
	info += sprintf(info, " arm64");
#elif __arm__
	info += sprintf(info, " arm");
#else
	info += sprintf(info, " unknow");
#endif

	info += sprintf(info, " %s %s in %s\n", __TIME__, __DATE__, __FILE__);
	info += sprintf(info, "[MARCO] : MEM_WHEEL = %dM, WHEEL_SIZE = %dk, SIEVE_SIZE = %dk, WHEEL = %d\n",
			MEM_BLOCK * WHEEL_SIZE >> 17 , WHEEL_SIZE >> 7, SIEVE_SIZE, WHEEL);
	info += sprintf(info, "[CACHE] : L1Size = %u, L2Size = %u, SieveSize = %u, Bucket = %u, Block = %u\n",
			Threshold.L1Size >> 10, Threshold.L2Size >> 10, Config.SieveSize >> 10, BucketInfo.LoopSize, BucketInfo.StockSize);
	info += sprintf(info, "[ARGS ] : L1Segs/L2Segs/Mseg = (%u,%u,%u)\n",
		Threshold.L1Segs, Threshold.L2Segs, Config.Msegs);
	info += sprintf(info, "[ARGS ] : L1Maxp/L2Maxp/Medium/Large/SieveSize = (%u,%u,%u,%u,%u)",
		 Threshold.L1Maxp, Threshold.L2Maxp, Threshold.Medium, isqrt(Threshold.BucketStart + 1), Config.SieveSize);
	*info = 0;
	puts(buff);
	puts(sepator);
}

static void doCompile(const char* flag)
{
	char programming[257];
	sprintf(programming,"%s", __FILE__);
	char* pdot = strchr(programming, '.');
	if (pdot) {
		sprintf(pdot, "_.exe");
	}

	const char* const cxxflag =
#if _MSC_VER
		"cl /O2 /Oi /Ot /Oy /GT /GL %s %s %s";
#elif X86_64
		"g++ -m64 -march=native %s -O3 -funroll-loops -s -pipe %s -o %s";
#else
		"g++ -m32 -march=native %s -O3 -funroll-loops -s -pipe %s -o %s";
#endif

	char ccmd[257] = {0};
	if (flag == NULL)
		flag = "";
	sprintf(ccmd, cxxflag, flag, __FILE__, programming);

	puts(ccmd);
	system(ccmd);
	system(programming);
}

//get the first digit number index
static int parseCmd(const char params[][60])
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
			case 'H': puts(Help); if (n) puts(Benchmark);   break;
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

static bool executeCmd(const char* cmd)
{
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

		if (end == 0)
			end = start, start = 0;
		else if (end < start)
			end += start;

		if (cmdc == 'B') {
			puts(Benchmark);
			if (isdigit(params[cmdi + 2][0])) {
				int powi = atoi(params[cmdi + 2]);
				uint64 range = powi > 12 ? ipow(2, powi) : ipow(10, powi);
				for (int i = 32; i < 64 && powi > 0; i ++) {
					uint64 start2 = ipow(2, i);
					doSieve(start2, start2 + range, NULL);
				}
			} 
			if (isdigit(params[cmdi + 1][0])) {
				int powi = atoi(params[cmdi + 1]);
				uint64 range = powi > 12 ? ipow(2, powi) : ipow(10, powi);
				for (int j = 11; j < 20 && powi > 0; j ++) {
					uint64 start2 = ipow(10, j);
					doSieve(start2, start2 + range, NULL);
				}
			} else
				startBenchmark();
		} else if (cmdc == 'P') {
			PrimeCall pcall = {PCALL_BACK, (uchar*)dumpPrime, 0};
			doSieve(start, end, &pcall);
		} else if (cmdc == 'Z') {
			doCompile(params[cmdi + 1]);
		} else if (cmdi >= 0 && end > 0) {
			doSieve(start, end, NULL);
		}

		if (pcmd) {
			cmd = pcmd + 1;
		} else {
			break;
		}
	}

	return true;
}

int main(int argc, char* argv[])
{
	initPrime(SIEVE_SIZE);

	if (argc > 1)
		executeCmd(argv[1]);

	if (argc > 2) {
	srand(time(0));
	Config.Flag &= ~PRINT_RET;
	for (int j = 1; j < 100; j++) {
		for (int i = 2; i <= 6; i++) {
			setCacheSegs(2, i + 1); setCacheSegs(1, rand() % 5 + 2), setCacheSegs(3, rand() % 5 + 2);
			const uint64 medium = Config.SieveSize * WHEEL30 / Config.Msegs;
			uint64 start = medium * medium  / i;
			start -= start % WHEEL210;
			const uint64 range = Config.SieveSize * WHEEL30 * (rand() % 8 + 1);
			assert(start > range);
			const uint64 r1 = doSieve(start - range, start + range, NULL);

			setSieveSize(Config.SieveSize >> 9);
			setCacheSegs(3, i); setCacheSegs(1, rand() % 5 + 2), setCacheSegs(2, rand() % 5 + 2);
			const uint64 r2 = doSieve(start - range, start + range, NULL);
			if (r1 != r2)
				printf(" r1 = %lld != %lld %lld %lld s%d\n", r1, r2, start - range, 2*range, Config.SieveSize >> 10);
		}
	}}

#if GCOV
	executeCmd("m5 s1024 l14; 1e12 1e12+1e9"); putchar('\n');
	executeCmd("1e18 1e9"); putchar('\n');
	executeCmd("c32 l22 1e16 10^9"); putchar('\n');
	executeCmd("l34 G 1e12 e8"); putchar('\n');
	executeCmd("da s256 1e5; 1e8"); putchar('\n');
	executeCmd("s4000 c2000 p 0 100"); putchar('\n');

	executeCmd("C32 s1024 H I 1e18 1e8"); putchar('\n');
	executeCmd("df p 1e12+100 2e2*2; c1222222; g 1e16 1e8; s256 l32 l24 l14 1e8"); putchar('\n');
	executeCmd("0-1e4 0-1; 2^30");
#endif

	if (Threshold.L2Size == 512 << 10 || Config.SieveSize == 4096 << 10)
		executeCmd("L41 L42 i 2^31;");

	executeCmd("1e12 1e10; e14 e10; e10+0");
	executeCmd("10^12 1e9; e16 e9; e18 e9*1; 0-e9 0-1;");

	while (true) {
		char ccmd[257];
		printf("\n>> ");
		if (!fgets(ccmd, sizeof(ccmd), stdin) || !executeCmd(ccmd))
			break;
	}

	return 0;
}
