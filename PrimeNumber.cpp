/***************************************************************
Prime Number or PI(x) : calculate number of prime number up to x

copyright (C) 2010-2012 by Huang Yuanbing
mail to: bailuzhou@163.com
free use for non-commercial purposes

Optimization:
1. segment size is adjustable to CPU second level cache size
   (which is about 1/4-2/3 L2 cache size for saving memory and best performance)
2. set segment size equal to multiples of the first 8th prime multication and
presieved numbers is multiple of prime list: 2,3,5,7,11,13,17,19
3. the crossing out of number (which is multiples of prime) mod 30 equal to
	:1,7,11,13,17,19,23,29
4. presieved the first multiple of current prime and prime gap
5. use bit operation/packing for the crossing-out flags which can
reduce loops and improve cache hit rate
6. use thread/openMP to parallel
7. sieve more numbers in same loop for improving loop cache hit rate
8. cache segment result in each segment for fast query
9. accelerates the hotpot by assemble code
10.caching some important table: number bit 1 table, sieve gap pattern table
11.for large prime, precache the next sieveindex
12.use bucket algorithm sieve big prime factor

Todo in plan:
	1. Add UI
	3. Add detail comment and ideal of this algorithm
	4. multiple thread for pi2
	5. wheel 210 for eratSieveMedium
	6. simplify api

Link:
http://primes.utm.edu/KthPrime/algorithm.php
http://en.wikipedia.org/wiki/Sieve_of_Atkin
http://code.google.com/p/primesieve/
http://aggregate.org/MAGIC/
http://www.ieeta.pt/~tos/software/prime_sieve.html
***********************************************************************/

#include <ctype.h>
#include <math.h>
#include <memory.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

# define WHEEL           30
# define WHEEL210        210
# define PRIME7          7
# define PRIME_PRODUCT   (WHEEL * 7 * 11 * 13 * 17 * 19)
# define MAX_THREADS     16
# define HALF_SEG        2
# define TEST_FILE       "prime.pi"
# define VERSION         "1.8"

//use of the SSE4.2/ SSE4a POPCNT instruction for fast bit counting.
#if _MSC_VER > 1400
	# define POPCNT      1
	# include <intrin.h>
#elif (__GNUC__ * 10 + __GNUC_MINOR__ > 44)
	# define POPCNT      1
	# include <popcntintrin.h>
#else
	# define POPCNT      0
#endif

#if __x86_64__ || _M_AMD64
	# define TREE2       1
	# define X86_64
#else
	# define TREE2       0
#endif

#ifdef _MSC_VER
	# pragma warning(disable: 4616 4244 4018 6328 6031)
#endif

#if defined _M_AMD64
	# define ASM_X86     0
#elif _MSC_VER >= 1200
	# define ASM_X86     1
#else
	# define ASM_X86     0
#endif

#if CPU == 0 //intel core i5/7
	# define LIANGBCH        1
	# define SIEVE_SIZE      1024
	# define MAX_CACHE       (1100 << 10)
#elif CPU == 1 //amd k8/10
	# define LIANGBCH        0
	# define SIEVE_SIZE      256
	# define MAX_CACHE       (540 << 10)
#else //intel old pentium4
	# define LIANGBCH        1
	# define SIEVE_SIZE      512
	# define MAX_CACHE       (540 << 10)
#endif

# define MASK_N(n)       (1 << (n & 7))
# define PACK_DWORD(bitarray, size) \
		*((uint*)bitarray + (size >> 5) + 0) |= ~((1u << (size & 31)) - 1);\
		*((uint*)bitarray + (size >> 5) + 1) = ~0;

typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;

#ifdef _WIN32
	typedef __int64 int64;
	typedef unsigned __int64 uint64;
	#define CONSOLE "CON"
	#include <windows.h>
#else //linux, unix
	typedef long long int64;
	typedef unsigned long long uint64;
	#define CONSOLE "/dev/tty"
	#include <sys/time.h>
	#include <pthread.h>
	#include <unistd.h>
#endif

#if 0
	typedef uint stype;
	# define SMOVE   5
#else
	typedef uint64 stype;//15% fast
	# define SMOVE   6
#endif

#if _MSC_VER == 1200
	typedef int64 ltype;
#else
	typedef uint64 ltype;
#endif

#if 0
	typedef ushort pdtype;
	typedef uint64 pstype;
	#define MAX_PDIFF         0x0ffff
	#define PADD_DIFF(p)
	#define NEXT_PRIME(p, j)  p += Prime[++j]
#else
	typedef uchar pdtype;
	typedef uint pstype;
	#define MAX_PDIFF         255
	#define PADD_DIFF(p)      if (p % 2 == 0) p += MAX_PDIFF
	#define NEXT_PRIME(p, j)  p += Prime[++j]
#endif

#define L1_DCACHE_SIZE 32
#define L2_CACHE_SIZE  256
#define L1_SIEVE_SEG   4

static struct
{
	uint L1Size;
	uint L1Maxp;
	uint L2Size;
	uint L2Maxp;
}
CpuCache =
{
	L1_DCACHE_SIZE * (WHEEL << 10),
	(L1_DCACHE_SIZE << 10) / L1_SIEVE_SEG,
	L2_CACHE_SIZE * (WHEEL << 10),
	(L2_CACHE_SIZE << 10) / 2,
};

static const char* const Help = "\
	[H: Help]\n\
	[B: Benchmark]\n\
	[D: Debug log]\n\
	[G: Progress of calculating]\n\
	[Q: Exit]\n\
	[V: Version]\n\
	[F: Save result to prime.pi]\n\
	[A: Result compare by two algorithm]\n\
	[S: Set sievesize (16 - 2048)]\n\
	[C: Cpu L1 data cache size (16, 32, 64, 128, 256)]\n\
	[U: Unit test prime.pi (cases) (cache) (rw flag)]\n\
	[O: Factor of n]\n\
	[T: Threads number (2 - 16)]\n\
	[Y: find maxp prime gap [start, end]]\n\
	[P: Print prime in [start, end]]\n\
	[K: Kth prime number (n 1 - e8)]\n\
	[L: List (start) (end/count) (step) (type 0 1 2)]\n\n\
Example:\n\
	PrimeNumber c32 s1024 1e16 1e16+1e10\n\
	PrimeNumber y 1e19 2^32\n\
	PrimeNumber p 1e10 100";

//config
static struct
{
	//display the result
	bool ShowRet;
	//cache prime table for performance
	bool ResultCached;
	//display calculating time
	bool ShowTime;
	//save result
	bool SaveResult;
	//show log
	bool Showlog;
	//
	bool CheckResult;
	//display calculating time
	uint Progress;

	//number of threads
	int Threads;
	//sieve size
	int SieveSize;
	//SSE4 popcnt instruction
	int Popcnt;

	int Maxp;
}
Config =
{
	true, false, true,
	false, true, false, 31,
	4, SIEVE_SIZE * (WHEEL << 10),
	0, 10000
};

enum PRIMEDATA_RESULT
{
	COUNTBITS = 1,
	COPYBYBIT,
	PRINTPRIME,
	SAVEPRIME,
	SAVEPRIMEGAP,
	PRIMESUM,
	FINDMAXGAP,
};

struct CmdData
{
	PRIMEDATA_RESULT Cmd;
	ltype Primes;
	ltype* Data;
};

//primes < pi(sqrt(1e19 + 1e13))
const uint SQRT_PRIMES = (151876936 + 10000);
/**
the first adjacent difference more than 256 which
Prime (i + 1) th - (i)th difference
**/
static pdtype Prime[SQRT_PRIMES];

struct WheelPrime
{
	//sieve index
	uint SieveIndex; //[0 - 5]: wheel index, [6 - 31]: sieve index / WHEEL
	//sieve prime
	uint p; //[0 - 5]: pattern index, [6 - 31]: p / sievesize
};

//the bucket data is thread ...
static struct _BucketInfo
{
	uint MaxIndex;
	uint BucketSize;
	uint Log2Size;
	uint SieveSize;
	uint CurIndex;
	bool Active;
} BucketInfo;

const uint BUCKET_SIZE = 1 << 12; //depend maxp/sievesize
const uint BLOCK_SIZE  = 1 << 12; //12: 32k, best in [10 - 13]
const uint STOCK_SIZE  = (1 << 21) / BLOCK_SIZE;
const uint FREEBLK_SIZE = (SQRT_PRIMES / BLOCK_SIZE) + STOCK_SIZE * 8;

struct Stock
{
	WheelPrime* WheelData;
	Stock* NextStock;
};

static struct _BucketStock
{
	WheelPrime* CurWheel;
	Stock* Head;
	uint WheelSize;
} Bucket[BUCKET_SIZE];

//free block list
static Stock* StockList[FREEBLK_SIZE];
static int StockSize = 0;

//each segment sieveindex: (SegIndex[i] + start) % p[j] == 0
static WheelPrime MediumWheel[(280738 / 256) * (MAX_CACHE >> 10)];

/*
presieved small prime number <=19 bit array.
the crossing out bit mod WHEEL, the first
16 bit of PreSieved map to
----------------------------------------
|01/1|07/1|11/1|13/1|17/1|19/1|23/1|29/1| = 0x1111 1111 = PreSieved[0]
----------------------------------------
|31/1|37/1|41/1|43/1|47/1|49/0|53/1|59/1| = 0x1101 1111 = PreSieved[1]
----------------------------------------
****/
static uchar PreSieved[PRIME_PRODUCT / WHEEL];
//position of least significant 1 bit of an integer
static uchar Lsb[1 << 16];
//number of bits 1 binary representation table in Range[0-2^16)
static uchar WordNumBit1[1 << 16];
//cache small segment primes
static int64 PiCache[10000];

//accelerates to save diff prime < 2^31 into Prime[]
//range = [n, n + 2 * WHEEL]
static struct _WheelGap
{
	pdtype Gap[13];
	//number of bits1 in wheel range
	uchar Bits;
	//fist bit index 1 in range
	uchar Beg;
	//last bit index 1 in range
	uchar End;
} WheelGap[1 << 16];

struct WheelFactorization
{
	uchar MultipleIndex;
	uchar WheelIndex;
	uchar Correct;
	uchar UnsetBit;
};

//static WheelFactorization NextWheelData30[8][8];
static WheelFactorization NextWheel210[48][48];

typedef WheelFactorization FirstWheel;

static FirstWheel FirstWheel30[WHEEL][8];
static FirstWheel FirstWheel210[WHEEL210][48];
static uchar MultipleFactor210[48 * 2];

struct _WheelData
{
	char Index;
	uchar Mask;
	short Leng;
};

static _WheelData WheelData210[WHEEL210];

static const uchar Pattern[ ] =
{
	1,  7,  11, 13, 17, 19, 23, 29,
	31, 37, 41, 43, 47, 49, 53, 59
};

static const _WheelData WheelData30[ ] =
{
	{-1, -0, 0}, { 0,  1, 0}, {-1, -0, 1}, {-1, -0, 1}, {-1, -0, 1},
	{-1, -0, 1}, {-1, -0, 1}, { 1,  2, 1}, {-2, -0, 2}, {-2, -0, 2},
	{-2, -0, 2}, { 2,  4, 2}, {-2, -0, 3}, { 3,  8, 3}, {-3, -0, 4},
	{-3, -0, 4}, {-3, -0, 4}, { 4, 16, 4}, {-4, -0, 5}, { 5, 32, 5},
	{-5, -0, 6}, {-5, -0, 6}, {-5, -0, 6}, { 6, 64, 6}, {-6, -0, 7},
	{-6, -0, 7}, {-6, -0, 7}, {-6, -0, 7}, {-6, -0, 7}, { 7, 128,7}
};

//adjacent element difference of pattern,
//MultipleFactor30[i] = Pattern[j] - Pattern[j - 1]
static const uchar MultipleFactor30[ ] =
{
	6, 4, 2, 4, 2, 4, 6, 2,
	6, 4, 2, 4, 2, 4, 6, 2,
	6, 4, 2, 4, 2, 4, 6, 2,
};

static struct ThreadInfo
{
	int Starti;
	int Threads;
	ltype Endi;
	ltype Primes;
} ThreadData[MAX_THREADS];

static int sievePrime(pdtype[], uint);
static ltype initPiCache(ltype, ltype, int, CmdData*);

#ifdef _WIN32
static DWORD WINAPI threadProc(void* ptinfo)
#else
static void* threadProc(void* ptinfo)
#endif
{
	struct ThreadInfo* pThreadInfo = (struct ThreadInfo*)ptinfo;
	pThreadInfo->Primes = initPiCache(pThreadInfo->Starti,
			pThreadInfo->Endi, pThreadInfo->Threads, NULL);
	if (Config.Showlog) {
		printf("thread %d is finished\n", pThreadInfo->Starti);
	}

	return 0;
}

static void divideTask(ltype Endi, int threads)
{
	for (int i = 0; i < threads; i++) {
		ThreadData[i].Starti = i + 1;
		ThreadData[i].Endi = Endi;
		ThreadData[i].Threads = threads;
		ThreadData[i].Primes = 0;
	}

	if (PiCache[2] < PiCache[3]) {
		Config.ResultCached = false;
	}
}

static ltype startWorkThread(ltype maxn, int threads)
{
	int i;
	ltype primes = 0;

	if (threads > MAX_THREADS) {
		threads = 4;
	}

	divideTask(maxn / Config.SieveSize + 1, threads);

#ifdef _WIN32
	HANDLE thandle[MAX_THREADS];
	DWORD tid[MAX_THREADS];
	for (i = 0; i < threads; i++) {
		thandle[i] = CreateThread(NULL, 0, threadProc,
				(LPVOID)(&ThreadData[i]), 0, &tid[i]);
		if (thandle[i] == NULL) {
			printf("create win32 thread error %ld\n", GetLastError( ));
		}
	}
	for (i = 0; i < threads; i++) {
		WaitForSingleObject(thandle[i], INFINITE);
		CloseHandle(thandle[i]);
		primes += ThreadData[i].Primes;
	}
#else
	pthread_t tid[MAX_THREADS];
	for (i = 0; i < threads; i++) {
		int error = pthread_create(&tid[i], NULL, threadProc, &ThreadData[i]);
		if (error != 0) {
			printf("create posix thread error %d\n", error);
		}
	}
	for (i = 0; i < threads; i++) {
		pthread_join(tid[i], NULL);
		primes += ThreadData[i].Primes;
	}
#endif

	for (i = 2; PiCache[i] > 0; i++) {
		if (PiCache[i] < PiCache[i - 1])
			PiCache[i] += PiCache[i - 1];
	}

	Config.ResultCached = true;

	return primes;
}

//get current time
static double getTime( )
{
#ifdef _WIN32
	LARGE_INTEGER s_freq;
	LARGE_INTEGER performanceCount = {0};
	QueryPerformanceFrequency(&s_freq);
	QueryPerformanceCounter(&performanceCount);
	return 1000. * performanceCount.QuadPart / s_freq.QuadPart;
#else
	struct timeval tmVal;
	gettimeofday(&tmVal, NULL);
	return tmVal.tv_sec * 1000. + tmVal.tv_usec / 1000.;
#endif
}

static uint ilog2(uint n)
{
	uint log2 = 0;
	if (n >= (1 << 16)) { n >>= 16; log2 |= 16; }
	if (n >= (1 <<  8)) { n >>=  8; log2 |=  8; }
	if (n >= (1 <<  4)) { n >>=  4; log2 |=  4; }
	if (n >= (1 <<  2)) { n >>=  2; log2 |=  2; }
	if (n >= (1 <<  1)) { n >>=  0; log2 |=  1; }
	return log2;
}

static ltype ipow(ltype x, uint n)
{
	ltype pown = 1;
	while (n != 0) {
		if (n & 1) {
			pown *= x;
			n -= 1;
		}
		x *= x;
		n /= 2;
	}

	return pown;
}

static uint isqrt(ltype n)
{
	uint rem = 0, root = 0, divisor = 0;

	for (int i = 0; i < 32; i++) {
		root <<= 1;
		rem = ((rem << 2) + (uint)(n >> 62));
		n <<= 2;
		divisor = (root << 1) + 1;
		if (divisor <= rem) {
			rem -= divisor;
			root++;
		}
	}

	return root;
}

//make sure no any overflow. it's 2 time fast than _lldiv
static inline uint asmllDiv(const uint highw, const uint loww, uint p)
{
#ifdef NX86
	p = (((ltype)highw << 32) + loww) % p;
#elif !defined _MSC_VER
	__asm
	(
		"divl %%ecx\n"
		: "=d" (p)
		: "d"(highw), "a"(loww), "c"(p)
	);
#else
	__asm
	{
		mov eax, loww
		mov edx, highw
		div p
		mov p, edx
	}
#endif

	return p;
}

//return min odd sieveindex : (start + sieveindex) % p == 0
//use one x86/old cpu/compiler
static inline uint fastDiv(const ltype start, const uint p)
{
#if 1 || _MSC_VER
	const uint edx = (uint)(start >> 32);
	uint sieveindex = p;
	if (edx < p) {
		sieveindex -= asmllDiv(edx, (uint)start, p);
	} else {
		sieveindex -= (uint)(start % p);
	}
#else
	uint sieveindex = p - (uint)(start % p);
#endif

	return sieveindex;
}

//the only invalid is [space][number][e][+-*][number]
//e9, 2^32, 1234, 10000*2, 2^30-1E2, 2e9+2^20 all are invalid
ltype atoint64(const char* str)
{
	ltype ret = 0;

	while (isspace(*str)) {
		str++;
	}

	while (isdigit(*str)) {
		ret = ret * 10 + *str++ - '0';
	}

	if (*str && isdigit(str[1])) {
		if (str[0] == '^') {//a*b^c
			ret = ipow(ret, atoi(str + 1));
		} else if (str[0] == 'e' || str[0] == 'E') {
			//a*'[Ee]'b
			if (ret == 0) {
				ret = 1;
			}
			ret *= ipow(10, atoi(str + 1));
		}
	}

	const char* ps = str;
	if (ps = strchr(str, '+')) {
		ret += atoint64(ps + 1);
	} else if (ps = strchr(str, '-')) {
		ret -= atoint64(ps + 1);
	} else if (ps = strchr(str, '*')) {
		ret *= atoi(ps + 1);
	}

	return ret;
}

//not safe
static inline void
crossOutWheelFactor(uchar* ppbeg[], const uchar* pend, const uint step)
{
	uchar* ps0 = ppbeg[0], *ps1 = ppbeg[1];
	uchar* ps2 = ppbeg[2], *ps3 = ppbeg[3];

	while (ps3 <= pend) {
		*ps0 |= 1 << 0, ps0 += step;
		*ps1 |= 1 << 1, ps1 += step;
		*ps2 |= 1 << 2, ps2 += step;
		*ps3 |= 1 << 3, ps3 += step;
	}
	*ps0 |= 1 << 0; *ps1 |= 1 << 1; *ps2 |= 1 << 2;

	ps0 = ppbeg[4], ps1 = ppbeg[5];
	ps2 = ppbeg[6], ps3 = ppbeg[7];
	while (ps3 <= pend) {
		*ps0 |= 1 << 4, ps0 += step;
		*ps1 |= 1 << 5, ps1 += step;
		*ps2 |= 1 << 6, ps2 += step;
		*ps3 |= 1 << 7, ps3 += step;
	}
	*ps0 |= 1 << 4; *ps1 |= 1 << 5; *ps2 |= 1 << 6;
}

//30% fast on old cpu pentium 4
static inline void
crossOutWheelFactor2(uchar *pbeg, const uchar *pend, const uint p)
{
	const uint o = p / WHEEL;
	uchar* ps1 = pbeg, *ps2 = pbeg;

	switch (WheelData30[p % WHEEL].Index) {
		case 0 :
			ps2 += o * 6 + 0;
			while (ps1 < pend) {
				*ps1 |= 1 << 0, ps1 += o * 10 + 0;
				*ps2 |= 1 << 1, ps2 += o *  6 + 0;
				*ps1 |= 1 << 2, ps1 += o *  6 + 0;
				*ps2 |= 1 << 3, ps2 += o *  6 + 0;
				*ps1 |= 1 << 4, ps1 += o *  6 + 0;
				*ps2 |= 1 << 5, ps2 += o * 10 + 0;
				*ps1 |= 1 << 6, ps1 += o *  8 + 1;
				*ps2 |= 1 << 7, ps2 += o *  8 + 1;
			}
		break;
		case 1 :
			ps2 += o * 4 + 0;
			while (ps1 < pend) {
				*ps1 |= 1 << 0, ps1 += o *  6 + 1;
				*ps2 |= 1 << 7, ps2 += o *  6 + 2;
				*ps1 |= 1 << 3, ps1 += o * 10 + 2;
				*ps2 |= 1 << 2, ps2 += o *  8 + 2;
				*ps1 |= 1 << 6, ps1 += o *  8 + 2;
				*ps2 |= 1 << 1, ps2 += o * 10 + 2;
				*ps1 |= 1 << 5, ps1 += o *  6 + 2;
				*ps2 |= 1 << 4, ps2 += o *  6 + 1;
			}
		break;
		case 2 :
			ps2 += o * 2 + 0;
			while (ps1 < pend) {
				*ps1 |= 1 << 0, ps1 += o *  6 + 2;
				*ps2 |= 1 << 6, ps2 += o *  6 + 2;
				*ps1 |= 1 << 1, ps1 += o *  6 + 2;
				*ps2 |= 1 << 7, ps2 += o * 10 + 4;
				*ps1 |= 1 << 3, ps1 += o *  8 + 3;
				*ps2 |= 1 << 5, ps2 += o *  8 + 3;
				*ps1 |= 1 << 2, ps1 += o * 10 + 4;
				*ps2 |= 1 << 4, ps2 += o *  6 + 2;
			}
		break;
		case 3 :
			ps2 += o * 4 + 1;
			while (ps1 < pend) {
				*ps1 |= 1 << 0, ps1 += o *  6 + 2;
				*ps2 |= 1 << 6, ps2 += o *  6 + 3;
				*ps1 |= 1 << 5, ps1 += o *  6 + 3;
				*ps2 |= 1 << 2, ps2 += o *  6 + 2;
				*ps1 |= 1 << 1, ps1 += o * 10 + 4;
				*ps2 |= 1 << 7, ps2 += o *  8 + 4;
				*ps1 |= 1 << 4, ps1 += o *  8 + 4;
				*ps2 |= 1 << 3, ps2 += o * 10 + 4;
			}
		break;
		case 4 :
			ps2 += o * 6 + 3;
			while (ps1 < pend) {
				*ps1 |= 1 << 0, ps1 += o *  8 + 4;
				*ps2 |= 1 << 3, ps2 += o *  8 + 4;
				*ps1 |= 1 << 4, ps1 += o * 10 + 6;
				*ps2 |= 1 << 7, ps2 += o *  6 + 4;
				*ps1 |= 1 << 1, ps1 += o *  6 + 3;
				*ps2 |= 1 << 2, ps2 += o *  6 + 3;
				*ps1 |= 1 << 5, ps1 += o *  6 + 4;
				*ps2 |= 1 << 6, ps2 += o * 10 + 6;
			}
		break;
		case 5 :
			ps2 += o * 4 + 2;
			while (ps1 < pend) {
				*ps1 |= 1 << 0, ps1 += o * 10 + 6;
				*ps2 |= 1 << 4, ps2 += o *  8 + 5;
				*ps1 |= 1 << 2, ps1 += o *  8 + 5;
				*ps2 |= 1 << 5, ps2 += o * 10 + 6;
				*ps1 |= 1 << 3, ps1 += o *  6 + 4;
				*ps2 |= 1 << 7, ps2 += o *  6 + 4;
				*ps1 |= 1 << 1, ps1 += o *  6 + 4;
				*ps2 |= 1 << 6, ps2 += o *  6 + 4;
			}
		break;
		case 6 :
			ps2 += o * 2 + 1;
			while (ps1 < pend) {
				*ps1 |= 1 << 0, ps1 += o *  6 + 4;
				*ps2 |= 1 << 4, ps2 += o * 10 + 8;
				*ps1 |= 1 << 5, ps1 += o *  8 + 6;
				*ps2 |= 1 << 1, ps2 += o *  8 + 6;
				*ps1 |= 1 << 6, ps1 += o * 10 + 8;
				*ps2 |= 1 << 2, ps2 += o *  6 + 4;
				*ps1 |= 1 << 3, ps1 += o *  6 + 5;
				*ps2 |= 1 << 7, ps2 += o *  6 + 5;
			}
		break;
		case 7 :
			ps2 += o * 2 + 1;
			while (ps1 < pend) {
				*ps1 |= 1 << 0, ps1 += o *  8 + 7;
				*ps2 |= 1 << 7, ps2 += o * 10 + 10;
				*ps1 |= 1 << 6, ps1 += o *  6 + 6;
				*ps2 |= 1 << 5, ps2 += o *  6 + 6;
				*ps1 |= 1 << 4, ps1 += o *  6 + 6;
				*ps2 |= 1 << 3, ps2 += o *  6 + 6;
				*ps1 |= 1 << 2, ps1 += o * 10 + 10;
				*ps2 |= 1 << 1, ps2 += o *  8 + 7;
			}
		break;
	}
}

#if ASM_X86 == 1
	#define ESP_OFFSET4 0x20
_declspec(naked)
#endif
static void
crossOut4Factor(uchar* ppbeg[], const uint mask,
		const uchar* pend, const uint step)
{
#if ASM_X86
	__asm
	{
		pushad //store all register into stack, the esp will decrease 0x20
		/*****
		ppbeg   esp + 32 + 04
		mask   esp + 32 + 08
		pend   esp + 32 + 12
		step   esp + 32 + 16
		*****/
		// define ebx, ebp, esi, edi as ppbeg[0], ppbeg[1], ppbeg[2], ppbeg[3]
		mov eax, dword ptr [esp + 4 + ESP_OFFSET4]
		mov ebx, dword ptr [eax + 0]
		mov ebp, dword ptr [eax + 4]
		mov esi, dword ptr [eax + 8]
		mov edi, dword ptr [eax + 12]

		//define edx, eax, ecx as mask, pend, step
		mov edx, dword ptr [esp + 8 + ESP_OFFSET4]
		mov eax, dword ptr [esp + 12 + ESP_OFFSET4]
		mov ecx, dword ptr [esp + 16 + ESP_OFFSET4]

		jmp a20
loop_10:
		or byte ptr [ebx], dl
		or byte ptr [ebp], dh
		bswap edx
//		ror edx, 16

		add ebx, ecx
		add ebp, ecx

		or byte ptr [esi], dh
		or byte ptr [edi], dl
		bswap edx
//		ror edx, 16

		add esi, ecx
		add edi, ecx
		//sse2 or mmx

a20:	cmp edi, eax
		jbe loop_10

		cmp ebx, eax
		jg a30
		or byte ptr [ebx], dl

		cmp ebp, eax
		jg a30
		or byte ptr [ebp], dh

		cmp esi, eax
		jg a30
		ror edx, 16
		or byte ptr [esi], dl
a30:
		popad
		ret
	}
#else
	union
	{
		uint dmask;
		uchar bmask[4];
	} umask;
	umask.dmask = mask;

	uchar* ps0 = ppbeg[0], *ps1 = ppbeg[1];
	uchar* ps2 = ppbeg[2], *ps3 = ppbeg[3];

	while (ps3 <= pend) {
#if 1
		*ps0 |= umask.bmask[0]; ps0 += step;
		*ps1 |= umask.bmask[1]; ps1 += step;
		*ps2 |= umask.bmask[2]; ps2 += step;
		*ps3 |= umask.bmask[3]; ps3 += step;
#else
		//uint mask32 = mask;
		*ps0 |= mask >> 0, ps0 += step;
		*ps1 |= mask >> 8, ps1 += step;
		*ps2 |= mask >>16, ps2 += step;
		*ps3 |= mask >>24, ps3 += step;
#endif
	}

	if (ps0 <= pend)
		*ps0 |= mask;
	if (ps1 <= pend)
		*ps1 |= mask >> 8;
	if (ps2 <= pend)
		*ps2 |= mask >> 16;

#endif
}

#if ASM_X86 == 1
	#define ESP_OFFSET2 8
_declspec(naked)
#endif
static inline void
crossOut2Factor(uchar* ps0, uchar* ps1, const uchar* pend,
		const ushort wordmask, const int step)
{
#if ASM_X86
	__asm
	{
		push esi
		push edi

		mov edi, dword ptr [esp + ESP_OFFSET2 + 4]	//ps0
		mov esi, dword ptr [esp + ESP_OFFSET2 + 8]	//ps1
		mov edx, dword ptr [esp + ESP_OFFSET2 + 12]	//pend
		mov eax, dword ptr [esp + ESP_OFFSET2 + 16]	//wordmask
		mov ecx, dword ptr [esp + ESP_OFFSET2 + 20]	//step

		jmp LCMP

LOOP2:
		or byte ptr [edi], al
		or byte ptr [esi], ah
		add edi, ecx
		add esi, ecx
LCMP:
		cmp esi, edx
		jle LOOP2

		cmp edi, edx
		jg RETP
		or byte ptr [edi], al
RETP:
		pop edi
		pop esi
		ret
	}
#else

	const uchar masks0 = (uchar)wordmask;
	const uchar masks1 = wordmask >> 8;

	while (ps1 <= pend) {
		*ps1 |= masks1; ps1 += step;
		*ps0 |= masks0; ps0 += step;
	}
	if (ps0 <= pend)
		*ps0 |= masks0;
#endif
}

static void copyFromTpl(uchar bitarray[], int bytes)
{
	while (bytes > sizeof(PreSieved)) {
		memcpy(bitarray, PreSieved, sizeof(PreSieved));
		bytes -= sizeof(PreSieved);
		bitarray += sizeof(PreSieved);
	}

	memcpy(bitarray, PreSieved, bytes);
}

/**
 * Pre-sieve multiples of small primes <= limit_ (e.g. 19).
 * Resets the sieve array (resets bits to 1) of SieveOfEratosthenes
 * objects after each sieved segment and removes the multiples of
 * small primes without sieving.
 */
static int preSieve(uchar bitarray[], ltype start, const int sievesize)
{
	const int segmentoffset = (int)(start % PRIME_PRODUCT) / WHEEL;
	int bits = sievesize + (int)(start % WHEEL); //pack sievesize
	bits = (bits / WHEEL * 8 + WheelData30[bits % WHEEL].Leng);
	const int bytes = (bits + 7) / 8;

	if (segmentoffset + bytes < sizeof(PreSieved)) {
		memcpy(bitarray, PreSieved + segmentoffset, bytes);
	} else {
		memcpy(bitarray, PreSieved + segmentoffset, sizeof(PreSieved) - segmentoffset);
		copyFromTpl(bitarray + sizeof(PreSieved) - segmentoffset,
				bytes + segmentoffset - sizeof(PreSieved));
	}

	//set bit position 1, 2
	if (start < WHEEL) {
		bitarray[0] = 0x3;
	}

	//pack the first byte with bit 1
	bitarray[0] |= (1 << WheelData30[start % WHEEL].Leng) - 1;
	//pack the last 2 dword with bit 1
	PACK_DWORD(bitarray, bits);
	PACK_DWORD(bitarray, (bits + 32));

	return bytes;
}

//3.Use of the SSE 4.2 POPCNT instruction for fast bit counting
#if POPCNT
static inline int countBitsPopcnt(const stype n)
{
	//popcnt instruction : INTEL ix/SSE4.2, AMD Phonem/SSE4A
#if x86_64
	return _mm_popcnt_u64(n);
#elif (SMOVE == 5)
	return _mm_popcnt_u32(n);
#else
	return _mm_popcnt_u32(n) + _mm_popcnt_u32((uint)(n >> 32));
#endif
}
#endif

static inline int countBitsTable(const stype n)
{
#if SMOVE == 5
	return WordNumBit1[(ushort)n] + WordNumBit1[n >> 16];
#else
	uint hig = (uint)(n >> 32), low = (uint)n;
	return WordNumBit1[(ushort)hig] + WordNumBit1[(ushort)low] +
		WordNumBit1[hig >> 16] + WordNumBit1[low >> 16];
#endif
}

#if TREE2
static inline int countBitsTree2(stype n)
{
#if SMOVE == 6
	n -= (n >> 1) & 0x5555555555555555;
	n = (n & 0x3333333333333333) + ((n >> 2) & 0x3333333333333333);
	n = (n + (n >> 4)) & 0x0F0F0F0F0F0F0F0F;
	n += n >> 8;
	n += n >> 16;
	n += n >> 32;
	return ((uint)n) & 0x00000000FF;
#elif 1
	n -= (n >> 1) & 0x55555555;
	n = (n & 0x33333333) + ((n >> 2) & 0x33333333);
	n = (n + (n >> 4)) & 0x0F0F0F0F;
	n += n >> 8;
	n += n >> 16;
	return n & 0x0000003F;
#else
	const uint mask7 = 0x77777777;
	uint tmp = (n >> 1) & mask7;
	n -= tmp;
	tmp = (tmp >> 1) & mask7;
	n -= tmp;
	tmp = (tmp >> 1) & mask7;
	n -= tmp;
	n = (n + (n >> 4)) & 0x0F0F0F0F;
	n *= 0x01010101;
	return n >> 24;
#endif
}
#endif

//count number of bit 0 in binary representation of array
//bit 0 mean it's a prime position
static int countZeroBitsArray(const stype bitarray[], const int bitsize)
{
	int bit1s = 0;
	int loops = bitsize >> SMOVE;

	while (loops-- >= 0) {
#if POPCNT
		bit1s += countBitsPopcnt(*bitarray++);
#elif TREE2
		bit1s += countBitsTree2(*bitarray++);
#else
		bit1s += countBitsTable(*bitarray++);
#endif
	}

	return ((1 + (bitsize >> SMOVE)) << SMOVE) - bit1s;
}

//get prime from bit buffer
static int savePrime(const ushort bitarray[], ltype sieveindex, const int wordsize, ltype* prime)
{
	int primes = 0;
	sieveindex -= sieveindex % WHEEL;

	for (int bi = 0; bi <= wordsize; bi++) {
		ushort mask = ~bitarray[bi];
		while (mask > 0) {
			prime[primes++] = sieveindex + Lsb[mask];
			mask &= mask - 1;
		}
		sieveindex += WHEEL * 2;
	}

	return primes;
}

//print prime from bit buffer
static int printPrime(const ushort bitarray[], ltype sieveindex, const int wordsize, ltype sumprime)
{
	int primes = 0;
	sieveindex -= sieveindex % WHEEL;

	for (int bi = 0; bi <= wordsize; bi++) {
		ushort mask = ~bitarray[bi];
		while (mask > 0) {
			printf("%llu %llu\n", ++primes + sumprime, sieveindex + Lsb[mask]);
			mask &= mask - 1;
		}
		sieveindex += WHEEL * 2;
	}

	return primes;
}

#define IS_BIT0(a, n) ((a & (1 << n)) == 0)

//get prime from bit buffer
static ltype countPrimeSum(const ushort bitarray[], ltype sieveindex, const int wordsize, ltype* prime)
{
	ltype primes = 0;
	sieveindex -= sieveindex % WHEEL;
	//1, 7, 11, 13, 17, 19, 23, 29
	static short sexprime[1 << 16];
	if (sexprime[0] == 0) {
		static int pair[] = {0, 1, 1, 3, 2, 4, 3, 5, 4, 6, 6, 7};
		for (int i = 0; i < 255; i++) {
			int sexn = 0;
			for (int j = 0; j < 12; j += 2) {
				if (IS_BIT0(i, pair[j]) && IS_BIT0(i, pair[j + 1]))
					sexn++;
			}
			sexprime[i] = sexn;
		}
		for (int j = 256; j < (1 << 16); j++) {
			sexprime[j] = sexprime[j >> 8] + sexprime[j & 255];
		}
	}

	ltype sum = sieveindex;
	for (int bi = 0; bi <= wordsize; bi++) {
		const ushort mask = bitarray[bi];
		primes += sexprime[mask] + sum * WordNumBit1[mask];
		sum += 2 * WHEEL;
	}

	return *prime = primes;
}

//get prime from bit buffer
//wordsize > 1e4 and
static int findPrimeGap(const ushort bitarray[], ltype sieveindex, const int wordsize, ltype *primegap)
{
	ushort pdiff = (ushort)primegap[0];
	ushort maxgap = (ushort)primegap[1];
	int skipwords = maxgap / (WHEEL * 2) - 1;
	if (skipwords < 1)
		skipwords = 1;

	for (int i = 0; i < wordsize; i++) {
		const uint msk = ~bitarray[i] & 0xffff;
		if (msk != 0) {
			ushort gap = pdiff + WheelGap[msk].Beg;
			if (gap > maxgap) { //
				maxgap = gap;
				primegap[2] = sieveindex + i * WHEEL * 2 + Lsb[msk];
			}
			if (bitarray[i + skipwords] != 0xffff) {
				i += skipwords - 1;
				pdiff = 0;
			} else {
				pdiff = WheelGap[msk].End;
			}
		} else {
			pdiff += WHEEL * 2;
		}
	}
	primegap[0] = WheelGap[~bitarray[wordsize - 1] & 0xffff].End;
	primegap[1] = maxgap;
	return 0;
}

static int savePrimeGap(ushort bitarray[], const int wordsize, pdtype primegap[])
{
	pdtype* prime = primegap;
	ushort pdiff = primegap[0];

	for (int i = 0; i < wordsize; i++) {
		const uint msk = ~bitarray[i] & 0xffff;
		if (msk != 0) {
			ushort gap = pdiff + WheelGap[msk].Beg;
			if (gap > MAX_PDIFF && MAX_PDIFF < 1000)
				gap -= MAX_PDIFF;
			primegap[0] = gap;
#ifdef X86_64
			*(ltype*)(primegap + 1) = *(ltype*)(WheelGap[msk].Gap + 0);
#else
			*(pstype*)(primegap + 1) = *(pstype*)(WheelGap[msk].Gap + 0);
			*(pstype*)(primegap + 5) = *(pstype*)(WheelGap[msk].Gap + 4);
#endif
			pdiff = WheelGap[msk].End;
			const int bits = WheelGap[msk].Bits;
			if (bits > 8)
				*(pstype*)(primegap + 9) = *(pstype*)(WheelGap[msk].Gap + 8);
			primegap += bits;
		} else {
			pdiff += WHEEL * 2;
		}
	}
	primegap[0] = pdiff;

	return primegap - prime;
}

static int initBucketInfo(const uint sievesize, const uint sqrtp, ltype rangesize)
{
	BucketInfo.CurIndex = 0;
	BucketInfo.MaxIndex = (uint)(rangesize / sievesize) + 1;

	BucketInfo.BucketSize = (sqrtp * 10.0) / sievesize + 2;
	if (BucketInfo.MaxIndex >= BUCKET_SIZE)
		assert(BucketInfo.BucketSize <= BUCKET_SIZE);

	BucketInfo.BucketSize = (2 << ilog2(BucketInfo.BucketSize)) - 1;
	BucketInfo.Log2Size = ilog2(sievesize / WHEEL);
	BucketInfo.SieveSize = (1 << BucketInfo.Log2Size) - 1;

	const uint pix = (uint)(sqrtp / log((double)sqrtp) * (1 + 1.200 / log((double)sqrtp)));

	StockSize = pix / BLOCK_SIZE + BucketInfo.BucketSize;
	if (rangesize * 3 < sqrtp)
		StockSize >>= 1;
	assert(StockSize < sizeof(StockList) / sizeof(StockList[0]));

	//dynamic memory use
	WheelPrime *pwheel =
		(WheelPrime*) malloc(StockSize * BLOCK_SIZE * sizeof(WheelPrime));
	assert(pwheel);

	Stock* pstock = (Stock*) malloc(StockSize * sizeof(Stock));
	for (int j = 0; j < StockSize; j++) {
		StockList[j] = pstock + j;
	}

	const uint aligned = (ltype)pwheel & (8192 - 1);
	StockList[0]->WheelData = pwheel;
	StockList[1]->WheelData = (WheelPrime*)((ltype)pwheel + 8192 - aligned);
	for (int i = 2; i < StockSize; i++) {
		StockList[i]->WheelData = StockList[i - 1]->WheelData + BLOCK_SIZE;
	}

	//printf("new memory %d MB\n", StockSize * BLOCK_SIZE * sizeof(WheelPrime) >> 20);
	return pix;
}

static inline void
pushBucket(const uint sieveindex, const uint wheelprime, const uchar wheelindex)
{
	const uint nextbucket = (sieveindex >> BucketInfo.Log2Size) + BucketInfo.CurIndex;
	if (nextbucket < BucketInfo.MaxIndex) {
		_BucketStock* pbucket = &Bucket[nextbucket & BucketInfo.BucketSize];
		if (pbucket->WheelSize++ % BLOCK_SIZE == 0) {
			Stock* head = StockList[-- StockSize];
			pbucket->CurWheel = head->WheelData;
			head->NextStock = pbucket->Head;
			pbucket->Head = head;
		}

		WheelPrime* nextwheel = pbucket->CurWheel++;
		nextwheel->SieveIndex = (sieveindex & BucketInfo.SieveSize) << 6 | wheelindex;
		nextwheel->p = wheelprime;
	}
}

static void initMediumWheel(const uint sievesize, const uint sqrtp, const ltype start)
{
	uint j = 8 + PRIME_PRODUCT / 9699690;
	uint p = Prime[0];
	uint sievesize2 = sievesize;

	if (BucketInfo.Active) {
		sievesize2 = CpuCache.L2Size;
		const uint bitsize = (sievesize / WHEEL) * 8 - 1;
		assert ((bitsize & (bitsize + 1)) == 0);
	}

	const uint minp = sqrtp < sievesize / HALF_SEG ? sqrtp : sievesize / HALF_SEG;
	//OPT_MEDI
	for (; p < minp; NEXT_PRIME(p, j)) {
		const ltype p2 = (ltype)p * p;
		uint sieveindex;

		if (start <= p2) {
			sieveindex = (p2 - start) % sievesize2;
		} else {
			sieveindex = fastDiv(start, p);
		}

		const FirstWheel cwn = FirstWheel30[sieveindex % WHEEL][WheelData30[p % WHEEL].Index];
		sieveindex += cwn.Correct * p;

		if (p2 > start + sievesize2)
			sieveindex %= sievesize2;

		MediumWheel[j + 0].SieveIndex = (sieveindex << 3) | cwn.MultipleIndex & 7;
		MediumWheel[j + 0].p = p;
	}
	MediumWheel[j].p = -1u;
	assert(j < sizeof(MediumWheel) / sizeof(MediumWheel[0]));
}

static void initBucketStock(const uint sievesize, const uint sqrtp, const ltype start)
{
	uint j = 8 + PRIME_PRODUCT / 9699690;
	uint p = Prime[0];

	while (p < sievesize / HALF_SEG) {
		NEXT_PRIME(p, j);
	}

	uint nextp = p + 40001;
	ltype remp = start / nextp;
	for (; p < sqrtp; NEXT_PRIME(p, j)) {
		PADD_DIFF(p);
#if 0 // x86 fast sqrtp > 1e9
		if (p > nextp)
			remp = start / (nextp = p + 40001);
		uint sieveindex = fastDiv(start - remp * p, p);
#else
		uint sieveindex = fastDiv(start, p);
#endif

		const uint wheelprime = (p / WHEEL210 * 64 + WheelData210[p % WHEEL210].Index);
		const uchar sievemodle = sieveindex % WHEEL210;
		const FirstWheel cwn = FirstWheel210[sievemodle][wheelprime & 63];

		sieveindex = sieveindex / WHEEL210 * 7 + (wheelprime >> 6) * cwn.Correct * 7;
		sieveindex += (cwn.Correct * (p % WHEEL210) + sievemodle) / WHEEL;
		pushBucket(sieveindex, wheelprime, cwn.WheelIndex);
	}
	assert(StockSize > 0);
}

static inline void
sieveSmall0(uchar bitarray[], const uchar* pend, const uint p, uint sieveindex, int skipindex)
{
#if CPU != 2
	uchar* ppbeg[8];
	for (int i = 0; i < 8; i++) {
		ppbeg[WheelData30[sieveindex % WHEEL].Index] = bitarray + sieveindex / WHEEL;
		sieveindex += MultipleFactor30[skipindex++] * p;
	}
	crossOutWheelFactor(ppbeg, pend, p);
#else //p4 fast
	for (int i = 0; i < 8; i++) {
		const uchar mask = WheelData30[sieveindex % WHEEL].Mask;
		bitarray[sieveindex / WHEEL] |= mask;
		if (mask == 1)
			break;
		sieveindex += MultipleFactor30[skipindex++] * p;
	}
	crossOutWheelFactor2(bitarray + sieveindex / WHEEL, pend + 1, p);
#endif
}

static inline void
sieveSmall1(uchar bitarray[], const uchar* pend, const uint p, uint sieveindex, int skipindex)
{
	for (int i = 0; i < 4; i++) {
		uchar* ps0 = bitarray + sieveindex / WHEEL;
		ushort wordmask = WheelData30[sieveindex % WHEEL].Mask;
		sieveindex += MultipleFactor30[skipindex++] * p;

		uchar* ps1 = bitarray + sieveindex / WHEEL;
		wordmask |= WheelData30[sieveindex % WHEEL].Mask << 8;
		sieveindex += MultipleFactor30[skipindex++] * p;
		crossOut2Factor(ps0, ps1, pend, wordmask, p);
	}
}

static inline void
sieveSmall2(uchar bitarray[], const uchar* pend, const uint p, uint sieveindex, int skipindex)
{
	uchar* ppbeg[8], dwordmask[8];
	for (int i = 0; i < 8; i++) {
		ppbeg[i] = bitarray + sieveindex / WHEEL;
		dwordmask[i] = WheelData30[sieveindex % WHEEL].Mask;
		sieveindex += MultipleFactor30[skipindex++] * p;
	}

	crossOut4Factor(ppbeg + 0, *(uint*)(dwordmask + 0), pend, p);
	crossOut4Factor(ppbeg + 4, *(uint*)(dwordmask + 4), pend, p);
}

static void eratSieveSmall(uchar bitarray[], const ltype start, const int sievesize, uint sqrtp)
{
	if ((start + sievesize) < ((ltype)sqrtp) * sqrtp) {
		sqrtp = (uint)sqrt((double)start + sievesize) + 2;
	}

	for (uint p = Prime[0], j = 8 + PRIME_PRODUCT / 9699690; p < sqrtp; NEXT_PRIME(p, j)) {
		uint sieveindex = fastDiv(start, p);
		if (start <= p) {
			sieveindex = p * p - (uint)start;
		}

		const FirstWheel cwn = FirstWheel30[sieveindex % WHEEL][WheelData30[p % WHEEL].Index];
		sieveindex += cwn.Correct * p;

#if LIANGBCH
		sieveSmall0(bitarray, bitarray + sievesize / WHEEL, p, sieveindex, cwn.MultipleIndex);
#else
		sieveSmall1(bitarray, bitarray + sievesize / WHEEL, p, sieveindex, cwn.MultipleIndex);
#endif
	}
}

static void eratSieveMedium2(uchar bitarray[], const ltype start, const uint sievesize, const uint minp, uint sqrtp)
{
	const uint bytes = sievesize / WHEEL + 0;
	uint j = 8 + PRIME_PRODUCT / 9699690, p = Prime[0];//why int more fast than uint ?
	if (minp >= 8192) {
		j = 1029, p = 8209;
	}

	while (p < minp) {
		NEXT_PRIME(p, j);
	}

	uint nextp = 0;
	ltype remp = 0;

	for (; p < sqrtp; NEXT_PRIME(p, j)) {
		if (p < sievesize / 2) {
			uint sieveindex = fastDiv(start, p);
			const FirstWheel cwn = FirstWheel30[sieveindex % WHEEL][WheelData30[p % WHEEL].Index];
			sieveindex += cwn.Correct * p;
			sieveSmall1(bitarray, bitarray + bytes, p, sieveindex, cwn.MultipleIndex);
		} else {
			PADD_DIFF(p);
			if (p > nextp)
				remp = start / (nextp = p + 20001);
			uint sieveindex = fastDiv(start - remp * p, p);
			if (sieveindex <= sievesize) {
				if (sieveindex % 2 == 0)
					sieveindex += p;
				if (sieveindex <= sievesize)
					bitarray[sieveindex / WHEEL] |= WheelData30[sieveindex % WHEEL].Mask;
			}
		}
	}
}

//core code of this algorithm for large range
//sieve prime multiples in [start, start + sievesize)
static void eratSieveMedium(uchar bitarray[], const ltype start, const uint sievesize, const uint minp, uint sqrtp)
{
	const uint bytes = sievesize / WHEEL + 1;
	if ((start + sievesize) < (ltype)sqrtp * sqrtp) {
		sqrtp = (uint)sqrt((double)start + sievesize) + 1;
	}

	int j = 8 + PRIME_PRODUCT / 9699690, p = Prime[0];//why int more fast than uint ?
	if (minp >= 131072) {
		j = 12252, p = 131101;
	} else if (minp >= 8192) {
		j = 1029, p = 8209;
	}

	while (p < minp) {
		p = MediumWheel[++j].p;
	}

	for (WheelPrime* pwheel = MediumWheel + j; p < sqrtp; p = pwheel->p) {
		if (pwheel->SieveIndex > (sievesize << 3)) {
			pwheel++->SieveIndex -= sievesize << 3;
			continue;
		}

		uint sieveindex = pwheel->SieveIndex >> 3;
		uint skipindex = pwheel->SieveIndex & 7;
		if (p > bytes / 2) {
			do {
				//250/1900 //150/1350
				bitarray[sieveindex / WHEEL] |= WheelData30[sieveindex % WHEEL].Mask;
				sieveindex += MultipleFactor30[skipindex++] * p;
			} while (sieveindex < sievesize);
		} else {
			//140/1900 ms //200/1350
			sieveSmall1(bitarray, bitarray + bytes, p, sieveindex, skipindex);
			sieveindex += ((sievesize - sieveindex) / (WHEEL * p)) * (WHEEL * p);
			do {
				sieveindex += MultipleFactor30[skipindex++] * p;
			} while (sieveindex < sievesize);
		}
		pwheel++->SieveIndex = (sieveindex - sievesize) << 3 | (skipindex & 7);
	}
}

/**
* This implementation uses a sieve array with 30 numbers per byte and
* a modulo 210 wheel that skips multiples of 2, 3, 5 and 7.
*/
static void eratSieveBucket(uchar bitarray[], const int bucketindex)
{
	_BucketStock* pbucket = &Bucket[bucketindex];
	uint wheelprimes = pbucket->WheelSize;
	pbucket->WheelSize = 0;

	for (int loops = wheelprimes % BLOCK_SIZE; wheelprimes > 0; loops = BLOCK_SIZE) {
		Stock* phead = StockList[StockSize++] = pbucket->Head;
		WheelPrime* pwheel = phead->WheelData;
		pbucket->Head = phead->NextStock;
		wheelprimes -= loops;

		while (loops--) {

			const uint p = pwheel->p;//int more fast on amd
			uint sieveindex = pwheel++->SieveIndex;

#if CPU != 0
			const uint wheeldata = *(uint*)&NextWheel210[sieveindex & 63][p & 63];
			bitarray[sieveindex >> 6] |= wheeldata >> 24;
			sieveindex = (uchar)(wheeldata >> 16) + (wheeldata & 255) * (p >> 6) + (sieveindex >> 6);
			pushBucket(sieveindex, p, wheeldata >> 8);
#else   //fast on core2 i5
			const WheelFactorization wheeldata = NextWheel210[sieveindex & 63][p & 63];
			bitarray[sieveindex >> 6] |= wheeldata.UnsetBit;
			sieveindex = wheeldata.Correct + wheeldata.MultipleIndex * (p >> 6) + (sieveindex >> 6);
			pushBucket(sieveindex, p, wheeldata.WheelIndex);
#endif
		}
	}
}

static int doSieveResult(uchar bitarray[], ltype start, const int bytes, CmdData* cmd)
{
	int primes = 0;
	if (cmd == NULL || cmd->Cmd == COUNTBITS) {
		primes = countZeroBitsArray((stype*)bitarray, bytes * 8);
	} else if (cmd->Cmd == SAVEPRIMEGAP) {
		primes = savePrimeGap((ushort*)bitarray, (bytes + 1) / 2, (pdtype*)cmd->Data + cmd->Primes);
	} else if (cmd->Cmd == SAVEPRIME) {
		primes = savePrime((ushort*)bitarray, start, bytes / 2, cmd->Data + cmd->Primes);
	} else if (cmd->Cmd == PRINTPRIME) {
		primes = printPrime((ushort*)bitarray, start, bytes / 2, cmd->Primes);
	} else if (cmd->Cmd == COPYBYBIT) {
		memcpy(cmd->Data, bitarray, bytes + 8);
		primes = bytes;
	} else if (cmd->Cmd == PRIMESUM) {
		primes = countPrimeSum((ushort*)bitarray, start, bytes / 2, (ltype*)cmd->Data);
	} else if (cmd->Cmd == FINDMAXGAP) {
		primes = findPrimeGap((ushort*)bitarray, start, (bytes + 1) / 2, cmd->Data);
	}

	if (cmd)
		cmd->Primes += primes;
	return primes;
}

static int segmentedSieve(const ltype start, const uint sievesize, uint wheeloffset, CmdData* cmd = NULL)
{
	uchar bitarray[MAX_CACHE + 64];

	const int bytes = preSieve(bitarray, start, sievesize);
	if (wheeloffset > 0) {
		memset(bitarray, -1u, wheeloffset / WHEEL);
		bitarray[wheeloffset / WHEEL] |= (1 << WheelData30[wheeloffset % WHEEL].Leng) - 1;
	}

	uint segsize = CpuCache.L1Size;
	#pragma omp parallel for num_threads(2)
	for (uint sieveindex = 0; sieveindex < sievesize; sieveindex += segsize) {
		if (segsize + sieveindex > sievesize)
			segsize = sievesize - sieveindex;
		eratSieveSmall(bitarray + sieveindex / WHEEL, start + sieveindex, segsize, CpuCache.L1Maxp);
	}

	const uint sqrtp = (uint)sqrt((double)start + sievesize) + 1;
	const uint maxp = sqrtp < Config.SieveSize / HALF_SEG ? sqrtp : Config.SieveSize / HALF_SEG;
	if (sievesize == CpuCache.L2Size || !BucketInfo.Active) {
		eratSieveMedium(bitarray, start, sievesize, CpuCache.L1Maxp, maxp);
	} else {
		segsize = CpuCache.L2Size;
		eratSieveMedium(bitarray, start, sievesize, CpuCache.L2Maxp, maxp);
		for (uint sieveindex = 0; sieveindex < sievesize; sieveindex += CpuCache.L2Size) {
			if (segsize + sieveindex > sievesize)
				segsize = sievesize - sieveindex;
			eratSieveMedium(bitarray + sieveindex / WHEEL, start + sieveindex,
					segsize, CpuCache.L1Maxp, CpuCache.L2Maxp);
		}
	}

	//bucket sieve = 256, 1051/1700[260:800:580], 1904, 2010 amd phoenm2 820
	//bucket sieve = 1024, 930/1180, 1324, 1360[130:700:500] intel i5 560m
	if (BucketInfo.Active) {
//		static double time_use = 0; double ts = getTime();
		eratSieveBucket(bitarray, BucketInfo.CurIndex & BucketInfo.BucketSize);
//		time_use += getTime() - ts;
//		if (sievesize != Config.SieveSize) { printf("eratSieveBucket time %.f ms\n", time_use); time_use = 0; }
		BucketInfo.CurIndex++;
	}

	return doSieveResult(bitarray, start, bytes, cmd);
}

//core code of this algorithm
//sieve prime multiples in [start, start + sievesize)
static int segmentedSieve(ltype start, int sievesize, CmdData* cmd = NULL)
{
	uchar bitarray[MAX_CACHE + 64];
	const uint bytes = preSieve(bitarray, start, sievesize);
	const uint sqrtp = (uint)sqrt((double)start + sievesize) + 1;
	sievesize += (int)(start % WHEEL);
	start -= start % WHEEL;
	for (uint sieveindex = 0, segsize = CpuCache.L1Size; sieveindex < sievesize; sieveindex += segsize) {
		if (segsize + sieveindex > sievesize)
			segsize = sievesize - sieveindex;
		eratSieveSmall(bitarray + sieveindex / WHEEL, start + sieveindex, segsize, CpuCache.L1Maxp);
	}

	eratSieveMedium2(bitarray, start, sievesize, CpuCache.L1Maxp, sqrtp);

	return doSieveResult(bitarray, start, bytes, cmd);
}

static int setSieveSize(int sievesize)
{
	if (sievesize <= 0) {
		memset(PiCache, 0, sizeof(PiCache));
		Config.ResultCached = false;
		return 0;
	}

	if (sievesize < 12) {
		sievesize = (WHEEL << 10) << sievesize;
	} if (sievesize < 2048) {
		sievesize *= (WHEEL << 10);
	}

	if (sievesize > MAX_CACHE * WHEEL)
		sievesize = MAX_CACHE * WHEEL - 8 * CpuCache.L1Size;

//	sievesize -= sievesize % (WHEEL210 * 8);
	sievesize = (1 << ilog2(sievesize / WHEEL + 1)) * WHEEL;

	if (sievesize != Config.SieveSize) {
		memset(PiCache, 0, sizeof(PiCache));
		Config.ResultCached = false;
	}

	return Config.SieveSize = sievesize;
}

static int checkSmall(const ltype start, const ltype end, ltype prime[], bool print = false)
{
	int primes = 0;
	const uchar smallp[] = {2, 3, 5, 7};
	for (int i = 0; i < sizeof(smallp) / sizeof(smallp[0]); i++) {
		if (start <= smallp[i] && smallp[i] <= end) {
			primes++;
			if (print)
				printf("%d %d\n", primes, (uint)smallp[i]);
			else if (prime) {
	//			prime[primes + 0] = prime[primes - 1];
	//			prime[primes - 1] = smallp[i];
			}
		}
	}

	return primes;
}

//calculate number of prime in Range[start, end] with start <= end
static ltype PI(ltype start, const ltype end, CmdData* cmd = NULL)
{
	assert (start <= end && start >= 0);

	ltype primes = 0;
	if (Config.SaveResult) {
		freopen("prime.txt", "w", stdout);
		CmdData cmdbuf = {PRINTPRIME, 0, NULL};
		cmd = &cmdbuf;
	}
	if (start <= PRIME7) {
		bool print = (cmd != NULL && cmd->Cmd == PRINTPRIME);
		if (cmd)
			primes = checkSmall(start, end, (ltype*)cmd->Data, print);
		else
			primes = checkSmall(start, end, NULL);
		if (cmd && (cmd->Cmd == SAVEPRIME || cmd->Cmd == PRINTPRIME))
			cmd->Primes += primes;
	}

	uint sievesize = Config.SieveSize;
	if (end - start <= sievesize) {
		primes += segmentedSieve(start, (int)(end - start) + 1, cmd);
		start += sievesize;
	} else {
		primes += segmentedSieve(start, sievesize - int(start % sievesize), cmd);
		primes += initPiCache(start / sievesize + 2, end / sievesize + 1, 1, cmd);
		primes += segmentedSieve(end - end % sievesize, int(end % sievesize) + 1, cmd);
	}

	if (Config.SaveResult) {
		freopen(CONSOLE, "w", stdout);
	}

	return primes;
}

static ltype PI2(ltype start, const ltype end, uint sievesize, CmdData* cmd)
{
	const double ts = getTime();
	const int modelstart = (int)(start % WHEEL210);
	start -= modelstart;
	const int segs = (end - start) / sievesize;
	ltype primes = segmentedSieve(start, sievesize, modelstart, cmd);
	start += sievesize;

	for (int ci = 1; start < end; start += sievesize) {
		if (start + sievesize > end)
			sievesize = end - start + (end & 1);
		const int seg2 = segmentedSieve(start, sievesize, 0, cmd);
		primes += seg2;
#if _DEBUG && 0
		const int seg = segmentedSieve(start, sievesize);
		if (seg != seg2) {
			printf("seg2 = %d != %d = seg\n", seg, seg2);
		}
#endif
		if ((ci++ & Config.Progress) == 1) {
			printf("%02d%%, sieve time ~ %.3f sec\r",
					ci * 100 / segs, (getTime() - ts) * segs / ci / 1000.0);
		}
	}

	return primes;
}

static void printPiResult(ltype start, ltype end, ltype primes, double ts)
{
	if (Config.ShowRet) {
		int sta10 = log((double)start + 1) / log(10.0) + 1e-4;
		int end10 = log((double)end + 1) / log(10.0) + 1e-4;
		int dif10 = log((double)(end - start) + 1) / log(10.0) + 1e-4;

		if (start > 0) {
			if (start % ipow(10, sta10) == 0)
				printf("PI[%de%d, ", (int)(start / ipow(10, sta10)), sta10);
			else
				printf("PI[%llu, ", (ltype)start);

			if (end % ipow(10, end10) == 0)
				printf("%de%d]", (int) (end / ipow(10, end10)), end10);
			else if ((end - start) % ipow(10, dif10) == 0) {
				if (start % ipow(10, sta10) == 0) {
					printf("%de%d+", (int)(start / ipow(10, sta10)), sta10);
					printf("%de%d]", (int)((end - start) / ipow(10, dif10)), dif10);
				} else
					printf("%llu+%de%d]", start, (int)((end - start) / ipow(10, dif10)), dif10);
			} else
				printf("%llu]", end);
		} else if (end % ipow(10, end10) == 0)
			printf("PI(%de%d)", (int) (end / ipow(10, end10)), end10);
		else
			printf("PI(%llu)", (ltype)end);

		printf(" = %llu", (ltype)primes);
		if (Config.ShowTime)
			printf(", time use %.3f sec", (getTime() - ts) / 1000.0);
		putchar('\n');
	}
}

ltype countPrime2(const ltype start, const ltype end, CmdData* cmd = NULL)
{
	assert (start <= end);

	ltype primes = 0;
	uint sievesize = Config.SieveSize;
	const uint sqrtp = (uint)sqrt((double)end) + 2;
	sievePrime(Prime + 1, sqrtp);
	double ts = getTime();

	//small range
	if (end - start < sievesize * 2) {
		primes = PI(start, end, cmd);
		printPiResult(start, end, primes, ts);
		return primes;
	}
	ts = getTime();

	primes = checkSmall(start, end, NULL);

	BucketInfo.Active = false;
	if (sqrtp > sievesize / HALF_SEG) {
		if (sievesize < CpuCache.L2Size && sqrtp > 1000000)
			sievesize = setSieveSize(SIEVE_SIZE);
		BucketInfo.Active = sqrtp > sievesize / HALF_SEG;
	}

	ltype new_start = start - start % WHEEL210;
	initMediumWheel(sievesize, sqrtp, new_start);

	if (BucketInfo.Active) {
		initBucketInfo(sievesize, sqrtp, end - new_start);
		initBucketStock(sievesize, sqrtp, new_start);
		if (Config.Showlog && 0) {
			printf("bucket info:\n");
			printf("	MaxIndex = %d, BucketSize = %d, Bucket Stock remain = %2.1f%%,"
					"\n	SieveSize = %d, First Bucket use = %2.1f%%\n",
					BucketInfo.MaxIndex, BucketInfo.BucketSize,
					StockSize * 100.0 / FREEBLK_SIZE, BucketInfo.SieveSize / 8 * WHEEL,
					Bucket[0].WheelSize * 100.0 / (BLOCK_SIZE * STOCK_SIZE));
			for (int i = 0; Bucket[i].WheelSize && 0; i++) {
				printf("%3d: wheelsize/stock = %7d -> [%3d]\n",
						i, Bucket[i].WheelSize, Bucket[i].WheelSize / BLOCK_SIZE);
			}
		}

		if (Config.Showlog) {
			printf("init bucket list time use %.3f sec and sieve size = %d k\n",
					(getTime() - ts) / 1000.0, sievesize / (WHEEL << 10));
		}
	}

	ts = getTime();
	primes += PI2(start, end, sievesize, cmd);
	
	if (!cmd || cmd->Cmd != FINDMAXGAP)
		printPiResult(start, end, primes, ts);

	if (BucketInfo.Active) {
		free(StockList[0]->WheelData);
		free(StockList[0]);
		BucketInfo.Active = false;
		Bucket[0].WheelSize = 0;
	}

	return primes;
}

ltype countPrime(const ltype start, const ltype end, CmdData* cmd = NULL)
{
	double ts = getTime();

	if (Config.Showlog) {
		printf("segment sievesize = %d k : %d \n",
				(Config.SieveSize / WHEEL) >> 10, Config.SieveSize);
	}

	sievePrime(Prime + 1, sqrt(double(end)) + 2);

	const ltype primes = PI(start, end, cmd);

	printPiResult(start, end, primes, ts);

	return primes;
}

static int showPrime(double ts, pdtype primegap[], const char* file)
{
	ts = getTime() - ts;

	int primes = 8 + PRIME_PRODUCT / 9699690;
	if (file)
		freopen(file, "w", stdout);

	uint p = Prime[0];
	if (Config.Showlog) {

		assert(primegap[1] + primegap[2] == 3);
#if 1
//		double sum = 0;
		for (; primegap[primes] > 0; primes++) {
			p += primegap[primes];
			PADD_DIFF(p);
#if 0
			//sum += 1.0 / p;
			if (primes < 25)
				printf("p = %d, sigma 1/p = %d\n", p, primegap[primes]);
#endif
		}

#endif

//		assert(PI(0, p + 1) == primes);
		printf("\nPrime[%u] = %u, ", primes, p);
//		printf("sum = %.6lf, ", sum);
		printf("save primes gap use %.3f sec\n", ts / 1000.0);
	}

	return primes;
}

static int sievePrime(pdtype primegap[], uint n)
{
	if (n <= Config.Maxp) {
		return 0;
	}
	Config.Maxp = n;

	double ts = getTime( );
	int primes = checkSmall(0, PRIME7, NULL);
	primegap[primes] = -PRIME7;
	CmdData cmd = {SAVEPRIMEGAP, primes, (ltype*)primegap};

	int sievesize = Config.SieveSize;
	Config.SieveSize = CpuCache.L1Size * 4;

//	initMediumWheel(0, sqrt(n) + 10, Config.SieveSize);
//	primes += PI2(0, n + 4 * WHEEL, sievesize, &cmd) - primes;
	primes += PI(0, n + 2 * WHEEL, &cmd) - primes;
	Config.SieveSize = sievesize;

	primegap[primes + 0] = 0;
	primegap[primes + 1] = primegap[primes + 2] = (pdtype)(1 << 31);

	showPrime(ts, Prime + 1, NULL);

	return primes;
}

//detect cpu L2 size
static int findBestSievesize( )
{
	int bestcache = 16;
	double mintime = 1E7;
	for (int sievesize = bestcache; sievesize < 512; sievesize *= 2) {
		double ts = getTime();
		setSieveSize(sievesize);
		Config.SieveSize = sievesize * WHEEL << 10;
		int primes = PI(0, 1E8);
		ts = getTime() - ts;
		if (ts < mintime) {
			mintime = ts;
			bestcache = sievesize;
		}
		if (Config.Showlog) {
			printf("cache %3dk: time %.2lf ms, primes = %d\n", sievesize, ts, primes);
		}
	}

	return bestcache;
}

static ltype initPiCache(ltype starti, ltype endi, int threads, CmdData* cmd)
{
	ltype pi = 0;
	const int sievesize = Config.SieveSize;

	for (ltype bi = starti; bi < endi; bi += threads) {
		if (bi < sizeof(PiCache) / sizeof(PiCache[0])) {
			if (PiCache[(int)bi] == 0 || cmd) {
				PiCache[(int)bi] = segmentedSieve((bi - 1) * sievesize, sievesize, cmd);
			}
			pi += PiCache[(int)bi];
		} else {
			pi += segmentedSieve((bi - 1) * sievesize, sievesize, cmd);
		}

		if (Config.Progress && (bi & 511) == 0) {
			printf("%2d%%\r", 100 * bi / endi);
		}
	}

	return pi;
}

//http://primes.utm.edu/nthprime/algorithm.php
static ltype findKthPrime(const ltype kth)
{
	double ts = getTime();

	if (kth <= 4) {
		uchar smallp[] = {0, 2, 3, 5, 7};
		printf("Prime[%llu] = %d\n", kth, (int)smallp[(int)kth]);
		return smallp[kth];
	}

	assert(Config.SieveSize % WHEEL == 0);
	/****
		Dusart99
		P. Dusart, "The kth prime is greater than k(ln k + ln ln k-1) for k> 2," Math. Comp.,
		n (log n + log log n - 1) < p(n) < n (log n + log log n - 0.9484) n > 39017
		(x/log x)(1 + 0.992/log x) < pi(x) <(x/log x)(1 + 1.2762/log x) x > 598
	*/
	uint bi = 0, ci = 0;
	const int smalltw = checkSmall(0, WHEEL, NULL);
	const ltype lowbound = kth * (log((double)kth) + log(log(double(kth))) - 1);
	const ltype upbound = lowbound + kth * (1 - 0.9484) + 2 * Config.SieveSize;

	if (kth > 1000000) {
		if (upbound / Config.SieveSize >= sizeof(PiCache) / sizeof(PiCache[0])) {
			return puts("number is too large");
		}
		bi = lowbound / Config.SieveSize;
		if (!Config.ResultCached || PiCache[2] >= PiCache[bi]) {
			startWorkThread(upbound, Config.Threads);
		}
		assert (kth > PiCache[bi] + smalltw);
	}

	if (!Config.ResultCached) {
		startWorkThread(upbound, Config.Threads);
	}

	for (; PiCache[bi + 1] && kth > PiCache[bi + 1] + smalltw; bi++) {

	}

	sievePrime(Prime + 1, (uint)sqrt((double)upbound));
	ushort* bitarray = (ushort*)malloc(MAX_CACHE + 8);
	bitarray[0] = COPYBYBIT;
//	segmentedSieve((ltype)(bi) * Config.SieveSize, Config.SieveSize, (uchar*)bitarray);

	ltype curth = PiCache[bi] + smalltw;
	while (true) {
		const int bit0 = 16 -
#if POPCNT == 0
			WordNumBit1[bitarray[ci++] | 0x0000ffff];
#else
		(Config.Popcnt ? _mm_popcnt_u32(bitarray[ci++]) : WordNumBit1[bitarray[ci++]]);
#endif
		if (curth + bit0 < kth) {
			curth += bit0;
		} else {
			ci--;
			break;
		}
	}

	ltype prime = (ltype)bi * Config.SieveSize + ci * WHEEL * 2;

	const ushort mask = bitarray[ci];
	for (int k = 0; k < 16; k++) {
		if (0 == (mask & (1 << k)) && (++curth == kth)) {
			prime += Pattern[k];
			break;
		}
	}

	free(bitarray);
	printf("Prime[%llu] = %llu, time use %.2lf ms\n",
			kth, (ltype)prime, (getTime() - ts));

	ltype pi2 = countPrime2(0, prime);
	if (pi2 != kth) {
		printf("error pi2 = %llu != %llu\n", pi2, kth);
	}

	return prime;
}

static void listFactor(const ltype dn, int count)
{
	for (int i = 0; i < count; i++) {
		ltype n = dn + i;

		if (n > 1)
			printf("%llu = ", n);

		for (uint j = 1, p = 2; p <= n / p && n > 1; p += Prime[++j]) {
			int factors = 0;
			while (n % p == 0) {
				n /= p;
				factors++;
			}

			if (factors > 0) {
				printf(" %d", p);
				if (factors > 1) {
					printf("^%d", factors);
				}
				if (n > 1) {
					printf(" *");
				}
			}
		}

		if (dn == n && dn > 1) {
			puts("is prime");
		} else if (n > 1) {
			printf(" %llu", n);
		}
		putchar('\n');
	}
}

//the sieve of Eratosthenes implementated by bit packing
//all prime less than 2^16 will be saved in prime buffer List
//Prime[0] is the first sieve prime, Prime[i] is the difference
//Prime[0] = 2, Prime[1] = 3 - 2, Prime[2] = 5 - 3;
//of the adjacent prime, Prime[i] = Prime[i] - Prime[i - 1];
static int eratoSieve(const uint maxp)
{
	int primes = 1;
	uchar bitarray[300006 >> 4];
	memset(bitarray, 0, sizeof(bitarray));

	Prime[0]= (PRIME_PRODUCT % 19 == 0) ? 23 : 19;
	uint lastp = Prime[1] = 2;

	for (uint p = 3; p <= maxp; p += 2) {
		if (0 == (bitarray[p >> 4] & MASK_N(p / 2))) {

			Prime[++primes] = p - lastp;
			lastp = p;
			if (p > (1 << 16)) {
				continue;
			}

			for (uint j = p * p / 2; j <= maxp / 2; j += p)
				bitarray[j >> 3] |= MASK_N(j);
		}
	}
//	assert(primes < sizeof(Prime) / sizeof(Prime[0]));

	//pack the last two byte for safety
	Prime[primes + 2] = Prime[primes + 1] = 255;

	return primes;
}

//The first presieved template
//sieve the first 8th prime multiples
static void initPreSieved( )
{
	const int smallprimes[ ] = {7, 11, 13, 17, 19, 23, 29};

	for (int i = 0; PRIME_PRODUCT % smallprimes[i] == 0; i++) {
		int start = smallprimes[i], p2 = 2 * smallprimes[i];
		for (; start < PRIME_PRODUCT; start += p2) {
			PreSieved[start / WHEEL] |= WheelData30[start % WHEEL].Mask;
		}
	}
}

//init WordNumBit1 table in 0-2^16
static void initBitTable( )
{
	int i;
	//1.
	if (Config.Popcnt == 0 || POPCNT == 0) {
		int nbitsize = sizeof(WordNumBit1) / sizeof(WordNumBit1[0]);
		WordNumBit1[0] = 0;
		for (i = 1; i < nbitsize; i++)
			WordNumBit1[i] = WordNumBit1[i >> 1] + (i & 1);
	}

	//2. init Left most bit table
	for (i = 0; i < (1 << 16); i += 2) {
		Lsb[i + 0] = Lsb[i >> 1] + 1;
		Lsb[i + 1] = 0;
	}

	for (i = 0; i < (1 << 16); i += 2) {
		Lsb[i + 0] = Pattern[Lsb[i]];
		Lsb[i + 1] = Pattern[0];
	}

	//4.
	for (uint k = 1; k < (1 << 16) - 1; k++) {
		int pattern = 0, bits = 0;
		for (int j = 0; j < 16; j++) {
			if (k & (1 << j)) {
				if (pattern == 0)
					WheelGap[k].Beg = Pattern[j];
				else
					WheelGap[k].Gap[bits - 1] = Pattern[j] - pattern;
				bits++;
				pattern = Pattern[j];
			}
		}
		WheelGap[k].End = 2 * WHEEL - pattern;
		WheelGap[k].Bits = bits;
	}
	//	WheelGap[0].Beg = WHEEL * 2; WheelGap[0].Bits = 0;
}

static void initWheel30()
{
	//wi/pi
	const uchar multipleIndex[][8] =
	{
		{0, 3, 2, 1, 6, 5, 4, 7},
		{1, 0, 4, 5, 2, 3, 7, 6},
		{2, 6, 0, 4, 3, 7, 1, 5},
		{3, 5, 6, 0, 7, 1, 2, 4},
		{4, 2, 1, 7, 0, 6, 5, 3},
		{5, 1, 7, 3, 4, 0, 6, 2},
		{6, 7, 3, 2, 5, 4, 0, 1},
		{7, 4, 5, 6, 1, 2, 3, 0},
	};

	for (int i = 0; i < WHEEL; i += 1) {
		for (int pi = 0; pi < 8; pi++) {
			int multiples = 0, sieveindex = i;
			if (i % 2 == 0) {
				multiples += 1;
				sieveindex += Pattern[pi];
			}
			int wi = WheelData30[sieveindex % WHEEL].Index;
			while (wi < 0) {
				sieveindex += Pattern[pi] * 2;
				wi = WheelData30[sieveindex % WHEEL].Index;
				multiples += 2;
			}
			FirstWheel firstWheel = {multipleIndex[wi][pi], wi, multiples};
			FirstWheel30[i][pi] = firstWheel;
		}
	}
}

static void initWheel210()
{
	const uchar Pattern210[ ] =
	{
		1, 11, 13, 17, 19, 23, 29, 31,
		37, 41, 43, 47, 53, 59, 61, 67,
		71, 73, 79, 83, 89, 97,101,103,
		107,109,113,121,127,131,137,139,
		143,149,151,157,163,167,169,173,
		179,181,187,191,193,197,199,209,
	};

	int wi = 0, i = 0;
	for (int j = 0; j < WHEEL210; j += 1) {
		WheelData210[j].Index = -1;
		if (WheelData30[j % WHEEL].Index >= 0 && j % 7 != 0)
			WheelData210[j].Index = wi++;
		WheelData210[j].Mask = WheelData30[j % WHEEL].Mask;
	}

	for (wi = 0; wi < 48; wi++) {
		for (int pi = 0; pi < 48; pi++) {
			int next = Pattern210[wi] + 2 * Pattern210[pi];
			int multiples = 2;
			while (WheelData210[next % WHEEL210].Index < 0) {
				next += Pattern210[pi] * 2;
				multiples += 2;
			}
			WheelFactorization* pdata = &NextWheel210[wi][pi];
			pdata->Correct = next / WHEEL - Pattern210[wi] / WHEEL;
			pdata->UnsetBit = WheelData210[Pattern210[wi]].Mask;
			pdata->WheelIndex = WheelData210[next % WHEEL210].Index;
			pdata->MultipleIndex = multiples * 7;
		}
	}

	for (i = 0; i < WHEEL210; i += 1) {
		for (int pi = 0; pi < 48; pi++) {
			int multiples = 0, next = i;
			if (i % 2 == 0) {
				multiples += 1;
				next += Pattern210[pi];
			}
			int wi = WheelData210[next % WHEEL210].Index;
			while (wi < 0) {
				next += Pattern210[pi] * 2;
				wi = WheelData210[next % WHEEL210].Index;
				multiples += 2;
			}
			FirstWheel210[i][pi].WheelIndex = wi;
			FirstWheel210[i][pi].Correct = multiples;
		}
	}

	for (wi = 0; wi < sizeof(MultipleFactor210); wi++) {
		if (wi % 48 != 47)
			MultipleFactor210[wi] = Pattern210[(wi + 1) % 48] - Pattern210[wi % 48];
		else
			MultipleFactor210[wi] = WHEEL210 + Pattern210[(wi + 1) % 48] - Pattern210[wi % 48];
	}

	for (wi = 0; wi < 48; wi++) {
		for (int pi = 0; pi < 48; pi++) {
			uchar skip[48], cwi = wi;
			for (i = 0; i < sizeof(skip); i++) {
				skip[i] = NextWheel210[cwi][pi].MultipleIndex / 7;
				cwi = NextWheel210[cwi][pi].WheelIndex;
			}

			for (i = 0; i < sizeof(skip) / sizeof(skip[0]); i++) {
				if (memcmp(skip, MultipleFactor210 + i, sizeof(skip)) == 0) {
					FirstWheel210[Pattern210[wi]][pi].MultipleIndex = i;
					break;
				}
			}
		}
	}
}

//init Prime, PreSieved and WordNumBit1 table
void initPrime(int sievesize)
{
	static bool initOnce = true;
	if (initOnce) {
		eratoSieve(101100);
		initPreSieved( );
		initBitTable( );
		initWheel30( );
		initWheel210( );
		initOnce = false;
	}
	setSieveSize(sievesize);
}

static void benchMarkTest( )
{
	uint primeCounts[] =
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
		//4118054813U,// pi(10^11)
		203280221, // pi(2^32)
		155428406,  // prime count of [10^12, 10^12 + 2^32]
		143482916,  // prime count of [10^13, 10^13 + 2^32]
		133235063,  // prime count of [10^14, 10^14 + 2^32]
		124350420,  // prime count of [10^15, 10^15 + 2^32]
		116578809,  // prime count of [10^16, 10^16 + 2^32]
		109726486,  // prime count of [10^17, 10^17 + 2^32]
		103626726,  // prime count of [10^18, 10^18 + 2^32]
		98169972,   // prime count of [10^19, 10^19 + 2^32]
	};

	double ts = getTime();
	Config.Showlog = Config.ShowRet = false;
	uint primes = 0;
	setSieveSize(CpuCache.L1Size * 8);
	Config.Progress = 0;
	for (int i = 1; i <= 10; i++) {
		primes = countPrime2(0, ipow(10, i));
		printf("pi(10^%02d) = %d\n", i, primes);
		assert(primes == primeCounts[i - 1]);
	}
	assert(primeCounts[10] == countPrime2(0, ipow(2, 32)));
	printf("pi(2^32)  = %d\n", primeCounts[10]);

	Config.Progress = 7;
	setSieveSize(L2_CACHE_SIZE * 4);
	for (int j = 12; j <= 19; j++) {
		ltype start = ipow(10, j);
		primes = countPrime2(start, start + ipow(2, 32));
		printf("pi(10^%d, 10^%d+2^32) = %d\n", j, j, primes);
		assert(primes == primeCounts[j - 1]);
	}

	printf("Time elapsed %.f sec\n", (getTime() - ts) / 1000);
	puts("All tests passed SUCCESSFULLY!");
//	Config.Showlog = Config.ShowRet = true;
}

//test case code, test data from third party
static int startRandTest(int tcases, int sievesize, int powbase)
{
	Config.Showlog = Config.ShowRet = false;

	printf("number case = %d, sievesize = %d, test file %s flag = %d\n",
			tcases, sievesize, TEST_FILE, powbase);

	if (powbase == 0) {
		if (!freopen(TEST_FILE, "r", stdin)) {
			puts("can not read test data file");
			Config.ShowRet = true;
			freopen(CONSOLE, "r", stdin);
			return -1;
		}
	} else {
		Config.Progress = 0;
		if (!freopen(TEST_FILE, "w", stdout)) {
			puts("can not write test data file");
			Config.ShowRet = true;
			freopen(CONSOLE, "w", stdout);
			return -2;
		}
	}

	srand(time(NULL));
	if (sievesize == 0) {
		setSieveSize(rand() % 257);
	} else if (sievesize < 6024) {
		setSieveSize(sievesize);
	} else {
		setSieveSize(32);
	}

	printf("sievesize = %d\n", Config.SieveSize);

	ltype maxn = (ltype)pow(10.0, powbase);
	if (maxn < 1E6) {
		maxn = 2e9;
	}

	sievePrime(Prime + 1, 100000008);

	double ts = getTime();
	const char* sformat1 = "%u PI[%u, %u] = %u\n";
	const char* sformat2 = "PI[%u, %u] = %u\n";
	const char* sformat3 = "PI(%u) = %u\n";

	if (sizeof(ltype) != sizeof(uint)) {
		sformat1 = "%u PI[%llu, %llu] = %u\n";
		sformat2 = "PI[%llu, %llu] = %u\n";
		sformat3 = "PI(%llu) = %u\n";
	}

	int failedcases = 0;
	for (int i = 1; i <= tcases; i++) {
		if (powbase) {
			ltype start = 0 + ((ltype)rand()) * rand() * rand() % maxn;
			ltype end = 1E6 + ((uint)rand() * rand()) % (1000000000);
			if (start > end)
				start ^= (end ^= (start ^= end));

			if (start < 10) {
				i--;
				continue;
			}

			if (sievesize == 0) {
				printf(sformat3, end, PI(0, end));
			} else {
				printf(sformat1, i, start, end, PI(start, end));
			}
		} else {
			char linebuf[420] = {0};
			ltype index, start = 0, end;
			int retInFile;
			gets(linebuf);
			if (sscanf(linebuf, sformat1, &index, &start, &end, &retInFile) != 4 &&
				sscanf(linebuf, sformat2, &start, &end, &retInFile) != 3 &&
				sscanf(linebuf, sformat3, &end, &retInFile) != 2) {
				printf("case %d with wrong data %s\n", i, linebuf);
				if (failedcases++ > 10)
					break;
			}
#if 1
			int primes = PI(start, end);
#else
			int primes = countPrime2(start, end);
#endif

			if (primes != retInFile) {
				printf(sformat1, i, start, end, primes);
				printf("ecprime %llu -n%llu\n", (ltype)end, (ltype)start);
			}
//			printf("case %d\r", i);
			if ((i & 255) == 0)
				printf("case pass %d%%\r", i * 100 / tcases);
		}
	}

	freopen(CONSOLE, "r", stdin);
	freopen(CONSOLE, "w", stdout);
	Config.ShowRet = true;
	printf("test case time use %.2lf ms\n", getTime() - ts);

	return 0;
}

//make sure the input value is valid
static void listPrime(char params[][60])
{
	ltype start = 0, end = 1e9, step = 1e6;
	ltype buf[8] = {0, start, end, step, 0};
	int ni = 1;

	double ts = getTime();

	for (int i = 1; params[i][0] && ni < sizeof(buf) / sizeof(buf[0]); i++) {
		if (isdigit(params[i][0]) || params[i][0] == 'e')
			buf[ni++] = atoint64(params[i]);
	}

	start = buf[1], end = buf[2], step = buf[3];

	if (start > end && end < 10000000)
		end = start + end * step;
	if (step <= 0)
		step = 1000000;
	printf("[%llu:%llu:%llu]\n", (ltype)start, (ltype)end, (ltype)step);

	if (Config.SaveResult) {
		Config.ShowTime = false;
		freopen(TEST_FILE, "w", stdout);
	} else if (step < 1000) {
		Config.ShowTime = false;
	}

	if (start <= end && end - start >= step - 1) {
		int flag = (int)buf[4];
		for (ltype i = start; i <= end - step + 1 && i >= start; i += step) {
			ltype tend = i + step - 1;
			if (flag == 0) {
				countPrime(i, tend);
			} else if (flag == 1) {
				countPrime(start, tend);
			} else {
				countPrime(0, tend);
			}
		}
	}

	if (Config.SaveResult) {
		Config.SaveResult = false;
		Config.ShowTime = true;
		freopen(CONSOLE, "w", stdout);
	}

	printf("time use %.2lf\n", getTime() - ts);
}

static void getCpuid(int cpuinfo[4], int id)
{
#if _MSC_VER > 1300
	__cpuid(cpuinfo, id);
#elif _MSC_VER == 1200
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
#else
	int deax, debx, decx, dedx;
	__asm
	(
		"cpuid\n"
		:"=a"(deax),"=b"(debx),"=c"(decx),"=d"(dedx)
		:"a"(id)
	);
	cpuinfo[0] = deax;
	cpuinfo[1] = debx;
	cpuinfo[2] = decx;
	cpuinfo[3] = dedx;
#endif
}

// http://msdn.microsoft.com/en-us/library/hskdteyh%28v=vs.80%29.aspx
static int getCpuInfo()
{
	char cpuName[255] = {-1};
	int (*pTmp)[4] = (int(*)[4])cpuName;
	getCpuid(*pTmp++, 0x80000002);
	getCpuid(*pTmp++, 0x80000003);
	getCpuid(*pTmp++, 0x80000004);

	for (int i = 0; cpuName[i]; i++) {
		if (cpuName[i] != ' ' || cpuName[i + 1] != ' ')
			putchar(cpuName[i]);
	}

	int cpuinfo[4];
	getCpuid(cpuinfo, 0x80000006);
	printf(", L2 cache = %d kB\n", cpuinfo[2] >> 16);

	//amd cpu
	if (cpuName[0] == 'A') {
		CpuCache.L1Size = 64 * (WHEEL << 10);
	} else{
		CpuCache.L1Size = 32 * (WHEEL << 10);
	}
	CpuCache.L1Maxp = (CpuCache.L1Size / WHEEL) / L1_SIEVE_SEG;
	return cpuinfo[2] >> 16 ;
}

#ifdef _WIN32
static LONG WINAPI filterFunc(DWORD dwExceptionCode)
{
	return ((dwExceptionCode == STATUS_ILLEGAL_INSTRUCTION)
			? EXCEPTION_EXECUTE_HANDLER : EXCEPTION_CONTINUE_SEARCH);
}

static int getSystemInfo( )
{
	SYSTEM_INFO si;
	GetSystemInfo(&si);

	if (2 * si.dwNumberOfProcessors > Config.Threads) {
		Config.Threads = 2 * si.dwNumberOfProcessors;
	}

#if (_MSC_VER && POPCNT)
	__try {
		Config.Popcnt = _mm_popcnt_u32(1);
	}
	__except (filterFunc(GetExceptionCode())) {
		Config.Popcnt = 0;
	}
#elif POPCNT
	try {
		Config.Popcnt = _mm_popcnt_u32(1);
	}
	catch (...) {
		Config.Popcnt = 0;
	}
#endif

	if (si.wProcessorArchitecture == PROCESSOR_ARCHITECTURE_INTEL)
		printf("Cpu arch = x86, ");
#if _MSC_VER > 1400 || __GNUC__ > 3
	else if (si.wProcessorArchitecture == PROCESSOR_ARCHITECTURE_AMD64)
		printf("Cpu arch = x64, ");
#endif

	printf("cores = %ld, SSE4_Popcnt = %d\n",
			si.dwNumberOfProcessors, Config.Popcnt);

	return si.dwNumberOfProcessors;
}
#endif

static void printInfo(int argc)
{
	const char* line =
		"--------------------------------------------------------------------";
	puts(line);
	printf("Count/Sieve number of primes (0, 1E19 + 1E13), version %s\n", VERSION);
	puts("Implemented by the segmented sieve of eratosthenes [wheel = 30/210]");
	puts("Copyright @ by Huang Yuanbing 2011 - 2012 bailuzhou@163.com");

	puts(line);
	puts(line);

#ifdef _MSC_VER
	printf("Compiled by MS/vc++ %d", _MSC_VER);
#else
	printf("Compiled by GNU/g++ %d.%d.%d",
			__GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#endif

#if x86_64
	printf(" on x64 bit");
#endif

	printf(" on %s %s\n", __TIME__, __DATE__);

#ifdef _WIN32
	Config.Threads = getSystemInfo();
#else
	Config.Threads = sysconf(_SC_NPROCESSORS_CONF);
#endif

	getCpuInfo();

	if (argc > 0) {
		puts(line);
		printf("[MARCO] : ASM_X86 = %d, LIANGBCH = %d\n", ASM_X86, LIANGBCH);
		printf("[MARCO] : SIEVE_SIZE = %dk, BLOCK_SIZE = %dk, PDIFF = %u\n",
			SIEVE_SIZE, BLOCK_SIZE >> 7, MAX_PDIFF);
	}
	puts(line);
	puts(line);
}

static void doCompile()
{
	char programming[255];
	strcpy(programming, __FILE__);
	char* pdot = strchr(programming, '.');
	if (pdot) {
		strcpy(pdot, "_.exe");
		puts(programming);
	}

	const char* const cxxflag =
#if _MSC_VER
		"cl /O2 /Oi /Ot /Oy /GT /GL %s %s";
#elif defined X86_64
		"g++ -m64 -msse3 -mpopcnt -march=native -O3 -funroll-loops -s -pipe %s -o %s";
#else
		"g++ -m32 -msse3 -mpopcnt -march=native -O3 -funroll-loops -s -pipe %s -o %s";
#endif

	char commpile[255] = {0};
	sprintf(commpile, cxxflag, __FILE__, programming);
	puts(commpile);
	system(commpile);
	system(programming);
}

//get the first digit number index
static int parseConfig(const char params[][60])
{
	int cmdi = -1;
	int cdata = 0;

	for (int i = 0; params[i][0]; i++) {
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
		//	case 'P':
		//		Config.ShowTime = !Config.ShowTime;
		//		break;
			case 'G':
				Config.Progress = (1 << cdata) - 1;
				if (Config.Progress == 0)
					Config.Progress = (1 << 30) - 1;
				break;
			case 'D':
				Config.Showlog = !Config.Showlog;
				break;
			case 'F':
				Config.SaveResult = !Config.SaveResult;
				break;
			case 'V':
				printf("version %s\n", VERSION);
				break;
			case 'A':
				Config.CheckResult = !Config.CheckResult;
				break;
			case 'S':
				setSieveSize(cdata);
				break;
			case 'C':
				if (cdata > 8 && cdata < 512 && ((cdata - 1) & cdata) == 0) {
					CpuCache.L1Size = cdata * (WHEEL << 10);
					CpuCache.L1Maxp = (cdata << 10) / L1_SIEVE_SEG;
				}
				break;
			case 'T':
				if (cdata < 64)
					Config.Threads = cdata;
				break;
			case 'M':
				doCompile();
				break;
			default:
				cmdi = i;
				break;
		}
	}

	return cmdi;
}

//split ccmd string to params array
static int splitCmd(const char* ccmd, char cmdparams[64][60])
{
	int nwords = 0;

	for (int i = 0; i < 64; i++) {
		while (isspace(*ccmd) || ',' == *ccmd) {
			ccmd++;
		}
		if (*ccmd == 0 || *ccmd == ';') {
			break;
		}
		char* pc = cmdparams[i];
		char c = *ccmd;
		bool isvalid = false;
		while (isalnum(c) || c == '^' ||
				c == '+' || c == '-' || c == '*') {
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

bool excuteCommand(const char* cmd)
{
	while (cmd) {

		char params[64][60] = {0};
		char* pcmd = (char*) strchr(cmd, ';');
		if (splitCmd(cmd, params) <= 0)
			return false;

		int cmdi = parseConfig(params);

		if (cmdi == -1) {
			return true;
		}

		char cmdc = toupper(params[cmdi][0]);
		ltype start = atoint64(params[cmdi]);
		ltype end = atoint64(params[cmdi + 1]);
		if (!isdigit(cmdc) && cmdc != 'E') {
			start = atoint64(params[cmdi + 1]);
			end = atoint64(params[cmdi + 2]);
		}

		if (end == 0)
			end = start, start = 0;
		else if (end < start)
			end += start;

		if (cmdc == 'H') {
			puts(Help);
		} else if (cmdc == 'B') {
			puts("---------------------start benchmark------------------------");
			if (isdigit(params[cmdi + 1][0])) {
				excuteCommand("1e12 e10; 1e13 1e10; 1e14 1e10; 1e15 1e10");
				excuteCommand("1e16 e10; 1e17 1e10; 1e18 1e10; 1e19 1e10");
			}
			benchMarkTest();
		} else if (cmdc == 'U') {
			puts("---------------------start unit test------------------------");
			int sievesize = 0, powbase = 0;
			int testcase = atoi(params[cmdi + 1]);
			if (isdigit(params[cmdi + 2][0])) {
				sievesize = (int)atoint64(params[cmdi + 2]);
				if (isdigit(params[cmdi + 3][0]))
					powbase = (int)atoint64(params[cmdi + 3]);
			}
//			for (int i = 0; i < 10; i++)
			startRandTest(testcase, sievesize, powbase);
		} else if (cmdc == 'L') {
			puts("------------------start list multi result----------------");
			listPrime(params);
		} else if (cmdc == 'K') {
			puts("-------------------start find kth prime------------------");
			ltype kth = start;
			findKthPrime(kth);
			if (isdigit(params[2][0])) {
				kth = end;
				for (int i = 0; i < kth; i++)
					findKthPrime((uint)rand() * rand() % 100000000 + 1);
			}
		} else if (cmdc == 'P') {
			puts("-------------------start print prime---------------------");
			CmdData cmd = {PRINTPRIME, 0, NULL};
			//PI(start, end, &cmd);
			countPrime2(start, end, &cmd);
		} else if (cmdc == 'O') {
			puts("-------------------start factor prime--------------------");
			if (end == start)
				end += 1;
			listFactor(start, end - start);
		} else if (cmdc == 'Q') {
			return false;
		} else if (cmdc == 'Z') {
			puts("-------------------start save prime----------------------");
			CmdData cmd = {SAVEPRIME, 0, (ltype*)malloc(sizeof(ltype) * ((end - start) / 12) + 10000)};
			int primes = PI(start, end, &cmd);
			printf("last prime[%d] = %llu\n", primes, cmd.Data[primes - 1]);
			free(cmd.Data);
		} else if (cmdc == 'Y') {
			//1425172824437699411 1476
			//http://www.ieeta.pt/~tos/gaps.html
			puts("-------------------start find max gap -------------------");
			ltype data[4] = {0};
			//start = (4100000 + rand()* rand() % 2000000) * 1e12;
			//end = start + atoint64("1e10");
			CmdData cmd = {FINDMAXGAP, 0, data};
			countPrime2(start, end, &cmd);
			printf("maxp prime gap = %d on %llu\n", (int)cmd.Data[1], cmd.Data[2] - cmd.Data[1]);
		} else if (cmdi >= 0 && end > 0) {
			puts("-------------------start count primes -------------------");
			countPrime2(start, end);
			if (Config.CheckResult)
				countPrime(start, end);
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
	if (argc < 2) {
		printInfo(argc);
//		puts(Help);
	}

	initPrime(1024);

	char ccmd[1023] = {0};
	for (int i = 1; i < argc; i++) {
		strcat(ccmd, argv[i]);
		strcat(ccmd, " ");
	}
	if (argc > 1)
		excuteCommand(ccmd);

	excuteCommand("1e16 1e9 s10");
//	excuteCommand("y 1e18 1e9");
//	excuteCommand("y 1e18 1e9");

	while (true) {
		printf("\n[command or number] : ");
		if (!gets(ccmd) || !excuteCommand(ccmd))
			break;
	}

	return 0;
}

/***
OS: windows 7 32 bit
Mingw: gcc 4.6.3
CPU: Intel core i5 560m 2.66G (L1 32k, L2 256k, L3 3M)
CXXFLAGS: -Ofast -msse4 -s -pipe -mtune=native -march=corei7 -fomit-frame-pointer

range                            primesieve   ecprime    primenumber     Oliveira(not adding init time)
[1E10, 1E10+1E10] = 427154205    3.65         3.52       3.09
[1E11, 1E11+1E10] = 394050419    5.04         5.44       3.85
[1E12, 1E12+1E10] = 361840208    6.15         6.38       4.92             6.0
[1E13, 1E13+1E10] = 334067230    7.78         8.25       6.41             7.4
[1E14, 1E14+1E10] = 310208140    9.68         10.5       7.93             8.6
[1E15, 1E15+1E10] = 289531946    10.9         12.8       9.42             9.9
[1E16, 1E16+1E10] = 271425366    13.0         14.9       11.2/10.8        11.2
[1E17, 1E17+1E10] = 255481287    15.0         17.8       13.3/12.3        12.6
[1E18, 1E18+1E10] = 241272176    20.1         21.9       16.0/13.6        14.0
[1E19, 1E19+1E10] = 228568014    31.4         31.7       21.6/15.0        15
[1E19, 1E19+1E11] = 2285693139   184          197        163./156.        159

[1E18, 1E18+1E7 ] = 241295       1.64         2.91       0.74
[1E18, 1E18+1E8 ] = 2414886      2.31         3.10       1.67
[1E18, 1E18+1E9 ] = 24127085     4.78         4.71       3.12
[1E19, 1E19+1E9 ] = 22854258     10.7         12.5       6.42
[1E18, 1E18+1E11] = 2412731214   168.         169        138.             140

benchmark for pi(10^n, 10^n + 10^9) by two programming
based on windows 7 32 bit, amd Phoenm 2 X4 820 2.8G
                                 primesieve  ktprime
PI[1E11, 1E11+1E9] = 39475591    0.41/0.64   0.58/0.58
PI[1E12, 1E12+1E9] = 36190991    0.61/0.85   0.73/0.75
PI[1E13, 1E13+1E9] = 33405006    0.83/1.08   0.92/0.91
PI[1E14, 1E14+1E9] = 31019409    1.04/1.36   1.28/1.12
PI[1E15, 1E15+1E9] = 28946421    1.24/1.53   1.49/1.40
PI[1E16, 1E16+1E9] = 27153205    1.54/1.72   1.64/1.65
PI[1E17, 1E17+1E9] = 25549226    1.85/2.06   1.78/1.84
PI[1E18, 1E18+1E9] = 24127085    2.22/2.15   1.88/1.98
PI[1E19, 1E19+1E9] = 22854258    2.32/2.35   1.97/2.10
PI[1E18, 1E18+1E10]= 241272176   24.1/       24.4/22.0
PI[1E18, 1E18+1E11]= 2412731214  228./212    218./222
PI[1E18, 1E18+1E12]= 24127637783 2270        2142/

vc++
	cl /O2 /Oi /Ot /Oy /GT /GL PrimeNumber.cpp
mingw
	g++ -mpopcnt -march=native -O3 -msse3 -s -pipe PrimeNumber.cpp -o PrimeNumber
gc++
if (cpu spupport popcnt instruction and gcc version > 4.3)
	g++ -mpopcnt -march=native -O3 -s -pipe -lpthread PrimeNumber.cpp -o PrimeNumber
else
	g++ -O3 -s -pipe -lpthread PrimeNumber.cpp -o PrimeNumber
***/

