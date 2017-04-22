/***
the most fast segmented sieving of twin prime before 2018
doc:
	http://sweet.ua.pt/tos/software/prime_sieve.html
	http://primesieve.org/
***/
const char* Benchmark =
"Mingw: g++ 5.1.0 bailuzhou@163\n"
":g++ -DSIEVE_SIZE=1024 -march=native -funroll-loops -O3 -s -pipe\n"
"Windows  10 x64  on x86_64 cpu  i3 350M, i5 3470, i7 7500u\n"
"pi2(1e11, 1e10)  =  20498568    4.25     2.25     1.90\n"
"pi2(1e12, 1e10)  =  17278660    5.01     2.50     2.30\n"
"pi2(1e13, 1e10)  =  14735239    5.86     3.00     2.56\n"
"pi2(1e14, 1e10)  =  12706059    7.28     3.62     3.20\n"
"pi2(1e15, 1e10)  =  11069670    8.63     4.25     3.73\n"
"pi2(1e16, 1e10)  =  9728024     10.1     4.93     4.32\n"
"pi2(1e17, 1e10)  =  8614943     11.9     5.72     4.91\n"
"pi2(1e18, 1e10)  =  7687050     15.1     7.05     5.95\n"
"pi2(1e19, 1e10)  =  6895846     20.7     10.0     8.60\n"
"pi2(2^64-1e9,1e9)=  670362      7.02     4.00     4.12\n"
"pi2(1e18, 1e6)   =  794         0.75     0.46     0.37\n"
"pi2(1e18, 1e8)   =  77036       1.30     0.80     0.70\n"
"pi2(1e18, 1e9)   =  769103      3.20     1.60     1.41\n"
"pi2(1e14, 1e12)  =  1270127074  725      354      244\n"
"pi2(1e16, 1e12)  =  972773783   980      490      362\n"
"pi2(1e18, 1e12)  =  768599834   1240     600      492\n"
"pi2(1e19, 1e12)  =  689816098   1300     670      522\n";

#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <assert.h>
#include <time.h>
#include <string.h>

//const
enum EBASIC
{
	WHEEL30       = 30,
	WHEEL210      = 210,
	PRIME_PRODUCT = 210 * 11 * 13 * 17 * 19,
	FIRST_INDEX   = PRIME_PRODUCT / 9699690 + 7,
	WHEEL_SKIP    = 0x799b799b,
	PATTERN_GAP   = 28,
	PRIME_SIZE    = 6542 + 10,//pi(2^16)
	MAX_SEGMENT   = (1 << 10) * 4,//max Level Cache

	//performance marco
	L1_DCACHE_SIZE = 64,
	SI_BIT         = 8, //8 - 10
	WP_BIT         = 6, //6 - SI_BIT
#if (L2_DCACHE_SIZE < 256 || L2_DCACHE_SIZE > 4096)
	L2_DCACHE_SIZE =  256,
#endif
#ifndef SIEVE_SIZE
	SIEVE_SIZE     =  2048,
#endif
#ifndef UINT_MAX
	UINT_MAX       = -1u,
#endif

	ERAT_SMALL     = 8, //4 - 16
	ERAT_MEDIUM    = 2, //2 - 6
};

#ifndef PRIME_GAP
#define PRIME_GAP 2
#endif

enum EBUCKET
{
	UINT_PIMAX = 203280221, //= pi(2^32)
	MAX_BUCKET = 5464 * 3, //= 28*2^32 / 2^18 * 30 + 4
	WHEEL_SIZE = 1 << 12, //=4096  11: 16k, [10 - 13]
	MEM_WHEEL  = WHEEL_SIZE * 8, //=32768
	MEM_BLOCK  = (1 << 19) / WHEEL_SIZE, //=256   1 << 20:8 MB
	MAX_STOCK  = UINT_PIMAX / WHEEL_SIZE + MAX_BUCKET + MEM_BLOCK, //=55221
	MAX_POOL   = UINT_PIMAX / (MEM_BLOCK * WHEEL_SIZE) + 100, //=487
};

enum EFLAG
{
	SLOW_TEST = 1 << ('A' - 'A'),
	PRINT_RET = 1 << ('R' - 'A'),
	PRINT_TIME= 1 << ('T' - 'A'),
	SAVE_DATA = 1 << ('F' - 'A'),
};

enum ECMD
{
	COUNT_PRIME = 1,
	COPY_BITS,
	SAVE_PRIME,
	SAVE_BYTE,
	PCALL_BACK
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
# define ERAT_BIG         6 //2 - 6
#endif

# define WHEEL        WHEEL210
# define WHEEL_MAP    Wheel210
# define WHEEL_INIT   WheelInit210
# define WHEEL_FIRST  WheelFirst210
# define PATTERN      Pattern210

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
	[C: Cpu L1/L2 data cache size (L1:16-64, L2:256-2048)]\n\
	[S: Set sieve segment size (32 - 4096)k]\n\
	[L: Set sieve cache segs L(1-12)1, L(1-8)2 L(2-6)3]\n\
	[I: Info of programming]\n\
	[Z: Compile self and run]\n\
	[P: Print prime in [start, end]]\n\n\
Example:\n\
	1e16 1e16+10^10 s1024\n\
	p 1e12+100 100";

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
	WheelPrime* Wheel; //read only
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
static Stock* StockHead;

static WheelPrime* WheelPool [MAX_POOL]; //2G vm
static uint Prime[PRIME_SIZE];
static WheelPrime* MediumWheel;

//presieved small prime number <= 19.
static uchar PreSieved[PRIME_PRODUCT / WHEEL30];

#if BIT_SCANF == 0
static uchar Lsb[1 << 16];
#endif

//number of bits 1 binary representation table in Range[0-2^16)
static uchar WordNumBit0[1 << 16];

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

static uchar TwBits[6] =
{
#if PRIME_GAP == 2
	1, 3, 6, 9, 11, 14 //PRIME_GAP = 2
#else
	0, 2, 4, 8, 10, 12 //PRIME_GAP = 4
#endif
};

static uchar Pattern30[64] =
{
	1, 7, 11, 13, 17, 19, 23, 29,
	31, 37, 41, 43, 47, 49, 53, 59, 61
};

static int Pattern210[WHEEL210];
//api
#ifndef PRIME_LIB
struct Cmd
{
	int Oper;
	uchar* Data;
	uint64 Primes;
};

typedef void (*sieve_call)(uint64, uint64);
void initPrime(int sieve_size);
bool executeCmd(const char* cmd);
uint64 doSieve(const uint64 start, const uint64 end, Cmd* cmd);
int setSieveSize(uint sieve_size);
void setLevelSegs(uint level);
void setCpuCache(uint cpusize);
#else
#include "PrimeNumber.h"
#endif

static int64 getTime( )
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

//30% fast on old cpu pentium 4, fam3
inline static void crossOffWheelFactor2(uchar* p, const uchar* pend, const uint step)
{
	const uint o = step / WHEEL30;
	switch (WheelInit30[step % WHEEL30].PrimeIndex) {
	case 0 :
		while (p <= pend) {
			p[o *  0 +  0] |= BIT0, p[o *  6 +  0] |= BIT1;
			p[o * 10 +  0] |= BIT2, p[o * 12 +  0] |= BIT3;
			p[o * 16 +  0] |= BIT4, p[o * 18 +  0] |= BIT5;
			p[o * 22 +  0] |= BIT6, p[o * 28 +  0] |= BIT7;
			p += step;
		}
		break;
	case 1 :
		while (p <= pend) {
			p[o *  0 +  0] |= BIT0, p[o *  4 +  0] |= BIT7;
			p[o *  6 +  1] |= BIT3, p[o * 10 +  2] |= BIT2;
			p[o * 16 +  3] |= BIT6, p[o * 18 +  4] |= BIT1;
			p[o * 24 +  5] |= BIT5, p[o * 28 +  6] |= BIT4;
			p += step;
		}
		break;
	case 2 :
		while (p <= pend) {
			p[o *  0 +  0] |= BIT0, p[o *  2 +  0] |= BIT6;
			p[o *  6 +  2] |= BIT1, p[o *  8 +  2] |= BIT7;
			p[o * 12 +  4] |= BIT3, p[o * 18 +  6] |= BIT5;
			p[o * 20 +  7] |= BIT2, p[o * 26 +  9] |= BIT4;
			p += step;
		}
		break;
	case 3 :
		while (p <= pend) {
			p[o *  0 +  0] |= BIT0, p[o *  4 +  1] |= BIT6;
			p[o *  6 +  2] |= BIT5, p[o * 10 +  4] |= BIT2;
			p[o * 12 +  5] |= BIT1, p[o * 16 +  6] |= BIT7;
			p[o * 22 +  9] |= BIT4, p[o * 24 + 10] |= BIT3;
			p += step;
		}
		break;
	case 4:
		while (p <= pend) {
			p[o *  0 +  0] |= BIT0, p[o *  6 +  3] |= BIT3;
			p[o *  8 +  4] |= BIT4, p[o * 14 +  7] |= BIT7;
			p[o * 18 + 10] |= BIT1, p[o * 20 + 11] |= BIT2;
			p[o * 24 + 13] |= BIT5, p[o * 26 + 14] |= BIT6;
			p += step;
		}
		break;
	case 5 :
		while (p <= pend) {
			p[o *  0 +  0] |= BIT0, p[o *  4 +  2] |= BIT4;
			p[o * 10 +  6] |= BIT2, p[o * 12 +  7] |= BIT5;
			p[o * 18 + 11] |= BIT3, p[o * 22 + 13] |= BIT7;
			p[o * 24 + 15] |= BIT1, p[o * 28 + 17] |= BIT6;
			p += step;
		}
		break;
	case 6 :
		while (p <= pend) {
			p[o *  0 +  0] |= BIT0, p[o *  2 +  1] |= BIT4;
			p[o *  6 +  4] |= BIT5, p[o * 12 +  9] |= BIT1;
			p[o * 14 + 10] |= BIT6, p[o * 20 + 15] |= BIT2;
			p[o * 24 + 18] |= BIT3, p[o * 26 + 19] |= BIT7;
			p += step;
		}
		break;
	case 7 :
		while (p <= pend) {
			p[o *  0 +  0] |= BIT0, p[o *  2 +  1] |= BIT7;
			p[o *  8 +  7] |= BIT6, p[o * 12 + 11] |= BIT5;
			p[o * 14 + 13] |= BIT4, p[o * 18 + 17] |= BIT3;
			p[o * 20 + 19] |= BIT2, p[o * 24 + 23] |= BIT1;
			p += step;
		}
		break;
	}
}

//fast on new cpu
inline static void crossOffWheelFactor(uchar* ppbeg[], const uchar* pend, const uint p)
{
	#define CMPEQ_OR(n) if (ps##n <= pend) *ps##n |= BIT##n

	uchar* ps0 = ppbeg[0], *ps1 = ppbeg[1];
	uchar* ps2 = ppbeg[2], *ps3 = ppbeg[3];
	while (ps3 <= pend) {
#if PRIME_GAP == 2
		*ps0 |= BIT0, ps0 += p;
#endif

#if PRIME_GAP == 4
		*ps1 |= BIT1, ps1 += p;
#endif
		*ps2 |= BIT2, ps2 += p;
		*ps3 |= BIT3, ps3 += p;
	}
	CMPEQ_OR(0);
	CMPEQ_OR(1);
 	CMPEQ_OR(2);

	uchar* ps4 = ppbeg[4], *ps5 = ppbeg[5];
	uchar* ps6 = ppbeg[6], *ps7 = ppbeg[7];
	while (ps7 <= pend) {
		*ps4 |= BIT4, ps4 += p;
		*ps5 |= BIT5, ps5 += p;
#if PRIME_GAP == 4
		*ps6 |= BIT6, ps6 += p;
#endif
#if PRIME_GAP == 2
		*ps7 |= BIT7,
#endif
		ps7 += p;
	}
	CMPEQ_OR(4);
	CMPEQ_OR(5);
	CMPEQ_OR(6);
}

inline static void crossOff2Factor(uchar* ps0, uchar* ps1, const uchar* pend, const ushort smask, const uint p)
{
	const uchar masks1 = smask >> 8, masks0 = (uchar)smask;
	while (ps1 <= pend) {
		*ps1 |= masks1, ps1 += p;
		*ps0 |= masks0, ps0 += p;
	}
	if (ps0 <= pend)
		*ps0 |= masks0;
}

inline static int countBitsTable(const uint64 n)
{
	int sum = 0;
	const uint hig = (uint)(n >> 32), low = (uint)n;
	sum =  WordNumBit0[(ushort)hig] + WordNumBit0[(ushort)low];
	sum += WordNumBit0[hig >> 16] + WordNumBit0[low >> 16];
	return sum;
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

static void printPrime(uint64 index, uint64 prime)
{
	printf("(%llu, %llu)\n", prime, prime + PRIME_GAP);
}

static int doCall(const ushort bitarray[], uint64 sieve_index, const int size, uint64 sum_prime, sieve_call func)
{
	int primes = 0;
	for (int bi = 0; bi <= size / 2; bi++) {
		ushort mask = (bitarray[bi] >> 1) | (bitarray[bi + 1] << 15);
		for (int i = 0; i < sizeof(TwBits); i++) {
			const uint bits = TwBits[i];
			if ((mask >> bits) % 4 == 0) {
				func(++primes + sum_prime, sieve_index + Pattern30[(1 + bits) % 8] + (1 + bits) / 8 * WHEEL30);
			}
		}
		sieve_index += WHEEL30 * sizeof(mask);
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
	WheelPrime *pwheel = (WheelPrime*) malloc((blocks * MEM_BLOCK + 1) * MEM_WHEEL);
	WheelPool[BucketInfo.PoolSize ++] = pwheel;
	//assert (BucketInfo.PoolSize < sizeof(WheelPool) / sizeof(WheelPool[0]));
	//assert (BucketInfo.StockSize + blocks < sizeof(StockCache) / sizeof(StockCache[0]));

	//align by MEM_WHEEL
	pwheel = (WheelPrime*)((size_t)pwheel + MEM_WHEEL - (uint)((size_t)pwheel) % MEM_WHEEL);

	for (int j = 0; j < blocks; j++) {
		for (uint i = 0; i < MEM_BLOCK; i ++) {
			Stock* pStock = StockCache + i + BucketInfo.StockSize;
			pStock->Wheel = pwheel + WHEEL_SIZE * i;
			pStock->Next = pStock + 1;
		}
		StockCache[BucketInfo.StockSize + MEM_BLOCK - 1].Next = StockHead;
		StockHead = StockCache + BucketInfo.StockSize;

		BucketInfo.StockSize += MEM_BLOCK;
		BucketInfo.CurStock +=  MEM_BLOCK;
		pwheel += WHEEL_SIZE * MEM_BLOCK;
	}
}

#define	PI(x, r) (x / log((double)x) * (1 + r / log((double)x)))
static double calR(double p, double mul)
{
	int rr[30] = {0};
	for (int j = 0; j < 36; j++)
	{
		int p1 = Pattern210[j], p2 = Pattern210[j + 1];
		if (p2 < p1)
			p2 = p1 + 2;

		int diff = p2 - p1;
		for (int j = 1; j <= diff; j++)
			rr[j] += diff - j + 1;
	}

	double sum = 0;
	double pir = PI(p, 1.1);
	for (int i = 1; rr[i + 1]; i++)
	{
		double r = PI(p / i, 1.1) / pir - PI(p / (i + 1), 1.1) / pir;
		printf("%d %3d, 1/%d * %d %2lf %lf\n", i, rr[i],  i*i + i, rr[i + 1], r, r * rr[i + 1]);
		sum += r * rr[i + 1];
	}
	printf("%u ration = %.3lf%%, size ration = %.3lf%%\n", (uint)p, 100*sum/210, Config.SieveSize * 100.0 / p);
	return sum / 210;
}

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
			MediumWheel[j +0].Wp = (p / WHEEL << SI_BIT) + pi;
			MediumWheel[j ++].Si = (sieve_index / WHEEL30 << SI_BIT) + wf.WheelIndex;
		}
	}

	MediumWheel[j + 0].Wp = UINT_MAX;
	MediumWheel[j + 1].Wp = UINT_MAX;
}

static void pushBucket(const uint sieve_index, const uint wp, const uchar wheel_index)
{
	const uint si = (sieve_index & BucketInfo.SieveSize) << SI_BIT | wheel_index;
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
		pbucket->Wheel = (wheel = phead->Wheel) + 1;
	}

	wheel->Wp = wp;
	wheel->Si = si;
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
				pbucket->Wheel = (wheel = phead->Wheel) + 1;
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
	crossOffWheelFactor(ppbeg, pend, p);
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
		sieveSmall3(bitarray, pend, p, sieve_index, multiples);
	}
}

static void eratSieveL0(uchar bitarray[], const uint64 start, const int segsize, uint maxp)
{
	const uchar* pend = bitarray + segsize / WHEEL30;
	if (start + segsize < ((uint64)maxp) * maxp) {
		maxp = isqrt(start + segsize) + 1;
	}

	for (uint p = Prime[FIRST_INDEX], j = FIRST_INDEX; p < maxp; p = Prime[++j]) {
		uint sieve_index = p - (uint)(start % p);
		if (start <= p) {
			sieve_index = p * p - (uint)start;
		}

		const WheelFirst wf = WheelFirst30[sieve_index % WHEEL30][WheelInit30[p % WHEEL30].PrimeIndex];
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

static void eratSieveSmall(uchar bitarray[], const uint64 start, const uint sieve_size, const int flag)
{
//	#pragma omp parallel for if (sieve_size > 2 * Threshold.L1Size)
//	preSieve(bitarray, start, sieve_size);
	for (uint sieve_index = 0; sieve_index < sieve_size; sieve_index += Threshold.L1Size) {
		int l1size = Threshold.L1Size;
		if (l1size + sieve_index > sieve_size)
			l1size = sieve_size - sieve_index;
		preSieve(bitarray + sieve_index / WHEEL30, start + sieve_index, l1size);
		if (flag == 1)
			eratSieveL1(bitarray + sieve_index / WHEEL30, start + sieve_index, l1size, Threshold.L1Maxp);
		else
			eratSieveL0(bitarray + sieve_index / WHEEL30, start + sieve_index, l1size, Threshold.L1Maxp);
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
inline static int sieveMedium(uchar bitarray[], const uint sieve_byte, WheelPrime* pwheel)
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
	while ((int)sieve_index1 < 0) {
		SAFE_SET(1);
	}
	while ((int)sieve_index2 < 0) {
		SAFE_SET(2);
	}
	wheel0 = sieve_index0 << SI_BIT | we0.WheelIndex;
	wheel1 = sieve_index1 << SI_BIT | we1.WheelIndex;
	wheel2 = sieve_index2 << SI_BIT | we2.WheelIndex;

	return 3;
}

//sieve prime multiples in [start, start + sieve_size) with medium algorithm
static void eratSieveMedium(uchar bitarray[], const uint64 start, const uint sieve_size, const uint whi, uint maxp)
{
	if (start + sieve_size < (uint64)maxp * maxp) {
		maxp = isqrt(start + sieve_size) + 1;
	}

	const uint sieve_byte = sieve_size / WHEEL30 + WheelInit30[sieve_size % WHEEL30].WheelIndex;
	WheelPrime* pwheel = MediumWheel + whi;

	maxp = (maxp / WHEEL << SI_BIT) + WHEEL_INIT[maxp % WHEEL].PrimeIndex;
	bitarray += sieve_byte;

	while (pwheel[2].Wp < maxp) {
		pwheel += sieveMedium(bitarray, sieve_byte, pwheel);
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
	if (sieve_index1 < sieve_size) {
		wheel1 = wh1 + wheel1->WheelIndex;
		bitarray[sieve_index1] |= wheel1->UnsetBit;
		sieve_index1 += wheel1->Correct + wheel1->NextMultiple * (wp1 >> WP_BIT);
	}

	bitarray[sieve_index2 >>= SI_BIT] |= wheel2->UnsetBit;
	sieve_index2 += wheel2->Correct + wheel2->NextMultiple * (wp2 >> WP_BIT);
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
		WheelPrime* cur_wheel = phead->Wheel;

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
	if (COUNT_PRIME == oper) {
		primes = countBit0sArray((uint64*)bitarray, bytes * sizeof(uint64));
	} else if (PCALL_BACK == oper) {
		primes = doCall((ushort*)bitarray, start, bytes, cmd->Primes, (sieve_call)cmd->Data);
	} else if (SAVE_BYTE == oper) {
		primes = savePrimeByte((stype*)bitarray, bytes, cmd->Data + cmd->Primes);
	}

	if (cmd)
		cmd->Primes += primes;
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
		eratSieveSmall(buffer, start + sieve_index, l2size, 0);
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

	eratSieveSmall(bitarray, start, sieve_size, 1);
	const uchar* pend = bitarray + sieve_size / WHEEL30;
	for (uint j = Threshold.L1Index, p = Threshold.L1Maxp; p <= sqrtp; p = Prime[++j]) {
		const uint sieve_index = p - start % p;
		const WheelFirst& wf = WheelFirst30[sieve_index % WHEEL30][WheelInit30[p % WHEEL30].PrimeIndex];
		const ushort multiples = WHEEL_SKIP >> (wf.NextMultiple * 2);
		sieveSmall1(bitarray, pend, p, sieve_index + wf.Correct * p, multiples);
	}

	*(uint64*)(bitarray + bytes) = ~0;
	return bytes;
}

static void setThreshold()
{
	for (uint p = Prime[1], j = 1; ; p = Prime[++j]) {
		if (p >= Threshold.L1Maxp) {
			Threshold.L1Index = j;
			Threshold.L1Maxp = p;
			return;
		}
	}
	assert(Threshold.L1Maxp <= 65536 && Threshold.L1Maxp < Threshold.L2Maxp);
}

void setCpuCache(uint cdata)
{
	cdata = 1 << ilog(cdata, 2);

	if (cdata >= 16 && cdata < L2_DCACHE_SIZE) { //L1
		Threshold.L1Size = cdata * (WHEEL30 << 10);
		Threshold.L1Maxp = Threshold.L1Size / (WHEEL30 * Threshold.L1Segs);
		setThreshold();
	} else if (cdata >= L2_DCACHE_SIZE && cdata <= MAX_SEGMENT) { //L2
		Threshold.L2Size = cdata * (WHEEL30 << 10);
		Threshold.L2Maxp = Threshold.L2Size / (WHEEL30 * Threshold.L2Segs);
	}
}

void setLevelSegs(uint cdata)
{
	const int level = cdata / 10, segs = cdata % 10;
	if (segs > 1) {
		if (level == 1) {
			Threshold.L1Segs = segs;
			Threshold.L1Maxp = Threshold.L1Size / (WHEEL30 * segs);
			setThreshold();
		} else if(level == 2) {
			Threshold.L2Segs = segs;
			Threshold.L2Maxp = Threshold.L2Size / (WHEEL30 * segs);
		} else if (segs <= 6) {
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
	const uint smallPrime[] = {3, 5, 7, 0, 0, 0};
	uint primes = 0;
	for (int i = 0; i < 5; i++) {
		uint p = smallPrime[i];
		if (start <= p && (p == smallPrime[i + 1] - PRIME_GAP || p == smallPrime[i + 2] - PRIME_GAP) && smallPrime[i + 1] <= end) {
			primes ++;
			if (cmd && cmd->Oper == PCALL_BACK) {
				(*(sieve_call)cmd->Data)(primes, p);
				cmd->Primes += 1;
			}
		}
	}

	return primes;
}

static uint64 pi2(uint64 start, uint64 end, Cmd* cmd)
{
	const int64 ts = getTime();
	uint64 primes = checkSmall(start, end, cmd);

	uint start_align = (uint)(start % WHEEL210);
	start -= start_align;

	if (++end == 0) end --; //watch overflow if end = 2^64-1
	const int64 range = (int64)(end - start);
	uchar* bitarray = (uchar*) malloc((Config.SieveSize + Threshold.L1Size) / WHEEL30) + 8;
	*(uint64*)(bitarray - 8) = ~0;
	uchar last_byte = ~0;

	for (uint si = 0, sieve_size = Config.SieveSize; start < end; start += sieve_size) {
		uint bytes = sieve_size / WHEEL30;
		if (sieve_size > end - start) {
			sieve_size = (uint)(end - start);
			bytes = sieve_size / WHEEL30 + (7 + WheelInit30[sieve_size % WHEEL30].WheelIndex) / 8;
		}
		
		segmentedSieve(bitarray, start, sieve_size);
		if (start_align > 0) {
			memset(bitarray, ~0, start_align / WHEEL30);
			bitarray[start_align / WHEEL30] |= (1 << WheelInit30[start_align % WHEEL30].WheelIndex) - 1;
			start_align = 0;
		}
  
		//pack adjacent segment byte
		bitarray[-1] = last_byte | 0x7f; last_byte = *(bitarray + bytes - 1);
		*(uint64*)(bitarray + bytes) = ~0;

		primes += segmentProcessed(bitarray - 8, start, bytes + 8, cmd);
		if ((si ++ & Config.Progress) == 15) {
			double ratio = 100 - 100.0 * (int64)(end - start - sieve_size) / range;
			double timeuse = (getTime() - ts) / (10 * ratio);
			printf(">> %.2f%%, sieve time ~= %.2f sec, primes ~= %llu\r", ratio, timeuse, (int64)((100 * (int64)primes) / ratio));
		}
	}

	free(bitarray - 8);
	return primes;
}

static void convertSci(uint64 n, char buff[20])
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
	} else if (begbuff[0] && endbuff[0]) {
		sprintf(buff, "%s,%s", begbuff, endbuff);
	} else if (rangebuff[0]) {
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
	printf("\rpi2(%s) = %llu", buff, primes);
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
	initPrime(SIEVE_SIZE);
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

	memset((char*)&BucketInfo, 0, sizeof(BucketInfo));
	//init bucket sieve
	uint64& buckets = Threshold.BucketStart;
	if (sqrtp > medium) {
		setBucketInfo(sieve_size, sqrtp, end - buckets);
//		printf("\t[0] all = %d, pool = %d, max = %d, cur = %d\n",
//				BucketInfo.StockSize, BucketInfo.PoolSize, BucketInfo.LoopSize, BucketInfo.MaxBucket);
		setWheelPrime(medium, sqrtp, buckets, end - buckets);
//		printf("\t[1] all = %d, pool = %d, remain = %d, init = %d ms\n",
//				BucketInfo.StockSize, BucketInfo.PoolSize, BucketInfo.CurStock, (int)(getTime() - ts));
		allocWheelBlock(1);
	} else {
		buckets = end + 1;
	}

	const int64 ti = getTime();
	const uint64 primes = pi2(start, end, cmd);

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

static void initBitTable( )
{
	int i = 0, n = 0;
	for (i = 2; Pattern30[i]; i++) {
		if (Pattern30[i] - Pattern30[i - 1] == PRIME_GAP) {
			TwBits[n++] = i - 2;
		}
	}

	for (i = 0; i < sizeof(WordNumBit0) / sizeof (WordNumBit0[0]); i ++) {
		for (int j = 0; j < sizeof(TwBits) / sizeof(TwBits[0]); j++) {
			WordNumBit0[i] += (i >> TwBits[j]) % 4 == 0;
		}
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

			WheelElement& wdata = Wheel30[pi][wi];
			wdata.NextMultiple = multiples * (WHEEL30 / WHEEL30);
			wdata.WheelIndex = WheelInit30[next % WHEEL30].WheelIndex;
			wdata.Correct = next / WHEEL30 - Pattern30[wi] / WHEEL30;
			wdata.UnsetBit = WheelInit30[Pattern30[wi]].UnsetBit;
		}
	}
}

static void initWheel210()
{
	int wi = 0, i = 0;
	const int psize = 48, wsize = 30;
	int wpattern[WHEEL210] = {0};

	for (i = 0; i < WHEEL210; i ++) {
		const uchar mask = WheelInit30[i % WHEEL30].UnsetBit;
		WheelInit210[i].UnsetBit = mask;
		WheelInit210[i].WheelIndex = -1;
		WheelInit210[i].PrimeIndex = wi;
		if (mask && i % (WHEEL210 / WHEEL30)) {
			Pattern210[wi ++] = i;
		}
	}

	wi = 0, i = 1;
	if (PRIME_GAP == 2) {
		WheelInit210[i].WheelIndex = wi;
		wpattern[wi ++] = i;
	}

	for (; i < WHEEL210; i += 2) {
		const uchar mask = WheelInit210[i % WHEEL30].UnsetBit;
		if (mask && i % (WHEEL210 / WHEEL30) &&
			(i + PRIME_GAP) % (WHEEL210 / WHEEL30) &&
			WheelInit210[(i + PRIME_GAP) % WHEEL30].UnsetBit) {
			WheelInit210[i].WheelIndex = wi;

			wpattern[wi ++] = i;
			if (i + PRIME_GAP < WHEEL210) {
				WheelInit210[i + PRIME_GAP].WheelIndex = wi;
				wpattern[wi ++] = (i + PRIME_GAP) % WHEEL210;
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

			WheelElement& wdata = Wheel210[pi][wi];
			wdata.Correct = next / WHEEL30 - wpattern[wi] / WHEEL30;
			wdata.UnsetBit = WheelInit210[wpattern[wi]].UnsetBit;
			wdata.WheelIndex = WheelInit210[next % WHEEL210].WheelIndex;
			wdata.NextMultiple = multiples * (WHEEL210 / WHEEL30);
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

static void fixRangeTest(uint64 lowerBound, const int64 range, uint64 Ret)
{
	const uint llog10 = ilog(lowerBound, 10), rlog10 = ilog(range, 10);
	printf("Sieving the primes within (10^%u, 10^%u+10^%u) randomly\n", llog10, llog10, rlog10);
	uint64 primes = 0, upperBound = lowerBound + range;

	while (lowerBound < upperBound) {
		uint64 rd = rand() * rand();
		uint64 end = lowerBound + (rd * rand() * rd) % ipow(2, 32) + ipow(10, 4);
		end = end - end % WHEEL210 + 6;
		if (end > upperBound)
			end = upperBound;
		setSieveSize(rand() % MAX_SEGMENT + 128);
		setLevelSegs(rand() % 5 + 12), setLevelSegs(rand() % 5 + 22);
		setLevelSegs(rand() % 3 + 34);
		primes += doSieve(lowerBound, end - 1);
		printf("Remaining chunk: %.2f%%\r", (int64)(upperBound - lowerBound) * 100.0 / range);
		lowerBound = end;
	}

	if (lowerBound + 1 == 0)
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
		primes = doSieve(0, ipow(10, i));
		printf("pi2(10^%02d) = %llu\n", i, primes);
	}

	srand(time(0));
	for (int j = 12; primeCounts[j]; j ++) {
		uint64 start = ipow(10, j), end = start + ipow(2, 32);
		setLevelSegs(rand() % 5 + 12), setLevelSegs(rand() % 5 + 22), setLevelSegs(rand() % 4 + 34);
		primes = doSieve(start, end);
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
		{ipow(10, 15), pow11, 110670758},
		{ipow(10, 17), pow11, 86176910},
		{ipow(10, 19), pow11, 68985092},
		{ipow(10, 14), pow12, 1270127074},
		{ipow(10, 16), pow12, 972773783},
		{ipow(10, 18), pow12, 768599834},
		{-1 - pow12,   pow12, 670910567},
	};

	for (int k = 0; k < sizeof(RangeData) / sizeof(RangeData[0]); k ++) {
		fixRangeTest(RangeData[k][0], RangeData[k][1], RangeData[k][2]);
	}

	Config.Flag |= PRINT_RET;
}

#if X86_64 || _M_IX86 || __i386__
static void cpuid(int cpuinfo[4], int id)
{
#if _MSC_VER > 1300
	//__cpuid(cpuinfo, id);
#elif _MSC_VER >= 1200
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
#elif __GNUC__
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

// http://msdn.microsoft.com/en-us/library/hskdteyh%28v=vs.100%29.aspx
static int getCpuInfo()
{
	char cpuName[257] = {-1};
	int (*pTmp)[4] = (int(*)[4])cpuName;
	cpuid(*pTmp ++, 0x80000002);
	cpuid(*pTmp ++, 0x80000003);
	cpuid(*pTmp ++, 0x80000004);

	for (int i = 0; cpuName[i]; i ++) {
		if (cpuName[i] != ' ' || cpuName[i + 1] != ' ')
			putchar(cpuName[i]);
	}

	int cpuinfo[4];
	cpuid(cpuinfo, 0x80000006);
	putchar('\n');

	if (cpuName[0] == 'A') { //amd cpu
		Threshold.L1Size = 64 * (WHEEL30 << 10);
	} else {
		Threshold.L1Size = 32 * (WHEEL30 << 10);
	}

	Threshold.L1Maxp = Threshold.L1Size / (WHEEL30 * Threshold.L1Segs);

	return cpuinfo[2] >> 16;
}
#endif

static void printInfo( )
{
	const char* sepator =
		"---------------------------------------------------------------------------------------";
	puts(sepator);
	puts("Fast implementation of the segmented sieve of Eratosthenes (2^64 - 1)\n"
	"Copyright (C) by 2010-2017 Huang YuanBing bailuzhou@163.com\n"
	"g++ -DSIEVE_SIZE=1024 -DL2_DCACHE_SIZE=256 -march=native -funroll-loops -O3 -s -pipe\n");
//	"Code: https://github.com/ktprime/ktprime/blob/master/TwinPrime.cpp\n"

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
			MEM_BLOCK * WHEEL_SIZE >> 17, WHEEL_SIZE >> 7, SIEVE_SIZE, WHEEL);
	info += sprintf(info, "[ARGS ] : L1Size = %d, L2Size = %d, SieveSize = %d, Loop = %d\n",
			Threshold.L1Size / WHEEL30 >> 10, Threshold.L2Size / WHEEL30 >> 10, Config.SieveSize / WHEEL30 >> 10, BucketInfo.LoopSize);
	info += sprintf(info, "[ARGS ] : L1Seg/L1Maxp/L2Seg/L2Maxp/Medium/Large = (%d,%d,%d,%d,%d,%u)",
		Threshold.L1Segs, Threshold.L1Maxp, Threshold.L2Segs, Threshold.L2Maxp, Threshold.Medium, isqrt(Threshold.BucketStart + 1));
	puts(buff);
	puts(sepator);
}

static void doCompile(const char* flag)
{
	char programming[257];
	strcpy(programming, __FILE__);
	char* pdot = strchr(programming, '.');
	if (pdot) {
		strcpy(pdot, "_.exe");
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
	int cmdi = -1;
	int cdata = 0;

	initPrime(SIEVE_SIZE);
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
				setCpuCache(cdata);
				break;
			case 'L':
				setLevelSegs(cdata);
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

bool executeCmd(const char* cmd)
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

		char cmdc = toupper(params[cmdi][0]);
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
		} else if (cmdc == 'P') {
			Cmd cmd = {PCALL_BACK, (uchar*)printPrime, 0};
			doSieve(start, end, &cmd);
		} else if (cmdc == 'Z') {
			doCompile(params[cmdi + 1]);
		} else if (cmdi >= 0 && end > 0) {
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
#if X86_64 || _M_IX86 || __i386__
		getCpuInfo();
#endif

		eratoSimple();
		initBitTable( );
		initWheel30( );
		initWheel210( );
		initPreSieved( );
		setThreshold();
		setSieveSize(sieve_size);
	}
}

#ifndef PRIME_LIB
int main(int argc, char* argv[])
{
	initPrime(SIEVE_SIZE);

	if (argc > 1)
		executeCmd(argv[1]);

#if MTEST
	srand(time(0));
	for (int j = 1; j < 100; j++) {
		for (int i = 2; i <= 6; i++) {
			setSieveSize(256);
//			setLevelSegs(i + 32); //setLevelSegs(rand() % 5 + 12), setLevelSegs(rand() % 5 + 22);
			const uint64 medium = Config.SieveSize / i;
			uint64 start = medium * medium;
			start -= start % WHEEL210;


			const uint64 range = ipow(10, 7) + (uint)rand() * rand();
			const uint64 r1 = doSieve(start - range, start + range, NULL);

			setSieveSize(2048);
//			setLevelSegs(36); setLevelSegs(rand() % 5 + 12), setLevelSegs(rand() % 5 + 22);
			const uint64 r2 = doSieve(start - range, start + range, NULL);
			if (r1 != r2)
				printf("%llu %llu  r1 = %ld, r2 = %ld\n", start - range, start + range, r1, r2);
		}
	}
#endif

#if GCOV_TEST
	executeCmd("m5 s1024 l14; 1e12 1e12+1e9"); putchar('\n');
	executeCmd("1e18 1e9"); putchar('\n');
	executeCmd("c32 l22 1e16 10^9"); putchar('\n');
	executeCmd("l34 G 1e12 e8"); putchar('\n');
	executeCmd("da s256 1e5; 1e8"); putchar('\n');
	executeCmd("s4000 c2000 p 0 100"); putchar('\n');

	executeCmd("da C32 s1024 H I 1e18 1e8"); putchar('\n');
	executeCmd("df p 1e12+100 2e2*2; c1222222; g 1e16 1e8; s256 l32 l24 l14 1e8"); putchar('\n');
	executeCmd("df s1024");
	executeCmd("0-1e4 0-1");
#endif

#ifndef B_R
	if (argc == 1)
	executeCmd("1e12 1e9; 1e18 1e9; 1e16 1e9; 0-1000000001 0-1;e16 e10;");
#else
	executeCmd("e16 e10");
#endif

	while (true) {
		char ccmd[257];
		printf("\n>> ");
#if _MSC_VER > 1600
		if (!gets_s(ccmd) || !executeCmd(ccmd))
#else
		if (!gets(ccmd) || !executeCmd(ccmd))
#endif
			break;
	}

	return 0;
}
#endif
