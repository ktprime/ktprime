//fast segmented sieving of twin prime
/***
MINGW: gcc 4.7.3
CXXFLAG:g++ -march=native [-DW210,-DSAFE=1] -funroll-loops -O3 -s -pipe;
Windows 7 x64,              I3 350M 2.26G / i5 3470 3.2G
pi2(1e11, 1e11+1e10) = 20498568      04.30/2.30
pi2(1e12, 1e12+1e10) = 17278660      04.96/2.63
pi2(1e13, 1e13+1e10) = 14735239      05.95/3.01
pi2(1e14, 1e14+1e10) = 12706059      07.37/3.67
pi2(1e15, 1e15+1e10) = 11069670      08.77/4.00
pi2(1e16, 1e16+1e10) = 9728024       10.21/5.09
pi2(1e17, 1e17+1e10) = 8614943       12.01/6.02
pi2(1e18, 1e18+1e10) = 7687050       14.92/7.56
pi2(1e19, 1e19+1e10) = 6895846       21.27/11.12

pi2(1e18, 1e18+1e6) = 794            0.72/0.49
pi2(1e18, 1e18+1e8) = 77306          1.54/0.80
pi2(1e18, 1e18+1e9) = 769103         3.80/1.72
pi2(1e18, 1e18+1e12)= 24127637783    1640/860
pi2(1e16, 1e16+1e12)= 27143405794    1310/690
pi2(1e14, 1e14+1e12)= 31016203073     950/500

doc:
	http://www.ieeta.pt/~tos/software/prime_sieve.html
	http://code.google.com/p/primesieve/wiki/Links
	http://primesieve.org/
	http://www.compileonline.com/
	http://liveworkspace.org/code/
	http://coliru.stacked-crooked.com/
***/

#include <ctype.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>

//const
# define WHEEL30          30
# define WHEEL210         210
# define PRIME_PRODUCT    (210 * 11 * 13 * 17 * 19)
# define FIRST_INDEX      PRIME_PRODUCT / 9699690 + 8
# define WHEEL_SKIP       0x799b799b
//2 for twin, 4 for Cousin
# define PRIME_GAP        2

//performance marco
# define L1_DCACHE_SIZE   64
# define SIEVE_BIT        8 //8 - 10
# define WHEEL_BIT        6 //6 - 7

# define ERAT_SMALL       6 //4 - 16
# define ERAT_MEDIUM      2 //2 - 6
# define ERAT_BIG         6 //2 - 6

//x86 cpu only
#ifndef POPCNT
# define POPCNT           0
#endif

#ifndef SAFE
# define SAFE             0
#endif

#ifndef L2_DCACHE_SIZE
# define L2_DCACHE_SIZE   256
#endif
#if MAX_SIEVE < 256
# define MAX_SIEVE        2048 //>= L2_DCACHE_SIZE < MAX_LEVEL_SIZE / 2
#endif

#ifdef W30 //some intel cpu
	# define WHEEL        WHEEL30
	# define WHEEL_MAP    Wheel30
	# define WHEEL_INIT   WheelInit30
	# define WHEEL_FIRST  WheelFirst30
#else
	# define WHEEL        WHEEL210
	# define WHEEL_MAP    Wheel210
	# define WHEEL_INIT   WheelInit210
	# define WHEEL_FIRST  WheelFirst210
#endif

#if _MSC_VER > 1300
	# include <intrin.h>
#endif

#if __x86_64__ || _M_AMD64 || __amd64__
	# define X86_64       1
#endif

#if X86_64 && _MSC_VER
	# define ASM_X86      0
	# define BIT_SCANF    1
#elif X86_64 || _M_IX86 || __i386__
	# define ASM_X86      1
	# define BIT_SCANF    1
#else
	# define ASM_X86      0
	# define BIT_SCANF    0
#endif

#ifdef _WIN32
	typedef unsigned __int64 uint64;
	typedef __int64 int64;
	#define CONSOLE "CON"
	#include <windows.h>
#else
	typedef unsigned long long uint64;
	typedef long long int64;
	#define CONSOLE "/dev/tty"
	#include <sys/time.h>
	#include <unistd.h>
#endif

typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;

#define MIN(a, b)         (a < b ? a : b)

#if BIT_SCANF == 0
	#define PRIME_OFFSET(mask)   Lsb[mask]
	typedef ushort stype;
#else
	#define PRIME_OFFSET(mask)   Pattern30[bitScanForward(mask)]
#if X86_64
	typedef uint64 stype;
#else
	typedef uint stype;
#endif
#endif

static const char* const Help = "\
	[B: Benchmark [0]]\n\
	[D: D[T,R,A,F] time, result, run, save]\n\
	[M: Progress of calculating (0 - 20)]\n\
	[C: Cpu L1/L2 data cache size (L1:16-128, L2:256-1024)]\n\
	[S: Set sieve segment size (16 - 1024)]\n\
	[L: Set sieve segs L1(2-6), L2(2-6) L3(2-6)]\n\n\
	[I: Info of programming]\n\
	[Z: Compile self and run]\n\
	[P: Print prime in [start, end]]\n\n\
Example:\n\
	1e16 1e16+10^10 s512\n\
	p 1e12+100 100";

enum EFLAG
{
	SLOW_TEST  = 1 << ('A' - 'A'),
	PRINT_RET  = 1 << ('R' - 'A'),
	PRINT_TIME = 1 << ('T' - 'A'),
	SAVE_DATA  = 1 << ('F' - 'A'),
};

enum ECMD
{
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

enum EBUCKET
{
	UINT_PIMAX = 203280221, // = pi(2^32)
	MAX_BUCKET = 5464 * 3, //= 28*2^32 / 2^18 * 30 + 4
	WHEEL_SIZE = 1 << 12, //11: 16k, [10 - 13]
	MEM_WHEEL  = WHEEL_SIZE * 8, // = sizeof(WheelPrime),
	MEM_BLOCK  = (1 << 20) / WHEEL_SIZE, //20: 8 MB
	MAX_CACHE  = (MAX_SIEVE + L1_DCACHE_SIZE) << 10, //1080 k
	MAX_STOCK  = UINT_PIMAX / WHEEL_SIZE + MAX_BUCKET, //104722
};

static struct
{
	uint L1Size; //cpu L1/L2 cache size
	uint L2Size;

	uint L1Maxp; //limitEratSmall
	uint L1Index;
	uint L1Segs;

	uint L2Maxp; //limitEratMedium
	uint L2Index;
	uint L2Segs;

	uint Medium; //limitEratBig
	uint Msegs;
}
Threshold =
{
	L1_DCACHE_SIZE * (WHEEL30 << 10), L2_DCACHE_SIZE * (WHEEL30 << 10),
	(L1_DCACHE_SIZE << 10) / ERAT_SMALL, 0, ERAT_SMALL,
	(L2_DCACHE_SIZE << 10) / ERAT_MEDIUM, 0, ERAT_MEDIUM,
	MAX_SIEVE * (WHEEL30 << 10) / ERAT_BIG, ERAT_BIG,
};

struct _Config
{
	uint Flag;
	//print calculating progress
	uint Progress;
	//sieve size
	uint SieveSize;
	//min offset bucket
	uint64 MinBucket;
};

_Config Config =
{
	PRINT_RET | PRINT_TIME,
	(1 << 6) - 1,
	MAX_SIEVE * (WHEEL30 << 10),
	0,
};

struct WheelPrime
{
	//[0 - 6]: p % wheel, [7 - 31]: p / sieve_size
	uint Wp;
	//[0 - 7]: index % sieve_size, [8 - 31]: index / 30
	uint SieveIndex;
};

struct Stock
{
	//read only
	WheelPrime* Wheel;
	Stock* Next;
};

struct _Bucket
{
	//write only
	WheelPrime* Wheel;
	//Head ->Stock1 ->Stock2 ->.... ->Stockn
	Stock* Shead;
};

struct _BucketInfo
{
	uint SieveSize;
	uint Log2Size;
	uint MaxBucket;
	uint CurBucket;

	uint CurStock;
	uint StockSize;
	uint PtrSize;
};

//thread ...
static _BucketInfo BucketInfo;
static _Bucket Bucket [MAX_BUCKET];
static Stock* StockPool [MAX_STOCK];
static WheelPrime* WheelPtr [(1 << 17) / MEM_BLOCK]; //2G vm

#if CHECK
	static uchar Prime[UINT_PIMAX + 1000];
#else
	static uchar Prime[MAX_CACHE];
#endif

static WheelPrime* MediumWheel;

//presieved small prime number <= 19.
//the crossing out bit module WHEEL30, the first
//16 bit of PreSieved map to
//----------------------------------------
//|01/1|07/1|11/1|13/1|17/1|19/1|23/1|29/1| = 0x1111 1111 = PreSieved[0]
//----------------------------------------
//|31/1|37/1|41/1|43/1|47/1|49/0|53/1|59/1| = 0x1101 1111 = PreSieved[1]
//----------------------------------------
static uchar PreSieved[PRIME_PRODUCT / WHEEL30];

#if BIT_SCANF == 0
static uchar Lsb[1 << 16];
#endif

//number of bits 1 binary representation table in Range[0-2^16)
static uchar WordNumBit0[1 << 16];

struct WheelElement
{
#ifndef _LITTLE_ENDIAN
	uchar WheelIndex;
	uchar UnsetBit;
	uchar Correct;
	uchar NextMultiple;
#else
	uchar NextMultiple;
	uchar Correct;
	uchar UnsetBit;
	uchar WheelIndex;
#endif
};

struct WheelInit
{
	char WheelIndex;
	uchar UnsetBit;
	uchar PrimeIndex;
};

typedef WheelElement WheelFirst;
static WheelInit WheelInit30[WHEEL30];
static WheelElement Wheel30[8][8];
static WheelFirst WheelFirst30[WHEEL30][8];

static WheelInit WheelInit210[WHEEL210];
static WheelElement Wheel210[48][64];
static WheelFirst WheelFirst210[WHEEL210][64];

static uchar TwBits[6] =
{
	1, 3, 6, 9, 11, 14
//	0, 2, 4, 8, 10, 12
};

static const uchar SmallPrime[] =
{
	2, 3, 5
};

static uchar Pattern30[64] =
{
	1, 7, 11, 13, 17, 19, 23, 29,
	31, 37, 41, 43, 47, 49, 53, 59, 61
};

//api
#ifndef PRIME_LIB
struct Cmd
{
	uint Oper;
	uchar* Data;
	uint64 Primes;
};

typedef void (*sieve_call)(uint64, uint64);
void initPrime(int sieve_size);
bool executeCmd(const char* cmd);
uint64 doSieve(const uint64 start, const uint64 end, Cmd* cmd);
int setSieveSize(uint sieve_size);
void setLevelSegs(uint level);
void setCpuSize(uint cpusize);
#else
#include "PrimeNumber.h"
#endif

static int64 getTime( )
{
#ifdef _WIN32
	LARGE_INTEGER freq, count;
	QueryPerformanceFrequency(&freq);
	QueryPerformanceCounter(&count);
	return 1000ul * count.QuadPart / freq.QuadPart;
#else
	struct timeval tmVal;
	gettimeofday(&tmVal, NULL);
	return tmVal.tv_sec * 1000ul + tmVal.tv_usec / 1000;
#endif
}

//base > 1
static int ilog(uint64 n, const uint base)
{
	int powbase = 0;
	while (n / base) {
		powbase ++;
		n /= base;
	}

	return powbase;
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

//invalid is [number][e][+-*][number]
//ex: 1000, e9, 10^7, 2^32*2, 2^30-1E2, 2e9+2^20
uint64 atoint64(const char* str)
{
	uint64 ret = 0;

	while (isspace(*str)) {
		str ++;
	}

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
	if (ps = strchr(str, '+')) {
		ret += atoint64(ps + 1);
	} else if (ps = strchr(str, '-')) {
		ret -= atoint64(ps + 1);
	} else if (ps = strchr(str, '*')) {
		ret *= atoi(ps + 1);
	}

	return ret;
}

//30% fast on old cpu pentium 4, fam3
inline static uchar* crossOffWheelFactor2(uchar* p, const uchar* pend, const uint step)
{
	const uint o = step / WHEEL30;
	switch (WheelInit30[step % WHEEL30].PrimeIndex) {
	case 0 :
		while (p < pend) {
			p[o *  0 +  0] |= BIT0, p[o *  6 +  0] |= BIT1;
			p[o * 10 +  0] |= BIT2, p[o * 12 +  0] |= BIT3;
			p[o * 16 +  0] |= BIT4, p[o * 18 +  0] |= BIT5;
			p[o * 22 +  0] |= BIT6, p[o * 28 +  0] |= BIT7;
			p += step;
		}
		break;
	case 1 :
		while (p < pend) {
			p[o *  0 +  0] |= BIT0, p[o *  4 +  0] |= BIT7;
			p[o *  6 +  1] |= BIT3, p[o * 10 +  2] |= BIT2;
			p[o * 16 +  3] |= BIT6, p[o * 18 +  4] |= BIT1;
			p[o * 24 +  5] |= BIT5, p[o * 28 +  6] |= BIT4;
			p += step;
		}
		break;
	case 2 :
		while (p < pend) {
			p[o *  0 +  0] |= BIT0, p[o *  2 +  0] |= BIT6;
			p[o *  6 +  2] |= BIT1, p[o *  8 +  2] |= BIT7;
			p[o * 12 +  4] |= BIT3, p[o * 18 +  6] |= BIT5;
			p[o * 20 +  7] |= BIT2, p[o * 26 +  9] |= BIT4;
			p += step;
		}
		break;
	case 3 :
		while (p < pend) {
			p[o *  0 +  0] |= BIT0, p[o *  4 +  1] |= BIT6;
			p[o *  6 +  2] |= BIT5, p[o * 10 +  4] |= BIT2;
			p[o * 12 +  5] |= BIT1, p[o * 16 +  6] |= BIT7;
			p[o * 22 +  9] |= BIT4, p[o * 24 + 10] |= BIT3;
			p += step;
		}
		break;
	case 4:
		while (p < pend) {
			p[o *  0 +  0] |= BIT0, p[o *  6 +  3] |= BIT3;
			p[o *  8 +  4] |= BIT4, p[o * 14 +  7] |= BIT7;
			p[o * 18 + 10] |= BIT1, p[o * 20 + 11] |= BIT2;
			p[o * 24 + 13] |= BIT5, p[o * 26 + 14] |= BIT6;
			p += step;
		}
		break;
	case 5 :
		while (p < pend) {
			p[o *  0 +  0] |= BIT0, p[o *  4 +  2] |= BIT4;
			p[o * 10 +  6] |= BIT2, p[o * 12 +  7] |= BIT5;
			p[o * 18 + 11] |= BIT3, p[o * 22 + 13] |= BIT7;
			p[o * 24 + 15] |= BIT1, p[o * 28 + 17] |= BIT6;
			p += step;
		}
		break;
	case 6 :
		while (p < pend) {
			p[o *  0 +  0] |= BIT0, p[o *  2 +  1] |= BIT4;
			p[o *  6 +  4] |= BIT5, p[o * 12 +  9] |= BIT1;
			p[o * 14 + 10] |= BIT6, p[o * 20 + 15] |= BIT2;
			p[o * 24 + 18] |= BIT3, p[o * 26 + 19] |= BIT7;
			p += step;
		}
		break;
	case 7 :
		while (p < pend) {
			p[o *  0 +  0] |= BIT0, p[o *  2 +  1] |= BIT7;
			p[o *  8 +  7] |= BIT6, p[o * 12 + 11] |= BIT5;
			p[o * 14 + 13] |= BIT4, p[o * 18 + 17] |= BIT3;
			p[o * 20 + 19] |= BIT2, p[o * 24 + 23] |= BIT1;
			p += step;
		}
		break;
	}

	return p;
}

//fast on new cpu
inline static void crossOffWheelFactor(uchar* ppbeg[], const uchar* pend, const uint p)
{
	#define CMPEQ_OR(n) if (ps##n <= pend) *ps##n |= BIT##n

	uchar* ps0 = ppbeg[0], *ps1 = ppbeg[1];
	uchar* ps2 = ppbeg[2], *ps3 = ppbeg[3];
	while (ps3 <= pend) {
		*ps0 |= BIT0, ps0 += p;
//		*ps1 |= BIT1, ps1 += p;
		*ps2 |= BIT2, ps2 += p;
		*ps3 |= BIT3, ps3 += p;
	}
	CMPEQ_OR(0); CMPEQ_OR(1); CMPEQ_OR(2);
#if 1
	#define CMPEQ_OR2(n, i) if (ps##n <= pend) *ps##n |= BIT##i
	ps0 = ppbeg[4], ps1 = ppbeg[5];
	ps2 = ppbeg[6], ps3 = ppbeg[7];
	while (ps3 <= pend) {
		*ps0 |= BIT4, ps0 += p;
		*ps1 |= BIT5, ps1 += p;
//		*ps2 |= BIT6, ps2 += p;
		*ps3 |= BIT7, ps3 += p;
	}
	CMPEQ_OR2(0, 4); CMPEQ_OR2(1, 5); CMPEQ_OR2(2, 6);
#else
	uchar* ps4 = ppbeg[4], *ps5 = ppbeg[5];
	uchar* ps6 = ppbeg[6], *ps7 = ppbeg[7];
	while (ps7 <= pend) {
		*ps4 |= BIT4, ps4 += p;
		*ps5 |= BIT5, ps5 += p;
//		*ps6 |= BIT6, ps6 += p;
		*ps7 |= BIT7, ps7 += p;
	}
	CMPEQ_OR(4); CMPEQ_OR(5); CMPEQ_OR(6);
#endif
}

inline static void crossOff8Factor(uchar* ppbeg[], const uchar* pend, const uint64 mask64, const uint p)
{
	#define CMQEQ_QRM(n, mask) if (ps##n <= pend) *ps##n |= mask
//	#define CMQEQ_QRM(n, mask) *ps##n |= mask

	uchar* ps0 = ppbeg[0], *ps1 = ppbeg[1];
	uchar* ps2 = ppbeg[2], *ps3 = ppbeg[3];
	uint mask = (uint)(mask64 >> 32);
	while (ps3 <= pend) {
		*ps0 |= mask >> 24, ps0 += p;
		*ps1 |= mask >> 16, ps1 += p;
		*ps2 |= mask >> 8,  ps2 += p;
		*ps3 |= mask >> 0,  ps3 += p;
	}
	CMQEQ_QRM(0, mask >> 24); CMQEQ_QRM(1, mask >> 16); CMQEQ_QRM(2, mask >> 8);

	ps0 = ppbeg[4], ps1 = ppbeg[5];
	ps2 = ppbeg[6], ps3 = ppbeg[7];
	mask = (uint)mask64;
	while (ps3 <= pend) {
		*ps0 |= mask >> 24, ps0 += p;
		*ps1 |= mask >> 16, ps1 += p;
		*ps2 |= mask >> 8,  ps2 += p;
		*ps3 |= mask >> 0,  ps3 += p;
	}
	CMQEQ_QRM(0, mask >> 24); CMQEQ_QRM(1, mask >> 16); CMQEQ_QRM(2, mask >> 8);
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
	const uint hig = (uint)(n >> 32), low = (uint)n;
	return WordNumBit0[(ushort)hig] + WordNumBit0[(ushort)low] + WordNumBit0[hig >> 16] + WordNumBit0[low >> 16];
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

//print prime from bit buffer
static void printPrime(uint64 index, uint64 prime)
{
	printf("(%I64u, %I64u)\n", prime, prime + PRIME_GAP);
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

		int curp = bi * WHEEL30 * sizeof(mask) + PRIME_OFFSET(mask);
		mask &= mask - 1;
		*prime++ = curp - lastp + (curp - lastp) / 256;
		primes++;
		lastp = curp;
	}

	//assert(bytes * WHEEL30 - lastp) < 256
	*prime = lastp - bytes * WHEEL30;

	return primes;
}

static void allocWheelBlock(const uint blocks)
{
	static Stock StockCache [MAX_STOCK];
	WheelPrime *pwheel = (WheelPrime*) malloc((blocks + 1) * MEM_WHEEL);
	WheelPtr[BucketInfo.PtrSize ++] = pwheel;
//	assert (pwheel && BucketInfo.PtrSize < sizeof(WheelPtr) / sizeof(WheelPtr[0]));
//	assert (BucketInfo.StockSize + blocks < sizeof(StockCache) / sizeof(StockCache[0]));

	//align by MEM_WHEEL
	pwheel = (WheelPrime*)((size_t)pwheel + MEM_WHEEL - (size_t)pwheel % MEM_WHEEL);
	for (uint i = 0; i < blocks; i ++) {
		Stock* pStock = StockPool[BucketInfo.CurStock + i] = StockCache + i + BucketInfo.StockSize;
		pStock ->Wheel = pwheel + WHEEL_SIZE * i;
		pStock ->Next = NULL;
	}

	BucketInfo.StockSize += blocks;
	BucketInfo.CurStock += blocks;
}

static int initBucketInfo(const uint sieve_size, const uint sqrtp, const uint64 range)
{
//	assert ((range >> 32) < sieve_size);
	BucketInfo.CurBucket = range / sieve_size + 1;
	if (range % sieve_size == 0)
		BucketInfo.CurBucket --;

	//wheel 210 pattern, max pattern difference is 14 //30 ->10, 210 ->26
	BucketInfo.MaxBucket = (uint64)sqrtp * (26 + 2) / sieve_size + 2;
//	assert (BucketInfo.MaxBucket < sizeof(Bucket) / sizeof(Bucket[0]));

	BucketInfo.Log2Size = ilog(sieve_size / WHEEL30, 2);
	BucketInfo.SieveSize = (1 << BucketInfo.Log2Size) - 1;
//	assert(BucketInfo.Log2Size <= (32 - SIEVE_BIT) && SIEVE_BIT >= WHEEL_BIT);

	uint blocks = MIN(BucketInfo.MaxBucket, BucketInfo.CurBucket);
	memset(Bucket, 0, sizeof(Bucket[0]) * blocks);

#ifdef BIG_RANGE
	assert(BucketInfo.MaxBucket <= BucketInfo.CurBucket);
	const uint pix = (uint)(sqrtp / log((double)sqrtp) * (1 + 1.200 / log((double)sqrtp)));
	blocks += pix / WHEEL_SIZE + BucketInfo.MaxBucket;
#endif

	if (BucketInfo.CurBucket == blocks)
		BucketInfo.MaxBucket = blocks;

	for (uint i = 0; i < blocks; i += MEM_BLOCK) {
		allocWheelBlock(MEM_BLOCK);
	}

	return blocks;
}

static void initWheelMedium(const uint sieve_size, const uint maxp, const uint64 start)
{
	uint j = Threshold.L1Index;
	const uint pix = (uint)(maxp / log((double)maxp) * (1 + 1.200 / log((double)maxp)));
	MediumWheel = (WheelPrime*) malloc(sizeof(WheelPrime) * pix + 1000);
	uint64 offset = start;

	for (uint p = Threshold.L1Maxp; p < maxp; p += Prime[++j]) {
		uint l2size = Threshold.L2Size;
		if ((sieve_size < l2size || p >= Threshold.L2Maxp) && l2size != sieve_size) {
			l2size = sieve_size;
			offset = start;
		}
		const uint64 p2 = (uint64)p * p;
		while (p2 >= offset + l2size && p2 > offset) {
			offset += l2size;
		}

		uint sieve_index = p - (uint)(offset % p);
		if (p2 > offset) {
			sieve_index = p2 - offset;
		}

		const uint pi = WHEEL_INIT[p % WHEEL].PrimeIndex;
		WheelFirst wf = WHEEL_FIRST[(sieve_index + offset) % WHEEL][pi];
		sieve_index += wf.Correct * p;
//		assert(sieve_index / WHEEL30 < (-1u >> SIEVE_BIT));
		MediumWheel[j].SieveIndex = (sieve_index / WHEEL30 << SIEVE_BIT) | wf.WheelIndex;
		MediumWheel[j].Wp = (p / WHEEL << SIEVE_BIT) + pi;
	}

	MediumWheel[j].Wp = -1u;
}

static void pushBucket(const uint sieve_index, const uint wp, const uchar wheel_index)
{
	const uint next_bucket = sieve_index >> BucketInfo.Log2Size;
#ifndef BIG_RANGE
	if (next_bucket >= BucketInfo.CurBucket) {
		return;
	}
#endif

	_Bucket* pbucket = Bucket + next_bucket;
	WheelPrime* wheel = pbucket ->Wheel ++;
	if ((size_t)wheel % MEM_WHEEL == 0) {
		Stock* pstock = StockPool[-- BucketInfo.CurStock];
		wheel = pstock ->Wheel;
		pbucket ->Wheel = pstock ->Wheel + 1;
		pstock ->Next = pbucket ->Shead;
		pbucket ->Shead = pstock;

#ifndef BIG_RANGE
		if (BucketInfo.CurStock == 0)
			allocWheelBlock(MEM_BLOCK);
#endif
	}

	wheel ->Wp = wp;
	wheel ->SieveIndex = ((sieve_index & BucketInfo.SieveSize) << SIEVE_BIT) | wheel_index;
}

static int segmentedSieve(uint64 start, uint sieve_size, Cmd*);
static void initBucketWheel(uint medium, uint sqrtp, const uint64 start, const uint64 range)
{
	uint nextp = 0; uint64 remp = 0;
	if (++sqrtp == 0) sqrtp = -1u; //overflow

	for (uint segsize = L2_DCACHE_SIZE * (WHEEL30 << 10); medium < sqrtp; medium += segsize) {
		if (segsize > sqrtp - medium)
			segsize = sqrtp - medium;

		Cmd cmd = {COPY_BITS, 0, 0};
		segmentedSieve(medium, segsize, &cmd);

		stype mask = 0, *bitarray = (stype*)cmd.Data;
		uint offset = medium - medium % WHEEL30 - sizeof(stype) * WHEEL30;
		const uint pmax = (uint64)medium * medium / (uint)(start >> 32);
		const uint pn = 2 + ((int)cmd.Primes) / sizeof(stype);

		for (uint j = 0; j < pn; ) {
			if (mask == 0) {
				mask = ~bitarray[j ++];
				offset += sizeof(stype) * WHEEL30;
				continue;
			}

			const uint p = offset + PRIME_OFFSET(mask);
			mask &= mask - 1;
#if ASM_X86
			if (p > nextp) {
				remp = start / (nextp = p + pmax);
				if (p > nextp)
					remp = start >> 32, nextp = -1u;
			}
			uint sieve_index = p - fastMod(start - remp * p, p);
#else
			uint sieve_index = p - (uint)(start % p);
#endif

#ifndef BIG_RANGE
			if (sieve_index > range) {
				continue;
			}
#endif
			const uint wp = (p / WHEEL210 << WHEEL_BIT) + WheelInit210[p % WHEEL210].PrimeIndex;
			const uint module_wheel = sieve_index % WHEEL210;
			const WheelFirst& wf = WheelFirst210[module_wheel][wp % (1 << WHEEL_BIT)];

#if X86_64
			sieve_index = (sieve_index + wf.Correct * (uint64)p) / WHEEL30;
#else
			sieve_index = (sieve_index / WHEEL210 + (wp >> WHEEL_BIT) * wf.Correct) * (WHEEL210 / WHEEL30);
			sieve_index += (wf.Correct * (p % WHEEL210) + module_wheel) / WHEEL30;
#endif
			pushBucket(sieve_index, wp, wf.WheelIndex);
		}
	}
}

inline static void
sieveSmall0(uchar bitarray[], const uchar* pend, const uint p, uint sieve_index, ushort multiples)
{
	uchar* ppbeg[8];
	for (int i = 0; i < 8; i ++) {
		const uchar windex = WheelInit30[sieve_index % WHEEL30].WheelIndex;
		ppbeg[windex] = bitarray + sieve_index / WHEEL30;
		sieve_index += (multiples % 4) * 2 * p; multiples /= 4;
	}
	crossOffWheelFactor(ppbeg, pend, p);
}

inline static void
sieveSmall1(uchar bitarray[], const uchar* pend, const uint p, uint sieve_index, ushort multiples)
{
	for (int i = 0; i < 4; i ++) {
		uchar* ps0 = bitarray + sieve_index / WHEEL30;
		uchar masks0 = WheelInit30[sieve_index % WHEEL30].UnsetBit;
		sieve_index += (multiples % 4) * 2 * p; multiples /= 4;

		uchar* ps1 = bitarray + sieve_index / WHEEL30;
		uchar masks1 = WheelInit30[sieve_index % WHEEL30].UnsetBit;
		sieve_index += (multiples % 4) * 2 * p; multiples /= 4;
		crossOff2Factor(ps0, ps1, pend, masks0 | (masks1 << 8), p);
	}
}

inline static void
sieveSmall2(uchar bitarray[], const uchar* pend, const uint p, uint sieve_index, ushort multiples)
{
	uchar* ppbeg[8];
	uint64 mask = 0;
	for (int i = 0; i < 8; i ++) {
		ppbeg[i] = bitarray + sieve_index / WHEEL30;
		mask = (mask << 8) | WheelInit30[sieve_index % WHEEL30].UnsetBit;
		sieve_index += (multiples % 4) * 2 * p; multiples /= 4;
	}
	crossOff8Factor(ppbeg, pend, mask, p);
}

inline static void
sieveSmall3(uchar bitarray[], const uchar* pend, const uint p, uint sieve_index, ushort multiples)
{
	for (int i = 0; i < 8; i ++) {
		const uchar mask = WheelInit30[sieve_index % WHEEL30].UnsetBit;
		if (mask == BIT0)
			break;
		bitarray[sieve_index / WHEEL30] |= mask;
		sieve_index += (multiples % 4) * 2 * p; multiples /= 4;
	}
	crossOffWheelFactor2(bitarray + sieve_index / WHEEL30, pend + 1, p);
}

static void eratSieveL1(uchar bitarray[], const uint64 start, const int segsize, uint maxp)
{
	const uchar* pend = bitarray + segsize / WHEEL30;
	if ((start + segsize) < ((uint64)maxp) * maxp) {
		maxp = isqrt(start + segsize) + 1;
	}

	for (uint p = Prime[0], j = FIRST_INDEX; p < maxp; p += Prime[++j]) {
		uint sieve_index = p - (uint)(start % p);
		if (start <= p) {
			sieve_index = p * p - (uint)start;
		}

		const WheelFirst wf = WheelFirst30[sieve_index % WHEEL30][WheelInit30[p % WHEEL30].PrimeIndex];
		const ushort multiples = WHEEL_SKIP >> (wf.NextMultiple * 2);
		sieve_index += wf.Correct * p;
		sieveSmall3(bitarray, pend, p, sieve_index, multiples);
	}
}

static void eratSieveL0(uchar bitarray[], const uint64 start, const int segsize, uint maxp)
{
	const uchar* pend = bitarray + segsize / WHEEL30;
	if ((start + segsize) < ((uint64)maxp) * maxp) {
		maxp = isqrt(start + segsize) + 1;
	}

	for (uint p = Prime[0], j = FIRST_INDEX; p < maxp; p += Prime[++j]) {
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

//Pre-sieve multiples of small primes <= limit (e.g. 19).
//Resets the sieve array (resets bits to 1) of SieveOfEratosthenes
//objects after each sieved segment and removes the multiples of small primes without sieving.
static void preSieve(uchar bitarray[], const uint64 start, const uint sieve_size)
{
	const uint offset = (uint)(start % PRIME_PRODUCT) / WHEEL30;
	const uint bits = sieve_size / WHEEL30 * 8 + WheelInit30[sieve_size % WHEEL30].WheelIndex;
	const uint bytes = (bits + 7) / 8;

	if (offset + bytes < sizeof(PreSieved)) {
		memcpy(bitarray, PreSieved + offset, bytes);
	} else {
		memcpy(bitarray, PreSieved + offset, sizeof(PreSieved) - offset);
		memcpy(bitarray + sizeof(PreSieved) - offset, PreSieved, bytes + offset - sizeof(PreSieved));
	}

	//1 is not prime, other pattern < 30 is prime
	if (start == 0) {
		bitarray[0] = BIT0;
	}
	//set the last byte bit 1
	if (bits & 7)
		bitarray[bits >> 3] |= ~((1 << (bits & 7)) - 1);
}

static void eratSieveSmall(uchar bitarray[], const uint64 start, const int sieve_size)
{
	#pragma omp parallel for if (sieve_size > 2 * Threshold.L1Size)
	for (int sieve_index = 0; sieve_index < sieve_size; sieve_index += Threshold.L1Size) {
		int segsize = Threshold.L1Size;
		if (segsize + sieve_index > sieve_size)
			segsize = sieve_size - sieve_index;
		preSieve(bitarray + sieve_index / WHEEL30, start + sieve_index, segsize);
		eratSieveL0(bitarray + sieve_index / WHEEL30, start + sieve_index, segsize, Threshold.L1Maxp);
	}
}

#if WHEEL == WHEEL30
inline static WheelElement*
sieveMediumWheel(uchar bitarray[], const uchar* pend, const uint wi, const uint pi, WheelElement* wheel)
{
	const uint p = wi * WHEEL30 + Pattern30[pi];
	WheelElement* wdata = Wheel30[pi];

	for (int i = 0; i < 4; i ++) {
		uchar* ps0 = bitarray;
		ushort smask = wheel ->UnsetBit;
		bitarray += wheel ->Correct + wheel ->NextMultiple * wi;
		wheel = wdata + wheel ->WheelIndex;

		uchar* ps1 = bitarray;
		smask |= ((ushort)wheel ->UnsetBit) << 8;
		bitarray += wheel ->Correct + wheel ->NextMultiple * wi;
		if (i < 3)
			wheel = wdata + wheel ->WheelIndex;

		crossOff2Factor(ps0, ps1, pend, smask, p);
	}

	return wheel;
}
#endif

#if SAFE
	#define SAFE_SET(n) \
		we##n = wdata##n[we##n.WheelIndex]; \
		bitarray[(int)sieve_index##n] |= we##n.UnsetBit; \
		sieve_index##n += we##n.Correct + (we##n.NextMultiple) * wi##n
#else
	#define SAFE_SET(n) \
		wheel##n = *(uint*)(wdata##n + (uchar)wheel##n); \
		bitarray[(int)sieve_index##n] |= wheel##n >> 8; \
		sieve_index##n += (uchar)(wheel##n >> 16) + (wheel##n >> 24) * wi##n
#endif

//sieve one medium prime from array
inline static void
sieveMediumOne(uchar bitarray[], const uint sieve_byte, WheelPrime* pwheel)
{
	uint& wheel = pwheel ->SieveIndex;
	int sieve_index = (wheel >> SIEVE_BIT) - sieve_byte;

	if (sieve_index >= 0) {
		wheel -= sieve_byte << SIEVE_BIT;
		return;
	}

	const uint wi = pwheel ->Wp >> SIEVE_BIT;
	WheelElement* wdata = WHEEL_MAP[pwheel ->Wp % (1 << SIEVE_BIT)];

#if SAFE
	WheelElement we; we.WheelIndex = wheel % (1 << SIEVE_BIT);
#endif

	do {
		SAFE_SET();
	} while (sieve_index < 0);

#if SAFE
	wheel = (sieve_index << SIEVE_BIT) | we.WheelIndex;
#else
	wheel = (sieve_index << SIEVE_BIT) | (uchar)wheel;
#endif
}

//sieve two medium prime from array
inline static void
sieveMediumTwo(uchar bitarray[], const uint sieve_byte, WheelPrime* pwheel)
{
	uint wheel1 = pwheel[0].SieveIndex, wp1 = pwheel[0].Wp;
	uint wi1 = wp1 >> SIEVE_BIT, sieve_index1 = (wheel1 >> SIEVE_BIT) - sieve_byte;
	WheelElement* wdata1 = WHEEL_MAP[wp1 % (1 << SIEVE_BIT)];

	uint wheel2 = pwheel[1].SieveIndex, wp2 = pwheel[1].Wp;
	uint wi2 = wp2 >> SIEVE_BIT, sieve_index2 = (wheel2 >> SIEVE_BIT) - sieve_byte;
	WheelElement* wdata2 = WHEEL_MAP[wp2 % (1 << SIEVE_BIT)];

#if SAFE
	WheelElement we1; we1.WheelIndex = wheel1 % (1 << SIEVE_BIT);
	WheelElement we2; we2.WheelIndex = wheel2 % (1 << SIEVE_BIT);
#endif

	while ((int)sieve_index1 < 0) {
		SAFE_SET(1);
		if ((int)sieve_index2 < 0) {
			SAFE_SET(2);
		}
	}

	while ((int)sieve_index2 < 0) {
		SAFE_SET(2);
	}

#if SAFE
	wheel1 = we1.WheelIndex, wheel2 = we2.WheelIndex;
#endif

	pwheel[0].SieveIndex = sieve_index1 << SIEVE_BIT | (uchar)wheel1;
	pwheel[1].SieveIndex = sieve_index2 << SIEVE_BIT | (uchar)wheel2;
}

//core code of this algorithm for large range
//sieve prime multiples in [start, start + sieve_size)
static void eratSieveMedium(uchar bitarray[], const uint64 start, const uint sieve_size, const uint whi, uint maxp)
{
	if ((start + sieve_size) < (uint64)maxp * maxp) {
		maxp = isqrt(start + sieve_size) + 1;
	}

	const uint sieve_byte = sieve_size / WHEEL30 + WheelInit30[sieve_size % WHEEL30].WheelIndex;
	WheelPrime* pwheel = MediumWheel + whi;

#if WHEEL == WHEEL30
	const uchar* pend = bitarray + sieve_byte + 1;
	uint minp = MIN(maxp, sieve_byte / 2);
	minp = (minp / WHEEL << SIEVE_BIT) + WHEEL_INIT[minp % WHEEL].PrimeIndex;

	for (uint wp1 = pwheel ->Wp; wp1 < minp; wp1 = pwheel ->Wp) {
		const uint wi = wp1 >> SIEVE_BIT, pi = wp1 % (1 << SIEVE_BIT);
		const uint p = wi * WHEEL30 + Pattern30[pi];
		uint sieve_index = pwheel ->SieveIndex >> SIEVE_BIT;
		WheelElement* wdata = Wheel30[pi];
		WheelElement* wheel = wdata + pwheel ->SieveIndex % (1 << SIEVE_BIT);
		wheel = sieveMediumWheel(bitarray + sieve_index, pend, wi, pi, wheel);

		sieve_index += (sieve_byte - sieve_index) / p * p;
		while (sieve_index < sieve_byte) {
			wheel = wdata + wheel ->WheelIndex;
			sieve_index += wheel ->Correct + wheel ->NextMultiple * wi;
		}
		pwheel ++ ->SieveIndex = (sieve_index - sieve_byte) << SIEVE_BIT | wheel ->WheelIndex;
	}
#endif

	maxp = (maxp / WHEEL << SIEVE_BIT) + WHEEL_INIT[maxp % WHEEL].PrimeIndex;
	bitarray += sieve_byte;

#if 1
	while (pwheel[1].Wp < maxp) {
		sieveMediumTwo(bitarray, sieve_byte, pwheel), pwheel += 2;
	}
	if (pwheel -> Wp < maxp) {
		sieveMediumOne(bitarray, sieve_byte, pwheel);
	}
#else
	while (pwheel -> Wp < maxp) {
		sieveMediumOne(bitarray, sieve_byte, pwheel ++);
	}
#endif
}

//sieve one big prime from bucket
static void sieveBigOne(uchar bitarray[], const uint sieve_size, WheelPrime* cur_wheel)
{
	uint sieve_index = cur_wheel ->SieveIndex, wp = cur_wheel ->Wp;
	const uint wi = wp >> WHEEL_BIT;
	WheelElement* wdata = Wheel210[wp % (1 << WHEEL_BIT)];

#if SAFE
	WheelElement* wheel = &wdata[sieve_index % (1 << SIEVE_BIT)];
	bitarray[sieve_index >>= SIEVE_BIT] |= wheel ->UnsetBit;
	sieve_index += wheel ->Correct + wheel ->NextMultiple * wi;
#if ERAT_BIG > 2
	if (sieve_index < sieve_size) {
		wheel = wdata + wheel ->WheelIndex;
		bitarray[sieve_index] |= wheel ->UnsetBit;
		sieve_index += wheel ->Correct + wheel ->NextMultiple * wi;
	}
#endif
	pushBucket(sieve_index, wp, wheel ->WheelIndex);
#else
	uint wheel = *(uint*)(wdata + sieve_index % (1 << SIEVE_BIT));
	bitarray[sieve_index >>= SIEVE_BIT] |= wheel >> 8;
	sieve_index += (uchar)(wheel >> 16) + (wheel >> 24) * wi;

#if ERAT_BIG > 2
	if (sieve_index < sieve_size) {
		wheel = *(uint*)(wdata + (uchar)wheel);
		bitarray[sieve_index] |= wheel >> 8;
		sieve_index += (uchar)(wheel >> 16) + (wheel >> 24) * wi;
	}
#endif
	pushBucket(sieve_index, wp, wheel);
#endif
}

//sieve two big prime from bucket, 15% improvment
static void sieveBigTwo(uchar bitarray[], const uint sieve_size, const WheelPrime* cur_wheel)
{
	uint sieve_index1 = cur_wheel[0].SieveIndex, wp1 = cur_wheel[0].Wp;
	uint sieve_index2 = cur_wheel[1].SieveIndex, wp2 = cur_wheel[1].Wp;

#if SAFE
	WheelElement* wheel1 = &Wheel210[wp1 % (1 << WHEEL_BIT)][sieve_index1 % (1 << SIEVE_BIT)];
	WheelElement* wheel2 = &Wheel210[wp2 % (1 << WHEEL_BIT)][sieve_index2 % (1 << SIEVE_BIT)];
	bitarray[sieve_index1 >>= SIEVE_BIT] |= wheel1 ->UnsetBit;
	sieve_index1 += wheel1 ->Correct + wheel1 ->NextMultiple * (wp1 >> WHEEL_BIT);

	bitarray[sieve_index2 >>= SIEVE_BIT] |= wheel2 ->UnsetBit;
	sieve_index2 += wheel2 ->Correct + wheel2 ->NextMultiple * (wp2 >> WHEEL_BIT);

#if ERAT_BIG > 2
	if (sieve_index1 < sieve_size) {
		wheel1 = &Wheel210[wp1 % (1 << WHEEL_BIT)][wheel1 ->WheelIndex];
		bitarray[sieve_index1] |= wheel1 ->UnsetBit;
		sieve_index1 += wheel1 ->Correct + wheel1 ->NextMultiple * (wp1 >> WHEEL_BIT);
	}
	if (sieve_index2 < sieve_size) {
		wheel2 = &Wheel210[wp2 % (1 << WHEEL_BIT)][wheel2 ->WheelIndex];
		bitarray[sieve_index2] |= wheel2 ->UnsetBit;
		sieve_index2 += wheel2 ->Correct + wheel2 ->NextMultiple * (wp2 >> WHEEL_BIT);
	}
#endif

	pushBucket(sieve_index1, wp1, wheel1 ->WheelIndex);
	pushBucket(sieve_index2, wp2, wheel2 ->WheelIndex);
#else
	uint wheel1 = *(uint*)&Wheel210[wp1 % (1 << WHEEL_BIT)][sieve_index1 % (1 << SIEVE_BIT)];
	uint wheel2 = *(uint*)&Wheel210[wp2 % (1 << WHEEL_BIT)][sieve_index2 % (1 << SIEVE_BIT)];
	bitarray[sieve_index1 >>= SIEVE_BIT] |= wheel1 >> 8;
	sieve_index1 += (uchar)(wheel1 >> 16) + (wheel1 >> 24) * (wp1 >> WHEEL_BIT);

	bitarray[sieve_index2 >>= SIEVE_BIT] |= wheel2 >> 8;
	sieve_index2 += (uchar)(wheel2 >> 16) + (wheel2 >> 24) * (wp2 >> WHEEL_BIT);

#if ERAT_BIG > 2
	if (sieve_index1 < sieve_size) {
		wheel1 = *(uint*)&Wheel210[wp1 % (1 << WHEEL_BIT)][(uchar)wheel1];
		bitarray[sieve_index1] |= wheel1 >> 8;
		sieve_index1 += (uchar)(wheel1 >> 16) + (wheel1 >> 24) * (wp1 >> WHEEL_BIT);
	}
	if (sieve_index2 < sieve_size) {
		wheel2 = *(uint*)&Wheel210[wp2 % (1 << WHEEL_BIT)][(uchar)wheel2];
		bitarray[sieve_index2] |= wheel2 >> 8;
		sieve_index2 += (uchar)(wheel2 >> 16) + (wheel2 >> 24) * (wp2 >> WHEEL_BIT);
	}
#endif

	pushBucket(sieve_index1, wp1, wheel1);
	pushBucket(sieve_index2, wp2, wheel2);
#endif
}

//This implementation uses a sieve array with WHEEL210 numbers per byte and
//a modulo 210 wheel that skips multiples of 2, 3, 5 and 7.
static void eratSieveBucket(uchar bitarray[], const uint sieve_size)
{
	uint loops = (size_t)Bucket[0].Wheel % MEM_WHEEL / sizeof(WheelPrime);
	if (loops % 2) {
//		sieveBigOne(bitarray, sieve_size, Bucket[0].Wheel - 1), loops --;
		*(Bucket[0].Wheel) = *(Bucket[0].Wheel - 1), loops++;
	} else if (loops == 0) {
		loops = WHEEL_SIZE;
	}

	for (Stock* phead = (Stock*)Bucket; phead = phead ->Next; loops = WHEEL_SIZE) {
		WheelPrime* cur_wheel = phead ->Wheel;
		while (loops) {
#ifndef SBIG2
			sieveBigTwo(bitarray, sieve_size, cur_wheel), loops -= 2, cur_wheel += 2;
#else
			sieveBigOne(bitarray, sieve_size, cur_wheel ++), loops --;
#endif
		}
		StockPool[BucketInfo.CurStock ++] = phead;
	}

	BucketInfo.CurBucket --;
	//rotate bucket list
	memmove(Bucket, Bucket + 1, BucketInfo.MaxBucket * sizeof(Bucket[0]));
}

static int segmentProcessed(uchar bitarray[], const uint64 start, const uint bytes, Cmd* cmd)
{
	if (NULL == cmd)
		return countBit0sArray((uint64*)bitarray, bytes * 8);

	const int oper = cmd ->Oper;
	int primes = 0, dwords = (sizeof(stype) + bytes - 1) / sizeof(stype);

	if (COPY_BITS == oper) {
		cmd ->Data = bitarray;
		primes = bytes;
	} else if (PCALL_BACK == oper) {
		primes = doCall((ushort*)bitarray, start, bytes, cmd ->Primes, (sieve_call)cmd ->Data);
	} else if (SAVE_BYTE == oper) {
		primes = savePrimeByte((stype*)bitarray, bytes, cmd ->Data + cmd ->Primes);
//	} else if (SAVE_PRIME == oper) {
//		primes = savePrime((stype*)bitarray, start, dwords, (uint64*)cmd ->Data + cmd ->Primes);
	}

	cmd ->Primes += primes;
	return primes;
}

//sieve prime multiples in [start, start + sieve_size)
static int segmentedSieve(uchar bitarray[], const uint64 start, const uint wheel_offset, const uint sieve_size)
{
	const uint sqrtp = isqrt(start + sieve_size);
	const uint medium = MIN(sqrtp, Threshold.Medium);
	const uint bytes = sieve_size / WHEEL30 + (7 + WheelInit30[sieve_size % WHEEL30].WheelIndex) / 8;

	for (uint sieve_index = 0; sieve_index < sieve_size; sieve_index += Threshold.L2Size) {
		uint segsize = Threshold.L2Size;
		if (segsize + sieve_index > sieve_size)
			segsize = sieve_size - sieve_index;

		eratSieveSmall(bitarray + sieve_index / WHEEL30, start + sieve_index, segsize);

#if L2_DCACHE_SIZE != MAX_SIEVE
		//static int64 time_use = 0; int64 ts = getTime();
		if (MediumWheel)
			eratSieveMedium(bitarray + sieve_index / WHEEL30, start + sieve_index, segsize, Threshold.L1Index, Threshold.L2Maxp);
		//time_use += getTime() - ts; if (sieve_size != Config.SieveSize) { printf("eratSieveMedium1 time %.f ms\n", time_use); time_use = 0; }
#endif
	}

	//static int64 time_use1 = 0; int64 ts1 = getTime(); //280 ms
#if L2_DCACHE_SIZE != MAX_SIEVE
	if (medium >= Threshold.L2Maxp)
		eratSieveMedium(bitarray, start, sieve_size, Threshold.L2Index, medium + 1);
#else
	if (medium >= Threshold.L1Maxp)
		eratSieveMedium(bitarray, start, sieve_size, Threshold.L1Index, medium + 1);
#endif

	//time_use1 += getTime() - ts1; if (sieve_size != Config.SieveSize) { printf("eratSieveMedium time %lld ms\n", time_use1); time_use1 = 0; }
	if (start >= Config.MinBucket) {
		//static int64 time_use2 = 0; int64 ts2 = getTime(); //750
		eratSieveBucket(bitarray, Config.SieveSize / WHEEL30);
		//time_use2 += getTime() - ts2; if (sieve_size != Config.SieveSize) { printf("eratSieveBucket time %lld ms\n", time_use2); time_use2 = 0; }
	}

	*(uint64*)(bitarray + bytes) = ~0;
	if (wheel_offset > 0) {
		memset(bitarray, ~0, wheel_offset / WHEEL30);
		bitarray[wheel_offset / WHEEL30] |= (1 << WheelInit30[wheel_offset % WHEEL30].WheelIndex) - 1;
	}

	return bytes;
}

static int segmentedSieve(uint64 start, uint sieve_size, Cmd* cmd = NULL)
{
#if CHECK
	static uchar bitarray[MAX_CACHE];
#else
	uchar bitarray[(L1_DCACHE_SIZE + L2_DCACHE_SIZE) << 10];
#endif

	const uint wheel_offset = (uint)(start % WHEEL30);
	const uint sqrtp = isqrt(start + sieve_size);
	start -= wheel_offset, sieve_size += wheel_offset;
	const uint bytes = sieve_size / WHEEL30 + (7 + WheelInit30[sieve_size % WHEEL30].WheelIndex) / 8;

#if CHECK
	const uint mins = MIN(sqrtp, sieve_size);
#else
	const uint mins = sqrtp;
#endif

//	#pragma omp parallel for if (sieve_size > 2 * Threshold.L1Size)
	for (int sieve_index = 0; sieve_index < sieve_size; sieve_index += Threshold.L1Size) {
		int segsize = Threshold.L1Size;
		if (segsize + sieve_index > sieve_size)
			segsize = sieve_size - sieve_index;
		preSieve(bitarray + sieve_index / WHEEL30, start + sieve_index, segsize);
		eratSieveL1(bitarray + sieve_index / WHEEL30, start + sieve_index, segsize, Threshold.L1Maxp);
	}

	//	eratSieveSmall(bitarray, start, sieve_size);

	const uchar* pend = bitarray + sieve_size / WHEEL30;
	uint j = Threshold.L1Index, p = Threshold.L1Maxp;
	for (; p <= mins; p += Prime[++j]) {
		uint sieve_index = p - (uint)(start % p);
		if (start <= p) {
			sieve_index = p * p - (uint)start;
		}

		const WheelFirst wf = WheelFirst30[sieve_index % WHEEL30][WheelInit30[p % WHEEL30].PrimeIndex];
		const ushort multiples = WHEEL_SKIP >> (wf.NextMultiple * 2);
		sieve_index += wf.Correct * p;
		if (sieve_index <= sieve_size + WHEEL30)
			sieveSmall1(bitarray, pend, p, sieve_index, multiples);
	}

#if CHECK
	uint64 remp = 0;
	for (uint nextp = 0; p <= sqrtp && Prime[j]; p += Prime[++j]) {
		if (p % 2 == 0) p += 255;
#if 1
		if (p > nextp) {
			remp = start / (nextp = p + (uint64)p * p / (uint)(start >> 32));
			if (p > nextp)
				remp = start >> 32, nextp = -1u;
		}
		uint sieve_index = p - fastMod(start - remp * p, p);
#else
		uint sieve_index = p - (uint)(start % p);
#endif
		if (sieve_index <= sieve_size) {
			bitarray[sieve_index / WHEEL30] |= WheelInit30[sieve_index % WHEEL30].UnsetBit;
		}
	}
#endif

	bitarray[0] |= (1 << WheelInit30[wheel_offset].WheelIndex) - 1;
	*(uint64*)(bitarray + bytes) = ~0;

	return segmentProcessed(bitarray, start, bytes, cmd);
}

static void initThreshold()
{
	Threshold.L1Index = 0;

	uint p = Prime[0], j = FIRST_INDEX;
	for (; p < Threshold.L2Maxp && Prime[j]; p += Prime[++j]) {
		if (p >= Threshold.L1Maxp && Threshold.L1Index == 0) {
			Threshold.L1Index = j, Threshold.L1Maxp = p;
		}
	}
	if (p >= Threshold.L2Maxp) {
		Threshold.L2Index = j, Threshold.L2Maxp = p;
	}
}

void setCpuSize(uint cdata)
{
	cdata = 1 << ilog(cdata, 2);

	if (cdata >= 16 && cdata < L2_DCACHE_SIZE) { //L1
		Threshold.L1Size = cdata * (WHEEL30 << 10);
		Threshold.L1Maxp = Threshold.L1Size / (WHEEL30 * Threshold.L1Segs);
		initThreshold();
	} else if (cdata >= L2_DCACHE_SIZE && cdata <= MAX_SIEVE) { //L2
		Threshold.L2Size = cdata * (WHEEL30 << 10);
		Threshold.L2Maxp = Threshold.L2Size / (WHEEL30 * Threshold.L2Segs);
		initThreshold();
	}
}

void setLevelSegs(uint cdata)
{
	const int level = cdata / 10 % 10;
	const int segs = cdata % 10;
	if (segs > 1) {
		if (level == 1) {
			Threshold.L1Segs = segs;
			Threshold.L1Maxp = Threshold.L1Size / (WHEEL30 * segs);
			initThreshold();
		} else if(level == 2) {
			Threshold.L2Segs = segs;
			Threshold.L2Maxp = Threshold.L2Size / (WHEEL30 * segs);
			initThreshold();
		} else if (level == 3 && ERAT_BIG > 2 && segs <= 6) {
			Threshold.Msegs = segs;
		}
	}
}

int setSieveSize(uint sieve_size)
{
	if (sieve_size <= MAX_SIEVE && sieve_size > 16) {
		sieve_size *= (WHEEL30 << 10);
	} else {
		sieve_size = L2_DCACHE_SIZE * (WHEEL30 << 10);
	}

	sieve_size = WHEEL30 << ilog(sieve_size / WHEEL30, 2);

	return Config.SieveSize = sieve_size;
}

static int checkSmall(const uint64 start, const uint64 end, Cmd* cmd)
{
	int primes = 0;
	const uchar SmallPrime[] = {3, 5, 7, 0, 0};
	for (int i = 0, p = SmallPrime[i]; i < sizeof(SmallPrime) / sizeof(SmallPrime[0]); p = SmallPrime[++i]) {
		if (start <= p &&
			(p == SmallPrime[i + 1] - PRIME_GAP || p == SmallPrime[i + 2] - PRIME_GAP) &&
			SmallPrime[i + 1] <= end) {
			primes ++;
			if (cmd && cmd ->Oper == PCALL_BACK) {
				(*(sieve_call)cmd ->Data)(primes, p);
				cmd ->Primes += 1;
			}
		}
	}

	return primes;
}

//calculate number of prime in Range[start, end] with start <= end
static uint64 pi(const uint64 start, const uint64 end, Cmd* cmd = NULL)
{
	const int64 ts = getTime();
	uint sieve_size = Config.SieveSize;
	int64 primes = checkSmall(start, end, cmd);
	if (cmd == NULL && end > (uint64)sieve_size * sieve_size) {
		sieve_size = MAX_SIEVE * (WHEEL30 << 10);
	}

	if (end - start <= sieve_size) {
		primes += segmentedSieve(start, (uint)(end - start) + 1, cmd);
	} else {
		//bugs
		primes += segmentedSieve(start, sieve_size - (uint)(start % sieve_size), cmd);
		uint64 newstart = start - start % sieve_size + sieve_size, newend = end - sieve_size;
		#pragma omp parallel for reduction(+:primes) schedule(static, 2) if(cmd == NULL)
		for (uint64 offset = newstart; offset < newend; offset += sieve_size) {
			primes += segmentedSieve(offset, sieve_size, cmd);
#if CHECK
			if (primes % 16 == 0) {
				double ratio = 1000.0 * (int64)(offset + sieve_size - start) / (int64)(end - start);
				printf(">> %.2f%%, time(%.2f sec), primes ~= %llu\r",
						ratio / 10, (getTime() - ts) / ratio, (int64)(1000 * primes / ratio));
			}
#endif
		}
		primes += segmentedSieve(end - (uint)(end % sieve_size), (uint)(end % sieve_size) + 1, cmd);
	}

	return primes;
}

static uint64 pi2(uint64 start, uint64 end, Cmd* cmd)
{
	const int64 ts = getTime();

	uchar* bitarray = (uchar*) malloc(MAX_CACHE);
	const int64 range = (int64)(end - start);

	uint wheel_offset = (uint)(start % WHEEL210);
	start -= wheel_offset;
	int64 primes = checkSmall(start, end, cmd);
	uint64 last_qword = ~0;

	if (++end == 0) end --; //fix overflow 2^64 - 1

	for (uint si = 0, sieve_size = Config.SieveSize; start < end; start += sieve_size) {
		if (sieve_size > end - start) {
			sieve_size = (uint)(end - start);
		}
		const uint bytes = segmentedSieve(bitarray + 8, start, wheel_offset, sieve_size);
		*(uint64*)bitarray = (last_qword << 56) | ~((uint64)1 << 63);
		primes += segmentProcessed(bitarray, start - 8 * WHEEL30, bytes + 8, cmd);
		last_qword = *(bitarray + 8 + bytes - 1);
		wheel_offset = 0;
		if ((si ++ & Config.Progress) == 1) {
			double ratio = 1000 - 1000.0 * ((int64)(end - start) - sieve_size) / (int64)range;
			double timeuse = (getTime() - ts) / ratio;
			printf(">> %.2f%%, time(%.2f sec), primes ~= %llu\r", ratio / 10, timeuse, (int64)(1000 * primes / ratio));
		}
	}

	free(bitarray);

	return primes;
}

static void printPiResult(const uint64 start, const uint64 end, uint64 primes)
{
	const int sta10 = ilog(start, 10);
	const int end10 = ilog(end, 10);
	const int dif10 = ilog(end - start + 1, 10);

#if 0
	printf("pi2(%llu, %llu) = %llu", start, end, primes);
#else
	if (start > 0) {
		if (start % ipow(10, sta10) == 0 && sta10 > 2)
			printf("pi2(%de%d,", (int)(start / ipow(10, sta10)), sta10);
		else if ((start & (start - 1)) == 0)
			printf("pi2(2^%d,", ilog(start, 2));
		else
			printf("pi2(%llu,", start);

		if (end % ipow(10, end10) == 0)
			printf("%de%d)", (int) (end / ipow(10, end10)), end10);
		else if ((end - start) % ipow(10, dif10) == 0 && dif10 > 2) {
			if (start % ipow(10, sta10) == 0) {
				printf("%de%d+%de%d)", (int)(start / ipow(10, sta10)), sta10, (int)((end - start) / ipow(10, dif10)), dif10);
			} else if ((start & (start - 1)) == 0) {
				printf("2^%d+%de%d)", ilog(start, 2), (int)((end - start) / ipow(10, dif10)), dif10);
			} else {
				printf("%llu+%de%d)", start, (int)((end - start) / ipow(10, dif10)), dif10);
			}
		} else {
			printf("%llu)", end);
		}
	} else if (end % ipow(10, end10) == 0 && end10 > 2) {
		printf("pi2(%de%d)", (int) (end / ipow(10, end10)), end10);
	} else {
		printf("pi2(%llu)", end);
	}
	printf(" = %llu", primes);
#endif
}

static int sievePrime(uchar prime[], uint n)
{
	static uint maxp = 1 << 16;
	if (n <= maxp) {
		return 0;
	}
	maxp = n;

	Cmd cmd = {SAVE_BYTE, prime, sizeof(SmallPrime) / sizeof(SmallPrime[0])};

	const int sieve_size = Config.SieveSize;
	Config.SieveSize = L2_DCACHE_SIZE * (WHEEL30 << 10);
	int primes = (int)pi(0, n, &cmd);
	Config.SieveSize = sieve_size;

	prime[primes + 0] = 0;
	prime[primes + 1] = 200;

#if _DEBUG
	printf("pi(%u) = %d\n", n, primes);
	int pi = 8 + PRIME_PRODUCT / 9699690;
	uint p = Prime[0];
	//assert(prime[1] + prime[2] == 3);
	for (; prime[pi] > 0; ) {
		p += prime[pi++];
		if (p % 2 == 0) p += 255;
	}
	assert(primes == pi);
	printf("Prime[%d] = %u\n", pi, p);
#endif

	initThreshold();

	return primes;
}

static uint64 doSieve2(const uint64 start, const uint64 end, Cmd* cmd = NULL)
{
	const int64 ts = getTime();

	sievePrime(Prime + 1, isqrt(end));

	const int64 tini = getTime();

	const uint64 primes = pi(start, end, cmd);

	if (Config.Flag & PRINT_RET) {
		printPiResult(start, end, primes);
		if (Config.Flag & PRINT_TIME)
			printf(" (%.2f + %.2f sec)", (getTime() - tini) / 1000.0, (tini - ts) / 1000.0);
		putchar('\n');
	}

	return primes;
}

static uint initMediumSieve(const uint64 start, const uint sqrtp)
{
	const uint sieve_size = Config.SieveSize;
	uint medium = sieve_size / Threshold.Msegs + 1;
	uint64 offset = (uint64)medium * medium;
	//adjust medium
	if (sqrtp > medium && offset > start) {
		offset = offset - offset % sieve_size + start % sieve_size + sieve_size;
		while (offset % WHEEL210 != start % WHEEL210)
			offset += sieve_size;
		medium = isqrt((offset - offset % WHEEL210 + sieve_size)) + 1;
	} else {
		offset = start;
	}

	medium = MIN(sqrtp, medium);
	sievePrime(Prime + 1, medium + 1476);

	if (sqrtp >= Threshold.L1Maxp) {
		//medium < 19358325 or (ERAT_BIG > 3 || MAX_SIEVE < 2048)
		assert(medium * (26 + 2) / WHEEL30 < (-1u >> SIEVE_BIT));
		initWheelMedium(sieve_size, medium + 256, start - start % WHEEL210);
	}

	Config.MinBucket = offset - offset % WHEEL210;
	Threshold.Medium = medium;

	return medium;
}

uint64 doSieve(const uint64 start, const uint64 end, Cmd* cmd = NULL)
{
	const int64 ts = getTime();

	const uint sqrtp = isqrt(end);
	if (Config.SieveSize < Threshold.L2Size && sqrtp > 10000000) {
		setSieveSize(MAX_SIEVE);
	}

	const uint medium = initMediumSieve(start, sqrtp);
	if (sqrtp > medium) {
		initBucketInfo(Config.SieveSize, sqrtp, end - Config.MinBucket);
		initBucketWheel(medium, sqrtp, Config.MinBucket, end - Config.MinBucket);
	} else {
		Config.MinBucket = end;
	}

	const int64 its = getTime();
	const uint64 primes = pi2(start, end, cmd);

	if (Config.Flag & PRINT_RET) {
		printPiResult(start, end, primes);
		if (Config.Flag & PRINT_TIME)
			printf(" (%.2f + %.2f sec)", (getTime() - its) / 1000.0, (its - ts) / 1000.0);
		putchar('\n');
	}
	if (MediumWheel) {
		free(MediumWheel); MediumWheel = NULL;
	}
	for (uint i = 0; i < BucketInfo.PtrSize; i ++) {
		free(WheelPtr[i]);
	}

#ifndef BIG_RANGE
	assert(BucketInfo.StockSize == BucketInfo.CurStock);
#endif

	memset(&BucketInfo, 0, sizeof(BucketInfo));

	return primes;
}

//the sieve of Eratosthenes implementated by bit packing
//all prime less than 2^16 will be saved in prime buffer List
//Prime[0] is the first sieve prime, Prime[i] is the difference
//Prime[0] = 2, Prime[1] = 3 - 2, Prime[2] = 5 - 3;
//of the adjacent prime, Prime[i] = Prime[i] - Prime[i - 1];
static int eratoSimple()
{
	int primes = 2, lastp = 2;
	const uint maxp = 1 << 16;
	uchar bitarray[(maxp >> 4) + 8] = {0};

	for (uint p = 3; p <= maxp; p += 2) {
		if (0 == (bitarray[p >> 4] & (1 << (p / 2 & 7)))) {
			Prime[primes ++] = p - lastp;
			lastp = p;
			for (uint j = p * p / 2; j <= maxp / 2; j += p) {
				bitarray[j >> 3] |= 1 << (j & 7);
			}
		}
	}

	//pack the last two byte for safety
	Prime[primes + 1] = Prime[primes] = 254;

	return primes;
}

//The first presieved template, sieve the first prime multiples
static void initPreSieved( )
{
	const uchar prime[ ] = {7, 11, 13, 17, 19, 23, 29};
	for (int i = 0; PRIME_PRODUCT % prime[i] == 0; i ++) {
		Prime[0] = prime[i + 1];
		int p = prime[i];
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
	for (i = 0; i < nbitsize2; i += 2) {
		Lsb[i + 0] = Lsb[i >> 1] + 1;
	}

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
	int pattern210[WHEEL210] = {0};
	int wpattern[WHEEL210] = {0};

	for (i = 0; i < WHEEL210; i ++) {
		const uchar mask = WheelInit30[i % WHEEL30].UnsetBit;
		WheelInit210[i].UnsetBit = mask;
		WheelInit210[i].WheelIndex = -1;
		WheelInit210[i].PrimeIndex = wi;
		if (mask && i % (WHEEL210 / WHEEL30)) {
			pattern210[wi ++] = i;
		}
	}

	wi = 0, i = 1;
#if PRIME_GAP == 2
	WheelInit210[i].WheelIndex = wi;
	wpattern[wi ++] = i;
#endif

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
			int next = wpattern[wi] + pattern210[pi] * 2;
			while (WheelInit210[next % WHEEL210].WheelIndex < 0) {
				next += pattern210[pi] * 2;
				multiples += 2;
			}

			WheelElement& wdata = Wheel210[pi][wi];
			wdata.NextMultiple = multiples * (WHEEL210 / WHEEL30);
			wdata.WheelIndex = WheelInit210[next % WHEEL210].WheelIndex;
			wdata.Correct = next / WHEEL30 - wpattern[wi] / WHEEL30;
			wdata.UnsetBit = WheelInit210[wpattern[wi]].UnsetBit;
		}
	}

	for (i = 0; i < WHEEL210; i ++) {
		for (int pi = 0; pi < psize; pi ++) {
			int multiples = 0, next = i;
			if (i % 2 == 0) {
				multiples = 1;
				next += pattern210[pi];
			}

			while (WheelInit210[next % WHEEL210].WheelIndex < 0) {
				next += pattern210[pi] * 2;
				multiples += 2;
			}

			WheelFirst210[i][pi].WheelIndex = WheelInit210[next % WHEEL210].WheelIndex;
			WheelFirst210[i][pi].Correct = multiples;
		}
	}
}

static void fixRangeTest(uint64 lowerBound, const int64 range, uint64 Ret)
{
	uint llog10 = ilog(lowerBound, 10), rlog10 = ilog(range, 10);
	printf("Sieving the primes within (10^%u, 10^%u+10^%u) randomly\n", llog10, llog10, rlog10);
	uint64 primes = 0, upperBound = lowerBound + range;

	while (lowerBound < upperBound) {
		uint64 rd = rand() * rand();
		uint64 end = lowerBound + (rd * rand() * rd) % ipow(2, 32) + ipow(10, 4);
		end = end - end % WHEEL210 + 6;
		if (end > upperBound)
			end = upperBound;
		setSieveSize(rand() % MAX_SIEVE + 128);
		setLevelSegs(rand() % 5 + 12), setLevelSegs(rand() % 5 + 22);
		setLevelSegs(rand() % 3 + 34);
		primes += doSieve(lowerBound, end - 1);
		printf("Remaining chunk: %.2f%%\r", (int64)(upperBound - lowerBound) * 100.0 / range);
		lowerBound = end;
	}
	printf("Pi[10^%u, 10^%u+10^%u] = %llu\n", llog10, llog10, rlog10, primes);
	assert(primes == Ret);
}

static void startBenchmark( )
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
		110670758  // pi2(10^15, 10^15 + 10^11)
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
	for (int j = 12; j <= 19; j ++) {
		uint64 start = ipow(10, j), end = start + ipow(2, 32);
		setLevelSegs(rand() % 5 + 12), setLevelSegs(rand() % 5 + 22);
		setLevelSegs(rand() % 4 + 34);
		primes = doSieve(start, end);
		printf("pi2(10^%d, 10^%d+2^32) = %llu                 \n", j, j, primes);
		assert (primes == primeCounts[j]);
	}

	Config.Progress = 0;
	printf("Time elapsed %.f sec\n\n", (getTime() - ts) / 1000.0);
	puts("All Big tests passed SUCCESSFULLY!");

	const uint64 RangeData[][3] =
	{
		{ipow(10, 11), ipow(10, 11), 199708605},
		{ipow(10, 15), ipow(10, 11), 110670758},
		{ipow(10, 17), ipow(10, 11), 86176910},
		{ipow(10, 19), ipow(10, 11), 68985092},
		{ipow(10, 14), ipow(10, 12), 1270127074},
		{ipow(10, 16), ipow(10, 12), 972773783},
		{ipow(10, 18), ipow(10, 12), 768599834},
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
	__cpuid(cpuinfo, id);
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
		"-------------------------------------------------------------------------";
	puts(sepator);
	puts("Fast implementation of the segmented sieve of Eratosthenes (0, 2^64 - 1)\n"
	"Copyright @ by Huang Yuanbing (2011 - 2014) bailuzhou@163.com\n"
//	"Code: https://github.com/ktprime/ktprime/blob/master/PrimeNumber.cpp\n"
	"C++: g++ -march=native [-DW30,-DSAFE,-DPOPCNT] -funroll-loops -O3 -s -pipe");

#if _MSC_VER
	printf("Compiled by vc ++ %d", _MSC_VER);
#elif __GNUC__
	printf("Compiled by gcc %s", __VERSION__);
#endif

#if X86_64
	printf(" x86-64");
#endif

	printf(" %s %s\n", __TIME__, __DATE__);

	puts(sepator);
	printf("[MARCO] : ASM_X86, POPCNT, BIT_SCANF, SAFE = (%d, %d, %d, %d)\n", ASM_X86, POPCNT, BIT_SCANF, SAFE);
	printf("[MARCO] : MAX_SIEVE = %dk, WHEEL_SIZE = %dk, ERAT_BIG = %d, WHEEL = %d\n",
			MAX_SIEVE, WHEEL_SIZE >> 7, Threshold.Msegs, WHEEL);
	printf("[ARG ]  : L1Size = %dk, L2Size = %dk, SieveSize = %dk\n",
			Threshold.L1Size / WHEEL30 >> 10, Threshold.L2Size / WHEEL30 >> 10, Config.SieveSize / WHEEL30 >> 10);
	printf("[ARG ]  : L1Seg/L1Maxp/L2Seg/L2Maxp/Medium = (%d,%d,%d,%d,%d)\n",
		Threshold.L1Segs, Threshold.L1Maxp, Threshold.L2Segs, Threshold.L2Maxp, Threshold.Medium);
	puts(sepator);
	puts(sepator);
}

static void doCompile(const char* flag)
{
	char programming[257];
	strcpy(programming, __FILE__);
	char* pdot = strchr(programming, '.');
	if (pdot) {
		strcpy(pdot, "_.exe");
		puts(programming);
	}

	const char* const cxxflag =
#if _MSC_VER
		"cl /O2 /Oi /Ot /Oy /GT /GL %s %s %s";
#elif X86_64
		"g++ -m64 -march=native %s -O3 -funroll-loops -static -pipe %s -o %s";
#else
		"g++ -m32 -march=native %s -O3 -funroll-loops -static -pipe %s -o %s";
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
				puts(Help);
				break;
			case 'S':
				setSieveSize(cdata);
				break;
			case 'C':
				setCpuSize(cdata);
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
		while (isalnum(c) || c == '^' ||
				c == '+' || c == '-' || c == '*' || c == '=') {
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

		char params[64][60] = {0};

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
				for (int j = 12; j < 20; j ++) {
					uint64 start = ipow(10, j), range = ipow(10, 10);
					doSieve(start, start + range);
				}
				for (int i = 32; i < 64; i ++) {
					uint64 start = ipow(2, i), range = ipow(10, powi);
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
			if (Config.Flag & SLOW_TEST)
				doSieve2(start, end);
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
	static bool initOnce = true;
	if (initOnce) {
#if X86_64 || _M_IX86 || __i386__
		getCpuInfo();
#endif
		eratoSimple();
		initBitTable( );
		initWheel30( );
		initWheel210( );
		initPreSieved( );
		initThreshold();
		setSieveSize(sieve_size);
		initOnce = false;
	}
}

#ifndef PRIME_LIB
int main(int argc, char* argv[])
{
	initPrime(MAX_SIEVE);

	if (argc > 1)
		executeCmd(argv[1]);

#if 0
	srand(time(0));
	for (int j = 1; j < 20; j++) {
		for (int i = 2; i <= 6; i++) {
			setSieveSize(MAX_SIEVE / 2);
			setLevelSegs(i + 30); setLevelSegs(rand() % 5 + 12), setLevelSegs(rand() % 5 + 22);
			const uint64 medium = Config.SieveSize / i;
			const uint64 start = medium * medium;
			const uint64 range = ipow(10, 9) + (uint)rand() * rand();
			const uint64 r1 = doSieve(start - range, start + range, NULL);

			setSieveSize(L2_DCACHE_SIZE * 2);
			setLevelSegs(6); setLevelSegs(rand() % 5 + 12), setLevelSegs(rand() % 5 + 22);
			assert(r1 == doSieve(start - range, start + range, NULL));
		}
	}
#endif

#if GCOV_TEST
	executeCmd("m5 s1024 l14; 1e12 1e12+1e9"); putchar('\n');
//	executeCmd("1e18 1e9"); putchar('\n');
	executeCmd("c32 l22 1e16 10^9"); putchar('\n');
	executeCmd("l34 G 1e12 e8"); putchar('\n');
	executeCmd("da s256 1e5; 1e8"); putchar('\n');
	executeCmd("s4000 c2000 p 0 100"); putchar('\n');

	executeCmd("da C32 s1024 H I 1e18 1e8"); putchar('\n');
	executeCmd("df p 1e12+100 2e2*2; c1222222; g 1e16 1e8; s256 l32 l24 l14 1e8"); putchar('\n');
	executeCmd("df s1024");
	executeCmd("0-1e4 0-1");
#endif

#ifndef BIG_RANGE
	executeCmd("1e16 1e6 s0248");
	executeCmd("1e12 1e9; 1e16 1e9; e18 e9");
#else
	executeCmd("e16 1e10");
#endif

	while (true) {
		char ccmd[257];
		printf("\n>> ");
		if (!gets(ccmd) || !executeCmd(ccmd))
			break;
	}

	return 0;
}
#endif
