//Copyright @ by Huang Yuanbing 2011 - 2012 bailuzhou AT 163.com
//
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>

# define WHEEL30         30
# define WHEEL210        210
# define PRIME_PRODUCT   (WHEEL30 * 7 * 11 * 13 * 17 * 19)
# define TEST_FILE       "prime.pi"
# define VERSION         "2.43"
# define CHECK           0

# define L1_SIEVE_SEG    4 //2 - 6
# define L2_SIEVE_SEG    2 //2 - 6
# define L2_DCACHE_SIZE  256

//SSE4.2/SSE4a POPCNT instruction for fast bit counting.
#if _MSC_VER > 1400
	# define POPCNT      0
	# include <intrin.h>
#elif (__GNUC__ * 10 + __GNUC_MINOR__ > 44)
	# define POPCNT      1
//	# include <popcntintrin.h>
#else
	# define POPCNT      0
#endif

#if __x86_64__ || _M_AMD64
	# define TREE2       1
	# define X86_64      1
#else
	# define TREE2       1
#endif

#if defined _M_AMD64
	# define ASM_X86     0
#elif _MSC_VER >= 1200
	# define ASM_X86     0
#else
	# define ASM_X86     0
#endif

//#define CPU 1
#if CPU == 0 //intel core ix
	# define L1_DCACHE_SIZE  32
	# define MAX_SIEVE       1024
	# define SEGS            2
#elif CPU == 1 //amd k8/10
	# define L1_DCACHE_SIZE  64
	# define MAX_SIEVE       512
	# define SEGS            4
#elif CPU == 2 //intel old pentium4
	# define L1_DCACHE_SIZE  32
	# define MAX_SIEVE       512
	# define SEGS            2
#endif

typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;

#ifdef _WIN32
	# pragma warning(disable: 4616 4244 4018 6328 6031)
	typedef __int64 int64;

#if _MSC_VER == 1200
	typedef __int64 uint64;
#else
	typedef unsigned __int64 uint64;
#endif

	#define CONSOLE "CON"
	#include <windows.h>
#else //linux/unix
	typedef long long int64;
	typedef unsigned long long uint64;
	#define CONSOLE "/dev/tty"
	#include <sys/time.h>
	#include <unistd.h>
#endif


#define PADD_DIFF(p, j)   if (p % 2 == 0) { p += -1 + ((uint)Prime[++j] << 8); }
#define NEXT_PRIME(p, j)  p += Prime[++j]
#define MIN(a, b)         (a < b ? a : b)

static const char* const Help = "\
	[B: Benchmark]\n\
	[D: Debug log]\n\
	[G: Progress of calculating]\n\
	[C: Cpu L1/L2 data cache size (L1:16-128, L2:256-1024)]\n\
	[I: Info of sieve]\n\
	[F: Save result to prime.pi]\n\
	[A: Result compare by two algorithm]\n\
	[S: Set sieve segment size (16 - 2048)]\n\
	[U: Unit test prime.pi (cases) (cache) (rw flag)]\n\
	[Y: Find maxp adjacent prime gap [start, end]]\n\
	[P: Print prime in [start, end]]\n\
	[L: List (start) (end/count) (step) (type 0 1 2)]\n\n\
Example:\n\
	c32 c256 s1024 1e16 1e16+1e10\n\
	y 1e19 2^32\n\
	p 1e10 100";

enum EFLAG
{
	PRINT_RET = 1 << 1,
	PRINT_TIME = 1 << 2,
	PRINT_LOG = 1 << 3,
	SAVE_RESUTL = 1 << 4,
	CHCECK_RESUTL = 1 << 5
};

enum BITMASK
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

static struct
{
	uint L1Size;
	uint L1Maxp;
	uint L2Size;
	uint L2Maxp;
}
CpuCache =
{
	32 * (WHEEL30 << 10),
	(32 << 10) / L1_SIEVE_SEG,
	L2_DCACHE_SIZE * (WHEEL30 << 10),
	(L2_DCACHE_SIZE << 10) / L2_SIEVE_SEG,
};

//config
static struct
{
	uint Flag;
	//print calculating time
	uint Progress;
	//sieve size
	uint SieveSize;
}
Config =
{
	PRINT_RET | PRINT_TIME | PRINT_LOG,
	64 - 1, MAX_SIEVE * (WHEEL30 << 10)
};

enum ECMD
{
	COUNT_BITS = 1,
	COPY_BITS,
	SAVE_PRIME,
	SAVE_BYTE,
	SAVE_BYTEGAP,
	FIND_MAXGAP,
	PCALL_BACK
};

struct Cmd
{
	ECMD cmd;
	uint64 Primes;
	uchar* Data;
};

struct WheelPrime
{
	uint Wp; //[0 - 5]: pattern index, [6 - 31]: p / sieve_size
	uint SieveIndex; //[0 - 5]: wheel index, [6 - 31]: sieve offset / WHEEL30
};

struct _BucketInfo
{
	uint CurIndex;
	uint MaxIndex;
	uint Log2Size;
	uint BucketSize;
	uint SieveSize;
};

static _BucketInfo BucketInfo;

struct Stock
{
	WheelPrime* Wheel;
	Stock* NextStock;
};

//each bucket contains multi Stock(by list)
//each stock contains multi wheelprime(by array)
struct Bucket
{
	WheelPrime* CurWheel;
	Stock* StockHead;
	uint WheelSize;
};

//primes < pi(sqrt(1e19 + 1e13))
//pi(2^32) = 203280221,
static const uint SQRT_PRIMES = (203280221 + 1000);
static const uint BUCKET_SIZE = 1 << 13; //depend maxp/sieve_size 10*2^32 / 2^18 * 30
static const uint BLOCK_SIZE = 1 << 10; //12: 32k, best in [10 - 13]
static const uint STOCK_SIZE = SQRT_PRIMES / BLOCK_SIZE + BUCKET_SIZE;

//free stock list
static Stock* StockArray [STOCK_SIZE];
static Bucket BucketArray[BUCKET_SIZE];
//static Bucket* BucketArray = NULL;//[BUCKET_SIZE];
static uint StockSize = 0;
static uint AllStocks = 0;

static const uint MAX_CACHE = (MAX_SIEVE + 2 * L1_DCACHE_SIZE) << 10;
static const uint MAX_PRIME = (280738 / L2_DCACHE_SIZE) * (MAX_CACHE >> 10);
//Prime (i + 1) th - (i)th difference
static uchar Prime[MAX_PRIME];
//each segment sieve_low: (SegIndex[i] + start) % p[j] == 0
static WheelPrime MediumWheel[MAX_PRIME];

//presieved small prime number <=17 bit array.
//the crossing out bit module WHEEL30, the first
//16 bit of PreSieved map to
//----------------------------------------
//|01/1|07/1|11/1|13/1|17/1|19/1|23/1|29/1| = 0x1111 1111 = PreSieved[0]
//----------------------------------------
//|31/1|37/1|41/1|43/1|47/1|49/0|53/1|59/1| = 0x1101 1111 = PreSieved[1]
//----------------------------------------

static uchar PreSieved[PRIME_PRODUCT / WHEEL30];

//position of least significant 1 bit of an integer
static uchar Lsb[1 << 16];

//number of bits 1 binary representation table in Range[0-2^16)
static uchar WordNumBit1[1 << 16];

//cache small segment primes
static int64 PiCache[2001];

//accelerates to save diff prime < 2^31 into Prime[]
//range = [n, n + 2 * WHEEL30]
struct PrimeGap
{
	uchar Gap[13];
	//number of bits1 in wheel range
	uchar Bits;
	//fist bit index 1 in range
	uchar Beg;
	//last bit index 1 in range
	uchar End;
};

PrimeGap* WheelGap = NULL;

struct WheelElement
{
	uchar NextMultiple;
	char WheelIndex;
	uchar UnsetBit;
	uchar Correct;
};

/**
struct WheelInit
{
	char WheelIndex;
	uchar UnsetBit;
	char Correct;
};

struct Skip
{
	uchar NextMultiple;
	char WheelIndex;
}

struct WheelFirst
{
	uchar NextMultiple;
	char WheelIndex;
	uchar UnsetBit;
}
**/

typedef WheelElement WheelInit;
typedef WheelElement WheelFirst;
typedef WheelElement Skip;

static WheelElement Wheel210[48][48];
static WheelFirst WheelFirst210[WHEEL210][48];
static WheelInit WheelInit210[WHEEL210];

//static WheelElement Wheel30[8][8];
static WheelFirst WheelFirst30[WHEEL30][8];
static Skip Skip30[WHEEL30][8];

static const uchar Pattern30[ ] =
{
	1,  7,  11, 13, 17, 19, 23, 29,
	31, 37, 41, 43, 47, 49, 53, 59
};

static const WheelInit WheelInit30[ ] =
{
	{0, -1, 0, 0}, {0,  0, 1, 0}, {0, -1, 0, 1}, {0, -1, 0, 1}, {0, -1, 0, 1},
	{0, -1, 0, 1}, {0, -1, 0, 1}, {0,  1, 2, 1}, {0, -1, 0, 2}, {0, -1, 0, 2},
	{0, -1, 0, 2}, {0,  2, 4, 2}, {0, -1, 0, 3}, {0,  3, 8, 3}, {0, -1, 0, 4},
	{0, -1, 0, 4}, {0, -1, 0, 4}, {0,  4, 16,4}, {0, -1, 0, 5}, {0,  5, 32,5},
	{0, -1, 0, 6}, {0, -1, 0, 6}, {0, -1, 0, 6}, {0,  6, 64,6}, {0, -1, 0, 7},
	{0, -1, 0, 7}, {0, -1, 0, 7}, {0, -1, 0, 7}, {0, -1, 0, 7}, {0,  7, 128,7}
};

//adjacent element difference of pattern,
//MultipleFactor30[i] = Pattern30[j] - Pattern30[j - 1]
static const uchar MultipleFactor30[ ] =
{
	6, 4, 2, 4, 2, 4, 6, 2,
	6, 4, 2, 4, 2, 4, 6, 2,
	6, 4, 2, 4, 2, 4, 6, 2,
	6, 4, 2, 4, 2, 4, 6, 2,
};

typedef void (*call_back)(uint64, uint64);
static int sievePrime(uchar[], uint);
static uint64 initPiCache(uint64, uint64, int, Cmd*);
static int segmentedSieve(uint64 start, int sieve_size, Cmd*);
static void initWheelGap();

//get current time
static double getTime( )
{
#ifdef _WIN32
	LARGE_INTEGER s_freq, performanceCount;
	QueryPerformanceFrequency(&s_freq);
	QueryPerformanceCounter(&performanceCount);
	return 1000. * performanceCount.QuadPart / s_freq.QuadPart;
#else
	struct timeval tmVal;
	gettimeofday(&tmVal, NULL);
	return tmVal.tv_sec * 1000. + tmVal.tv_usec / 1000.;
#endif
}

static int ilog10(uint64 n)
{
	int powbase = 0;
	while (n / 10 > 0) {
		powbase ++;
		n /= 10;
	}

	return powbase;
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

static uint64 ipow(uint64 x, uint n)
{
	uint64 pown = 1;
	while (n != 0) {
		if (n & 1) {
			pown *= x; //overflow !!!
			n -= 1;
		}
		x *= x;
		n /= 2;
	}

	return pown;
}

static uint isqrt(uint64 n)
{
	uint64 rem = 0, root = 0, divisor = 0;

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

//return n % p
static uint fastMod(const uint64 n, uint p)
{
#if defined X86_64 && defined _MSC_VER
	p = (uint)(n % p);
#elif !defined _MSC_VER
	uint loww = n, higw = n >> 32;
	__asm
	(
		"divl %%ecx\n"
		: "=d" (p)
		: "d"(higw), "a"(loww), "c"(p)
	);
#else
	uint loww = n, higw = n >> 32;
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

//the only invalid is [space][number][e][+-*][number]
//e9, 2^32, 1234, 10000*2, 2^30-1E2, 2e9+2^20 all are invalid
uint64 atoint64(const char* str)
{
	uint64 ret = 0;

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

inline static void
crossOffWheelFactor(uchar* ppbeg[], const uchar* pend, const uint p)
{
	uchar* ps0 = ppbeg[0], *ps1 = ppbeg[1];
	uchar* ps2 = ppbeg[2], *ps3 = ppbeg[3];

	while (ps3 <= pend) {
		*ps0 |= BIT0, ps0 += p;
		*ps1 |= BIT1, ps1 += p;
		*ps2 |= BIT2, ps2 += p;
		*ps3 |= BIT3, ps3 += p;
	}
	//not safe
	*ps0 |= BIT0; *ps1 |= BIT1; *ps2 |= BIT2;

	ps0 = ppbeg[4], ps1 = ppbeg[5];
	ps2 = ppbeg[6], ps3 = ppbeg[7];
	while (ps3 <= pend) {
		*ps0 |= BIT4, ps0 += p;
		*ps1 |= BIT5, ps1 += p;
		*ps2 |= BIT6, ps2 += p;
		*ps3 |= BIT7, ps3 += p;
	}
	*ps0 |= BIT4; *ps1 |= BIT5; *ps2 |= BIT6;
}

//30% fast on old cpu pentium 4
static uchar* crossOffWheelFactor2(uchar *pbeg, const uchar *pend, const uint step)
{
	const uint o = step / WHEEL30;
	const uchar windex = WheelInit30[step % WHEEL30].WheelIndex;
	uchar* p = pbeg;

	switch (windex) {
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

#if ASM_X86 == 1
#define ESP_OFFSET4 0x20
_declspec(naked)
#endif
	inline static void
crossOff4Factor(uchar* ppbeg[], const uchar* pend, const uint mask, const uint p)
{
#if ASM_X86 && _MSC_VER
	__asm
	{
		pushad //store all register into stack, esp will decrease 0x20
		//ppbeg  esp + 32 + 04
		//p   esp + 32 + 16
		//define ebx, ebp, esi, edi as ppbeg[0], ppbeg[1], ppbeg[2], ppbeg[3]
		mov eax, dword ptr [esp + 4 + ESP_OFFSET4]
		mov ebx, dword ptr [eax + 0]
		mov ebp, dword ptr [eax + 4]
		mov esi, dword ptr [eax + 8]
		mov edi, dword ptr [eax + 12]

		//define eax, edx, ecx as pend, mask p
		mov eax, dword ptr [esp + 8 + ESP_OFFSET4]
		mov edx, dword ptr [esp + 12 + ESP_OFFSET4]
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
#if 0
		*ps0 |= umask.bmask[0]; ps0 += p;
		*ps1 |= umask.bmask[1]; ps1 += p;
		*ps2 |= umask.bmask[2]; ps2 += p;
		*ps3 |= umask.bmask[3]; ps3 += p;
#else
		*ps0 |= mask >> 0, ps0 += p;
		*ps1 |= mask >> 8, ps1 += p;
		*ps2 |= mask >>16, ps2 += p;
		*ps3 |= mask >>24, ps3 += p;
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
inline static void
crossOff2Factor(uchar* ps0, uchar* ps1, const uchar* pend,
		const ushort smask, const int p)
{
#if ASM_X86 && _MSC_VER
	__asm
	{
		push esi
		push edi

		mov edi, dword ptr [esp + ESP_OFFSET2 + 4]	//ps0
		mov esi, dword ptr [esp + ESP_OFFSET2 + 8]	//ps1
		mov edx, dword ptr [esp + ESP_OFFSET2 + 12]	//pend
		mov eax, dword ptr [esp + ESP_OFFSET2 + 16]	//smask
		mov ecx, dword ptr [esp + ESP_OFFSET2 + 20]	//p

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

	const uchar masks1 = smask >> 8;
	const uchar masks0 = (uchar)smask;

	while (ps1 <= pend) {
		*ps1 |= masks1; ps1 += p;
		*ps0 |= masks0; ps0 += p;
	}
	if (ps0 <= pend)
		*ps0 |= masks0;
#endif
}

//Pre-sieve multiples of small primes <= limit (e.g. 19).
//Resets the sieve array (resets bits to 1) of SieveOfEratosthenes
//objects after each sieved segment and removes the multiples of
//small primes without sieving.
static int preSieve(uchar bitarray[], const uint64 start, const int sieve_size)
{
	assert(start % WHEEL30 == 0);
	const int offset = (int)(start % PRIME_PRODUCT) / WHEEL30;
	const int bits = sieve_size / WHEEL30 * 8 + WheelInit30[sieve_size % WHEEL30].Correct;
	const int bytes = (bits + 7) / 8;

	if (offset + bytes < sizeof(PreSieved)) {
		memcpy(bitarray, PreSieved + offset, bytes);
	} else {
		memcpy(bitarray, PreSieved + offset, sizeof(PreSieved) - offset);
		memcpy(bitarray + sizeof(PreSieved) - offset, PreSieved, bytes + offset - sizeof(PreSieved));
	}

	if (start < WHEEL30) {
		//set bit position 1, 2 for prime 1, 7
		bitarray[0] = 0x3;
	}

	//pack the last 2 dword with bit 1
	uint* pdata = (uint*)bitarray + (bits >> 5);
	pdata[0] |= ~((1u << (bits & 31)) - 1), pdata[1] = ~0, pdata[2] = ~0;

	return bytes;
}

static int sieveBytes(const uint64 start, const int sieve_size)
{
	int bits = sieve_size + (int)(start % WHEEL30);
	bits = (bits / WHEEL30 * 8 + WheelInit30[bits % WHEEL30].Correct);
	return (bits + 7) / 8;
}

//popcnt instruction : INTEL ix/SSE4.2, AMD Phonem/SSE4A
//Use of the SSE 4.2 POPCNT instruction for fast bit counting
#if POPCNT
inline static int countBitsPopcnt(const uint64 n)
{
#ifdef X86_64
	return _mm_popcnt_u64(n);
#else
	return _mm_popcnt_u32(n) + _mm_popcnt_u32((uint)(n >> 32));
#endif
}
#endif

inline static int countBitsTable(const uint64 n)
{
	uint hig = (uint)(n >> 32), low = (uint)n;
	return WordNumBit1[(ushort)hig] + WordNumBit1[(ushort)low] +
		WordNumBit1[hig >> 16] + WordNumBit1[low >> 16];
}

inline static int countBitsTree2(uint64 n)
{
#if 1
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

//count number of bit 0 in binary representation of array
//bit 0 mean it's a prime position
static int countBit0sArray(const uint64 bitarray[], const int bitsize)
{
	int bit1s = 0;
	int loops = bitsize >> 6;

	while (loops-- >= 0) {
#if POPCNT
		bit1s += countBitsPopcnt(*bitarray++);
#elif TREE2
		bit1s += countBitsTree2(*bitarray++);
#else
		bit1s += countBitsTable(*bitarray++);
#endif
	}

	return ((1 + (bitsize >> 6)) << 6) - bit1s;
}

//print prime from bit buffer
static void printPrime(uint64 primes, uint64 prime)
{
	printf("%llu %llu\n", primes, prime);
}

static int callBack(const ushort bitarray[], uint64 sieve_low,
		const int word_size, uint64 sum_prime, call_back func)
{
	int primes = 0;
	sieve_low -= sieve_low % WHEEL30;

	for (int bi = 0; bi < word_size; bi++) {
		ushort mask = ~bitarray[bi];
		while (mask > 0) {
			func(++primes + sum_prime, sieve_low + Lsb[mask]);
			mask &= mask - 1;
		}
		sieve_low += WHEEL30 * 2;
	}

	return primes;
}

//get prime from bit buffer
//word_size > 1476 / 60
static int findPrimeGap(const ushort bitarray[], uint64 sieve_low, const int word_size, uchar *prime)
{
	uint64* result = (uint64*)prime;
	ushort pdiff = (ushort)result[0];
	ushort max_gap = (ushort)result[1];
	int i = 0, skip_words = max_gap / (WHEEL30 * 2) - 1;
	if (skip_words < 1)
		skip_words = 1;

	//find first
	for (i = 0; pdiff == 0; i++) {
		const ushort msk = ~bitarray[i] & 0xffff;
		pdiff = WheelGap[msk].End;//bug for small range
		initWheelGap();
	}

	for (; i < word_size; i++) {
		const ushort msk = ~bitarray[i] & 0xffff;
		if (msk == 0) {
			pdiff += WHEEL30 * 2;
			continue;
		}

		const ushort gap = pdiff + WheelGap[msk].Beg;
		if (gap > max_gap) { //
			max_gap = gap;
			result[2] = sieve_low + i * WHEEL30 * 2 + Lsb[msk];
		}
		if (bitarray[i + skip_words] != 0xffff && i + skip_words < word_size) {
			i += skip_words - 1;
			pdiff = 0;
		} else {
			pdiff = WheelGap[msk].End;
		}
	}

	result[0] = pdiff;
	result[1] = max_gap;
	return 0;
}

//get prime from bit buffer
static int savePrime(const ushort bitarray[], uint64 sieve_low, const int word_size, uint64* prime)
{
	int primes = 0;
	sieve_low -= sieve_low % WHEEL30;

	for (int bi = 0; bi < word_size; bi++) {
		ushort mask = ~bitarray[bi];
		while (mask > 0) {
			prime[primes++] = sieve_low + Lsb[mask];
			mask &= mask - 1;
		}
		sieve_low += WHEEL30 * 2;
	}

	return primes;
}

static int savePrimeByte(const ushort bitarray[], uint sieve_low, const int word_size, uchar* prime)
{
	int primes = 0, lastp = sieve_low + prime[0];
	sieve_low -= sieve_low % WHEEL30;

	for (int bi = 0; bi < word_size; bi++) {
		ushort mask = ~bitarray[bi];
		while (mask > 0) {
			int curp = sieve_low + Lsb[mask];
			prime[primes++] = curp - lastp;
			lastp = curp;
			mask &= mask - 1;
		}
		sieve_low += WHEEL30 * 2;
	}

	prime[primes] = lastp - sieve_low;

	return primes;
}

static int savePrimeGap(ushort bitarray[], const uint64 start, const int word_size, uchar* prime)
{
	ushort pdiff = *(char*)prime;
	int primes = 0;
	initWheelGap();

	for (int i = 0; i < word_size; i++) {
		const uint msk = ~bitarray[i] & 0xffff;
		if (msk == 0) {
			pdiff += WHEEL30 * 2;
			continue;
		}

		ushort gap = pdiff + WheelGap[msk].Beg;
		*prime = gap;
		if (gap > 255) {
			*(ushort*)(prime++) = gap + 1;
		}
#if 1
		*(uint64*)(prime + 1) = *(uint64*)(WheelGap[msk].Gap + 0);
#else
		*(uint*)(prime + 1) = *(uint*)(WheelGap[msk].Gap + 0);
		*(uint*)(prime + 5) = *(uint*)(WheelGap[msk].Gap + 4);
#endif
		pdiff = WheelGap[msk].End;
		uchar bits = WheelGap[msk].Bits;
		if (bits > 8)
			*(uint*)(prime + 9) = *(uint*)(WheelGap[msk].Gap + 8);
		prime += bits;
		primes += bits;
	}

	prime[0] = pdiff;
	prime[1] = 0;

	return primes;
}

static int initBucketInfo(const uint sieve_size, const uint sqrtp, uint64 range)
{
	BucketInfo.CurIndex = 0;
	BucketInfo.MaxIndex = (range / sieve_size) + 1;//overflow for large range

	//maxp wheel 210 pattern, max pattern difference is 10
	BucketInfo.BucketSize = ((uint64)sqrtp * 10) / sieve_size + 2;
	assert (range / sieve_size < -1u);
	assert (BucketInfo.MaxIndex < BUCKET_SIZE || BucketInfo.BucketSize <= BUCKET_SIZE);

	BucketInfo.BucketSize = (2 << ilog2(BucketInfo.BucketSize)) - 1;

	assert (BucketInfo.BucketSize <= BUCKET_SIZE);
	BucketInfo.Log2Size = ilog2(sieve_size / WHEEL30);
	BucketInfo.SieveSize = (1 << BucketInfo.Log2Size) - 1;

	const uint pix = (uint)(sqrtp / log((double)sqrtp) * (1 + 1.200 / log((double)sqrtp)));

	//203280221 >> 12
	StockSize = pix / BLOCK_SIZE + MIN(BucketInfo.BucketSize, BucketInfo.MaxIndex) / 2;
	if (range < sqrtp / 4)
		StockSize >>= 2;
	else if (range < sqrtp)
		StockSize >>= 1;

	assert(StockSize < STOCK_SIZE);

	//dynamic memory use
	WheelPrime *pwheel = (WheelPrime*) malloc(StockSize * BLOCK_SIZE * sizeof(WheelPrime));
	assert(pwheel);

	Stock* pstock = (Stock*) malloc(StockSize * sizeof(Stock));
	for (int j = 0; j < StockSize; j++) {
		StockArray[j] = pstock + j;
	}

	const uint aligned = (uint64)pwheel % (8 * BLOCK_SIZE);
	StockArray[0]->Wheel = pwheel;
	StockArray[1]->Wheel = (WheelPrime*)((uint64)pwheel + (8 * BLOCK_SIZE) - aligned);
	for (int i = 2; i < StockSize; i++) {
		StockArray[i]->Wheel = StockArray[i - 1]->Wheel + BLOCK_SIZE;
	}

	AllStocks = StockSize;

	if (!BucketArray) {
//		BucketArray = (Bucket*)malloc(sizeof(Bucket) * BUCKET_SIZE);
	}

//	printf("new memory %d MB\n", StockSize * BLOCK_SIZE * sizeof(WheelPrime) >> 20);
	return pix;
}

inline static void
pushBucket(const uint sieve_low, const uint wp, const uchar wheel_index)
{
	const uint next_bucket = (sieve_low >> BucketInfo.Log2Size) + BucketInfo.CurIndex;
	if (next_bucket < BucketInfo.MaxIndex) {
		Bucket* pbucket = BucketArray + (next_bucket & BucketInfo.BucketSize);
		if (pbucket->WheelSize++ % BLOCK_SIZE == 0) {
			Stock* pstock = StockArray[-- StockSize];
			pbucket->CurWheel = pstock->Wheel;
			pstock->NextStock = pbucket->StockHead;
			pbucket->StockHead = pstock;
		}

		WheelPrime* nextwheel = pbucket->CurWheel++;
		nextwheel->Wp = wp;
		nextwheel->SieveIndex = (sieve_low & BucketInfo.SieveSize) << 6 | wheel_index;
	}
}

static void initMediumWheel(const uint sieve_size, const uint maxp, const uint64 start, bool check)
{
	uint j = 5, p = 11;
	for (; p < maxp; NEXT_PRIME(p, j)) {
		uint medsieve = CpuCache.L2Size;
		if (check || p >= CpuCache.L2Maxp)
			medsieve = sieve_size;

		const uint64 pp = (uint64)p * p;
		uint sieve_low;
		if (start <= pp) {
			sieve_low = (pp - start) % medsieve;
		} else {
			sieve_low = p - (uint)(start % p);
		}

		const WheelFirst fw = WheelFirst30[sieve_low % WHEEL30][WheelInit30[p % WHEEL30].WheelIndex];
		sieve_low += fw.Correct * p;

		if (pp > start + medsieve) {
			sieve_low %= medsieve;
		}

		MediumWheel[j + 0].Wp = p;
		MediumWheel[j + 0].SieveIndex = (sieve_low << 3) | fw.NextMultiple & 7;
	}

	MediumWheel[j].Wp = -1u;
	assert(j < sizeof(MediumWheel) / sizeof(MediumWheel[0]));
}

static void initBucketStock(const uint sieve_size, const uint sqrtp, const uint64 start, const uint64 range)
{
	initBucketInfo(sieve_size, sqrtp, range);

	//static ushort bitarray[MAX_CACHE / 4];//prime = 252515
	uint pmax = (start >> 32) + 1, nextp = 0;
	uint64 remp = 0;

	for (uint offset = sieve_size / SEGS, segsize = CpuCache.L2Size / 2; offset < sqrtp; offset += segsize) {

		if (segsize > sqrtp - offset)
			segsize = sqrtp - offset;

		//Cmd cmd = {COPY_BITS, 0, (uchar*)bitarray};
		Cmd cmd = {COPY_BITS, 0, 0}; //no copy but no save
		segmentedSieve(offset, segsize, &cmd);

		ushort* wbitarray = (ushort*)cmd.Data;
		ushort mask = 0;
		uint p2 = offset - offset % WHEEL30 - sizeof(mask) * WHEEL30;

		for (int j = 0; j <= 1 + ((int)cmd.Primes) / sizeof(mask); ) {
			if (mask == 0) {
				mask = ~wbitarray[j++];
				p2 += WHEEL30 * sizeof(mask);
				continue;
			}

			const uint p = p2 + Lsb[mask];
			mask &= mask - 1;
//			printf("p = %d\n", p);

			if (p > nextp) {
#if MAX_SIEVE < 1024
				remp = start / (nextp = p + ((uint64)p * p) / pmax);
#else
				remp = start / (nextp = p + 10000);
#endif
				//for large prime overflow
				if (p > nextp) { remp = start / (nextp = -1u); }
			}

			uint sieve_low = p - fastMod(start - remp * p, p);
//			uint sieve_low = p - start % p;
			if (sieve_low >= range) // || (sieve_low % 2 == 0 && p > range))
				continue;

			const uint wp = p / WHEEL210 * 64 + WheelInit210[p % WHEEL210].WheelIndex;
			const uchar module_wheel = sieve_low % WHEEL210;
			const WheelFirst& fw = WheelFirst210[module_wheel][wp & 63];

			sieve_low = sieve_low / WHEEL210 * 7 + (wp >> 6) * fw.Correct * 7;
			sieve_low += (fw.Correct * (p % WHEEL210) + module_wheel) / WHEEL30;
			pushBucket(sieve_low, wp, fw.WheelIndex);
		}
	}

	assert(StockSize > 0);
}

inline static void
sieveSmall0(uchar bitarray[], const uchar* pend, const uint p, uint sieve_low, int skip_index)
{
#if _MSC_VER
	uchar* ppbeg[8];
	for (int i = 0; i < 8; i++) {
		ppbeg[WheelInit30[sieve_low % WHEEL30].WheelIndex] = bitarray + sieve_low / WHEEL30;
		sieve_low += MultipleFactor30[skip_index++] * p;
	}
	crossOffWheelFactor(ppbeg, pend, p);
#else //pentium4 fast
	for (int i = 0; i < 8; i++) {
		const uchar mask = WheelInit30[sieve_low % WHEEL30].UnsetBit;
		bitarray[sieve_low / WHEEL30] |= mask;
		if (mask == 1)
			break;
		sieve_low += MultipleFactor30[skip_index++] * p;
	}
	crossOffWheelFactor2(bitarray + sieve_low / WHEEL30, pend + 1 - 0, p);
#if 0
	sieve_low = 1;
	const int sieve_size = (1 + pend - ps) * WHEEL30;
	while ((int)sieve_low <= sieve_size) {
		ps[sieve_low / WHEEL30] |= WheelInit30[sieve_low % WHEEL30].UnsetBit;
		sieve_low += MultipleFactor30[skip_index++] * p;
	}
#endif
#endif
}

inline static void
sieveSmall1(uchar bitarray[], const uchar* pend, const uint p, uint sieve_low, int skip_index)
{
	for (int i = 0; i < 4; i++) {
		uchar* ps0 = bitarray + sieve_low / WHEEL30;
		ushort smask = WheelInit30[sieve_low % WHEEL30].UnsetBit;
		sieve_low += MultipleFactor30[skip_index++] * p;

		uchar* ps1 = bitarray + sieve_low / WHEEL30;
		smask |= WheelInit30[sieve_low % WHEEL30].UnsetBit << 8;
		sieve_low += MultipleFactor30[skip_index++] * p;
		crossOff2Factor(ps0, ps1, pend, smask, p);
	}
}

inline static void
sieveSmall2(uchar bitarray[], const uchar* pend, const uint p, uint sieve_low, int skip_index)
{
#if 1
	uchar* ppbeg[8], mask[8];
	for (int i = 0; i < 8; i++) {
		ppbeg[i] = bitarray + sieve_low / WHEEL30;
		mask[i] = WheelInit30[sieve_low % WHEEL30].UnsetBit;
		sieve_low += MultipleFactor30[skip_index++] * p;
	}
	crossOff4Factor(ppbeg + 0, pend, *(uint*)(mask + 0), p);
	crossOff4Factor(ppbeg + 4, pend, *(uint*)(mask + 4), p);
#else
	for (int i = 0; i < 8; i++) {
		uchar* pbeg = bitarray + sieve_low / WHEEL30;
		uchar mask = WheelInit30[sieve_low % WHEEL30].UnsetBit;
		sieve_low += MultipleFactor30[skip_index++] * p;
		for (; pbeg < pend; pbeg += p)
			*pbeg |= mask;
	}
#endif
}

inline static void
eratSieveL1(uchar bitarray[], const uint64 start, const int segsize, uint maxp)
{
	const uchar* pend = bitarray + segsize / WHEEL30;

	if ((start + segsize) < ((uint64)maxp) * maxp) {
		maxp = (uint)sqrt((double)start + segsize) + 2;
	}

	for (uint p = Prime[0], j = 8 + PRIME_PRODUCT / 9699690; p < maxp; NEXT_PRIME(p, j)) {
		uint sieve_low = p - (uint)(start % p);
		if (start <= p) {
			sieve_low = p * p - (uint)start;
		}

		const WheelFirst fw = WheelFirst30[sieve_low % WHEEL30][WheelInit30[p % WHEEL30].WheelIndex];
#if _MSC_VER
		sieveSmall1(bitarray, pend, p, sieve_low + fw.Correct * p, fw.NextMultiple);
#else
		sieveSmall0(bitarray, pend, p, sieve_low + fw.Correct * p, fw.NextMultiple);
#endif
	}
}

static void eratSieveSmall(uchar bitarray[], const uint64 start, const uint sieve_size)
{
	const uint maxp = CpuCache.L1Maxp;
	uint segsize = CpuCache.L1Size;

	//preSieve(bitarray, start, sieve_size);
	//#pragma omp parallel for num_threads(2)
	for (uint sieve_low = 0; sieve_low < sieve_size; sieve_low += segsize) {
		if (segsize + sieve_low > sieve_size)
			segsize = sieve_size - sieve_low;
		preSieve(bitarray + sieve_low / WHEEL30, start + sieve_low, segsize);
		eratSieveL1(bitarray + sieve_low / WHEEL30, start + sieve_low, segsize, maxp);
	}
}

static void eratSieveMedium2(uchar bitarray[], const uint64 start,
		const uint sieve_size, const uint minp, uint maxp)
{
	uint j = 8 + PRIME_PRODUCT / 9699690, p = Prime[0];//why int more fast than uint ?
	if (minp >= 8192) {
		j = 1029, p = 8209;
	}

	while (p < minp) {
		NEXT_PRIME(p, j);
	}

	const uint mins = MIN(maxp, sieve_size);
	for (; p < mins; NEXT_PRIME(p, j)) {
		uint sieve_low = p - (uint)(start % p);
		const WheelFirst fw = WheelFirst30[sieve_low % WHEEL30][WheelInit30[p % WHEEL30].WheelIndex];
		sieve_low += fw.Correct * p;
		sieveSmall1(bitarray, bitarray + sieve_size / WHEEL30, p, sieve_low, fw.NextMultiple);
	}

	const uint pmax = (start >> 32) + 1;
	uint64 remp = 0;
	for (uint nextp = 0; p < maxp; NEXT_PRIME(p, j)) {

		PADD_DIFF(p, j);
		if (p > nextp) {
			remp = start / (nextp = p + ((uint64)p * p) / pmax);
			if (p > nextp) remp = start / (nextp = -1u);
			//remp = start / (nextp = p + 2000);
		}
		uint sieve_low = p - fastMod(start - remp * p, p);
		if (sieve_low <= sieve_size) {
		//	if (sieve_low % 2)
			bitarray[sieve_low / WHEEL30] |= WheelInit30[sieve_low % WHEEL30].UnsetBit;
		}
	}
}

//core code of this algorithm for large range
//sieve prime multiples in [start, start + sieve_size)
static void eratSieveMedium(uchar bitarray[], const uint64 start,
		const uint sieve_size, const uint minp, uint maxp)
{
	if ((start + sieve_size) < (uint64)maxp * maxp) {
		maxp = (uint)sqrt((double)start + sieve_size) + 1;
	}

	if (maxp < minp)
		return ;

	uint j = 8 + PRIME_PRODUCT / 9699690, p = Prime[0];//CPU, why int more fast than uint ?
	if (minp >= 131072) {
		j = 12252, p = 131101;
	} else if (minp >= 8192) {
		j = 1029, p = 8209;
	}

	while (p < minp) {
		p = MediumWheel[++j].Wp;
	}

	WheelPrime* pwheel = MediumWheel + j;
	const uint bytes = sieve_size / WHEEL30 + 1;
	const uint pmin = MIN(maxp, bytes / 1);

	for (; p < pmin; p = pwheel->Wp) {
		const uint sieve_low = pwheel->SieveIndex >> 3;
		const uchar skip_index = pwheel->SieveIndex & 7;

		const uint offset = sieve_size - sieve_low;
		const uint rem = offset / p;

		const Skip& skipd = Skip30[rem % WHEEL30][skip_index];
		pwheel++->SieveIndex = ((skipd.NextMultiple + rem) * p - offset) * 8 | skipd.WheelIndex;
		//300/1840
		sieveSmall1(bitarray, bitarray + bytes, p, sieve_low, skip_index);
	}

	const uint maxp2 = MIN(maxp, sieve_size / 6 + 1);
	for (; p < maxp; p = pwheel->Wp) {
		uint sieve_low = pwheel->SieveIndex >> 3;
#if SEGS < 6
		if (sieve_low > sieve_size) {
			pwheel++->SieveIndex -= sieve_size << 3; //15%
			continue;
		}
#endif
		uint skip_index = pwheel->SieveIndex & 7;
		do { //200/1840//
			bitarray[sieve_low / WHEEL30] |= WheelInit30[sieve_low % WHEEL30].UnsetBit;
			sieve_low += MultipleFactor30[skip_index++] * p;
		} while (sieve_low < sieve_size);

		pwheel++->SieveIndex = (sieve_low - sieve_size) << 3 | (skip_index & 7);
	}

#if 0
	for (; p < maxp; p = pwheel->Wp) {
		uint sieve_low = pwheel->SieveIndex >> 3;
#if SEGS < 6
		if (sieve_low > sieve_size) {
			pwheel++->SieveIndex -= sieve_size << 3; //15%
			continue;
		}
#endif
		uint skip_index = pwheel->SieveIndex & 7;
		bitarray[sieve_low / WHEEL30] |= WheelInit30[sieve_low % WHEEL30].UnsetBit;
		sieve_low += MultipleFactor30[skip_index++] * p;
		if (sieve_low < sieve_size) {
			bitarray[sieve_low / WHEEL30] |= WheelInit30[sieve_low % WHEEL30].UnsetBit;
			sieve_low += MultipleFactor30[skip_index++] * p;
		}
		pwheel++->SieveIndex = (sieve_low - sieve_size) << 3 | (skip_index & 7);
	}
#endif
}

//This implementation uses a sieve array with WHEEL30 numbers per byte and
//a modulo 210 wheel that skips multiples of 2, 3, 5 and 7.
static void eratSieveBucket(uchar bitarray[], Bucket* pbucket, uint loops)
{
	uint wheel_prime = pbucket->WheelSize;
	uint sieve_size = BucketInfo.SieveSize;
	pbucket->WheelSize = 0;

	for (; wheel_prime > 0; loops = BLOCK_SIZE) {
		wheel_prime -= loops;

		Stock* phead = StockArray[StockSize++] = pbucket->StockHead;
		pbucket->StockHead = phead->NextStock;
		WheelPrime* pwheel = phead->Wheel;

		while (loops--) {

			uint sieve_low = pwheel->SieveIndex;
			const uint wp = pwheel++->Wp;//int more fast on amd

#if CPU == 1
			uint wheel = *(uint*)&Wheel210[sieve_low & 63][wp & 63];
			bitarray[sieve_low >>= 6] |= wheel >> 16;
			sieve_low += (wheel >> 24) + (wheel & 255) * (wp >> 6);
	#if SEGS > 2
			if (sieve_low <= sieve_size) {
				wheel = *(uint*)&Wheel210[wheel >> 8 & 63][wp & 63];
				bitarray[sieve_low] |= wheel >> 16;
				sieve_low += (wheel >> 24) + (wheel & 255) * (wp >> 6);
			}
	#endif
			pushBucket(sieve_low, wp, wheel >> 8);
#else
			WheelElement& wheel = Wheel210[sieve_low & 63][wp & 63];
			/// 135/1840
			bitarray[sieve_low >>= 6] |= wheel.UnsetBit;
			sieve_low += wheel.Correct;
			sieve_low += wheel.NextMultiple * (wp >> 6);
			pushBucket(sieve_low, wp, wheel.WheelIndex);
#endif
		}
	}

	BucketInfo.CurIndex++;
}

static int doSieveResult(ushort bitarray[], const uint64 start, const int sieve_size, Cmd* cmd)
{
	const uint bytes = sieveBytes(start, sieve_size);
	if (cmd == NULL || cmd->cmd == COUNT_BITS)
		return countBit0sArray((uint64*)bitarray, bytes * 8);

	int primes = 0, words = (1 + bytes) / 2;
	switch (cmd->cmd)
	{
		case COPY_BITS:
			cmd->Data = (uchar*)bitarray;
			primes = bytes; //if (cmd->Data) memcpy(cmd->Data, bitarray, bytes + 8); else
			break;
		case SAVE_BYTE:
			primes = savePrimeByte(bitarray, start, words, cmd->Data + cmd->Primes);
			break;
//		case SAVE_BYTEGAP:
//			primes = savePrimeGap(bitarray, start, words, cmd->Data + cmd->Primes);
//			break;
		case FIND_MAXGAP:
			primes = findPrimeGap(bitarray, start, words, cmd->Data + cmd->Primes);
			break;
		case PCALL_BACK:
			primes = callBack(bitarray, start, words, cmd->Primes, (call_back)cmd->Data);
			break;
		case SAVE_PRIME:
			primes = savePrime(bitarray, start, words, (uint64*)cmd->Data + cmd->Primes);
			break;
	}

	cmd->Primes += primes;

	return primes;
}

//sieve prime multiples in [start, start + sieve_size)
static int segmentedSieve(const uint64 start, const uint sieve_size, const uint wheel_offset, Cmd* cmd = NULL)
{
	//1.sieve small/medium prime factor
	uchar bitarray[MAX_CACHE];
	for (uint sieve_low = 0, segsize = CpuCache.L2Size; sieve_low < sieve_size; sieve_low += segsize) {
		if (segsize + sieve_low > sieve_size)
			segsize = sieve_size - sieve_low;

		//sieve p in [23, L1_DCACHE_SIZE / 30]
		eratSieveSmall(bitarray + sieve_low / WHEEL30, start + sieve_low, segsize);
		//sieve p in [L1_DCACHE_SIZE / 30, L2_DCACHE_SIZE / 30]
		eratSieveMedium(bitarray + sieve_low / WHEEL30, start + sieve_low, segsize, CpuCache.L1Maxp, CpuCache.L2Maxp);
	}

	const uint sqrtp = (uint)sqrt((double)start + sieve_size) + 1;
	const uint minp = MIN(sqrtp, Config.SieveSize / SEGS);
	eratSieveMedium(bitarray, start, sieve_size, CpuCache.L2Maxp, minp);

	//600/1840
	//2.sieve p in [L2_DCACHE_SIZE / 2, sieve_size / 2]
	if (minp != sqrtp /*&& sqrtp - 1 != Config.SieveSize / SEGS */) {
//	if (start + sieve_size > (uint64)(Config.SieveSize / SEGS) * Config.SieveSize / SEGS) {
		//static double time_use = 0; double ts = getTime();
		//3.sieve p [sieve_size / 2, sqrtp]
		Bucket* pbucket = BucketArray + (BucketInfo.CurIndex & BucketInfo.BucketSize);
		int loops = pbucket->WheelSize % BLOCK_SIZE;
		if (loops == 0) loops = BLOCK_SIZE;//why slow moveto loop
		eratSieveBucket(bitarray, pbucket, loops);
		//time_use += getTime() - ts;
		//if (sieve_size != Config.SieveSize) { printf("eratSieveBucket time %.f ms\n", time_use); time_use = 0; }
	}

	if (wheel_offset > 0) {
		memset(bitarray, -1u, wheel_offset / WHEEL30);
		bitarray[wheel_offset / WHEEL30] |= (1 << WheelInit30[wheel_offset % WHEEL30].Correct) - 1;
	}

	return doSieveResult((ushort*)bitarray, start, sieve_size, cmd);
}

static int segmentedSieve(uint64 start, int sieve_size, Cmd* cmd = NULL)
{
	const uint sqrtp = (uint)sqrt((double)start + sieve_size) + 1;
	const uint wheel_offset = start % WHEEL30;
	start -= wheel_offset;
	sieve_size += wheel_offset;

	uchar bitarray[MAX_CACHE / 2];
	eratSieveSmall(bitarray, start, sieve_size);
	eratSieveMedium2(bitarray, start, sieve_size, CpuCache.L1Maxp, sqrtp);

	bitarray[0] |= (1 << WheelInit30[wheel_offset].Correct) - 1;
	return doSieveResult((ushort*)bitarray, start, sieve_size, cmd);
}

static void setCpuSize(uint cdata)
{
	if ((cdata & (cdata - 1)) != 0) {
		cdata = 1 << ilog2(cdata);
	}

	if (cdata >= 16 && cdata < 256) { //L1
		CpuCache.L1Size = cdata * (WHEEL30 << 10);
		CpuCache.L1Maxp = CpuCache.L1Size / (WHEEL30 * L1_SIEVE_SEG);
	} else if (cdata >= 256 && cdata <= 1024) { //L2
		CpuCache.L2Size = cdata * (WHEEL30 << 10);
		CpuCache.L2Maxp = CpuCache.L2Size / (WHEEL30 * L2_SIEVE_SEG);
	}
}

static int setSieveSize(uint sieve_size)
{
	if (sieve_size <= 0) {
		memset(PiCache, 0, sizeof(PiCache));
		return 0;
	}

	if (sieve_size < 12) {
		sieve_size = (WHEEL30 << 10) << sieve_size;
	} if (sieve_size <= 2048) {
		sieve_size *= (WHEEL30 << 10);
	}

	if (sieve_size > MAX_CACHE * WHEEL30)
		sieve_size = MAX_CACHE * WHEEL30 - 8 * CpuCache.L1Size;

//	sieve_size -= sieve_size % (WHEEL210 * 8);
	sieve_size = (1 << ilog2(sieve_size / WHEEL30 + 1)) * WHEEL30;

	if (sieve_size != Config.SieveSize) {
		memset(PiCache, 0, sizeof(PiCache));
	}

	return Config.SieveSize = sieve_size;
}

static int checkSmall(const uint64 start, const uint64 end, uint64 prime[], bool print = false)
{
	int primes = 0;
	const uchar smallp[] = {2, 3, 5, 7};
	for (int i = 0; i < sizeof(smallp) / sizeof(smallp[0]); i++) {
		if (start <= smallp[i] && smallp[i] <= end) {
			primes++;
			if (print)
				printPrime(primes, smallp[i]);
			else if (prime) {
	//			prime[primes + 0] = prime[primes - 1];
	//			prime[primes - 1] = smallp[i];
			}
		}
	}

	return primes;
}

//calculate number of prime in Range[start, end] with start <= end
static uint64 pi2(const uint64 start, const uint64 end, Cmd* cmd = NULL)
{
	assert (start <= end && start >= 0);

	uint64 primes = 0;
	if (Config.Flag & SAVE_RESUTL) {
		freopen("prime.txt", "w", stdout);
		Cmd cmdbuf = {PCALL_BACK, 0, (uchar*)printPrime};
		cmd = &cmdbuf;
	}
	if (start <= 7) {
		if (cmd) {
			primes = checkSmall(start, end, (uint64*)cmd->Data, cmd->cmd == PCALL_BACK);
		} else
			primes = checkSmall(start, end, NULL);
		if (cmd && (cmd->cmd == SAVE_PRIME || cmd->cmd == PCALL_BACK))
			cmd->Primes += primes;
	}

	const uint sieve_size = Config.SieveSize;
	if (end - start <= sieve_size) {
		primes += segmentedSieve(start, (int)(end - start) + 1, cmd);
	} else {
		primes += segmentedSieve(start, sieve_size - (int)(start % sieve_size), cmd);
		primes += initPiCache(start / sieve_size + 2, end / sieve_size + 1, 1, cmd);
		primes += segmentedSieve(end - (int)(end % sieve_size), (int)(end % sieve_size) + 1, cmd);
	}
	if (Config.Flag & SAVE_RESUTL) {
		freopen(CONSOLE, "w", stdout);
	}

	return primes;
}

static uint64 pi(uint64 start, const uint64 end, uint sieve_size, Cmd* cmd)
{
	double ts = getTime();
	double segs = 1000.0 * sieve_size / (end - start);
	uint64 primes = 0;

	int module_start = (int)(start % WHEEL210);
	start -= module_start;

	for (size_t si = 0; start < end; start += sieve_size) {
		if (start + sieve_size > end)
			sieve_size = int(end - start) + (end & 1);

		const int seg = segmentedSieve(start, sieve_size, module_start, cmd);
		primes += seg;

#if CHECK
		const int ok = segmentedSieve(start + module_start, sieve_size - module_start + 1);
		if (ok != seg) {
			printf("start = %I64d, seg = %d != %d = ok\n", start, seg, ok);
		}
#endif
		if ((si++ & Config.Progress) == Config.Progress - 1) {
			double ratio = si * segs;
			printf(">> %02d%%, total time ~ %.2f sec, primes ~= %lld\r",
					((int)ratio) / 10, (getTime() - ts) / ratio, (uint64)(1000 * primes / ratio));
		}
		module_start = 0;
	}

	return primes;
}

static void printResult(const uint64 start, const uint64 end, uint64 primes, double ts)
{
	int sta10 = ilog10(start);
	int end10 = ilog10(end);
	int dif10 = ilog10(end - start + 1);

	if (start > 0) {
		if (start % ipow(10, sta10) == 0 && sta10 > 2)
			printf("PI[%de%d, ", (int)(start / ipow(10, sta10)), sta10);
		else
			printf("PI[%llu, ", (uint64)start);

		if (end % ipow(10, end10) == 0)
			printf("%de%d]", (int) (end / ipow(10, end10)), end10);
		else if ((end - start) % ipow(10, dif10) == 0 && dif10 > 2) {
			if (start % ipow(10, sta10) == 0) {
				printf("%de%d+", (int)(start / ipow(10, sta10)), sta10);
				printf("%de%d]", (int)((end - start) / ipow(10, dif10)), dif10);
			} else
				printf("%llu+%de%d]", start, (int)((end - start) / ipow(10, dif10)), dif10);
		} else {
			printf("%llu]", end);
		}
	} else if (end % ipow(10, end10) == 0 && end10 > 2) {
		printf("PI(%de%d)", (int) (end / ipow(10, end10)), end10);
	} else {
		printf("PI(%llu)", (uint64)end);
	}

	printf(" = %llu", (uint64)primes);
	if (Config.Flag & PRINT_TIME)
		printf(", time use %.2lf sec\t", (getTime() - ts) / 1000.0);
	putchar('\n');
}

uint64 doSievePrime2(const uint64 start, const uint64 end, Cmd* cmd = NULL)
{
	double ts = getTime();
	const uint sqrtp = (uint)sqrt((double)end);

	sievePrime(Prime + 1, sqrtp);

	const uint64 primes = pi2(start, end, cmd);

	if ((!cmd || cmd->cmd != FIND_MAXGAP) && (Config.Flag & PRINT_RET))
		printResult(start, end, primes, ts);

	return primes;
}

uint64 doSievePrime(const uint64 start, const uint64 end, Cmd* cmd = NULL)
{
	assert (start <= end);

	if (Config.Flag & SAVE_RESUTL)
		freopen("prime.txt", "w", stdout);

	double ts = getTime();

	uint64 primes = checkSmall(start, end, NULL);

	uint sieve_size = Config.SieveSize;
	const uint sqrtp = (uint)sqrt((double)end) + 1;
	if (sqrtp > sieve_size / SEGS && sieve_size < CpuCache.L2Size && sqrtp > 1000000) {
		sieve_size = setSieveSize(CpuCache.L2Size);
	}

	const uint64 segment_low = start - (int)(start % WHEEL210);
	const uint minp = MIN(sqrtp, sieve_size / SEGS);

#if CHECK
	sievePrime(Prime + 1, sqrtp);
#else
	sievePrime(Prime + 1, minp);
#endif

	initMediumWheel(sieve_size, minp, segment_low, sieve_size < CpuCache.L2Size);

	bool bucket_sieve = sqrtp > sieve_size / SEGS && sqrtp >= CpuCache.L2Maxp;
	if (bucket_sieve) {
		initBucketStock(sieve_size, sqrtp, segment_low, end - segment_low);
		if (Config.Flag & PRINT_LOG) {
			printf("init bucket time use %.2f sec and sieve size = %d k\n",
					(getTime() - ts) / 1000.0, sieve_size / (WHEEL30 << 10));
		}
		ts = getTime();
	}

	primes += pi(start, end, sieve_size, cmd);
	if (bucket_sieve) {
		assert(AllStocks == StockSize);
		free(StockArray[0]->Wheel);
		free(StockArray[0]);
	}

	if ((!cmd || cmd->cmd != FIND_MAXGAP) && (Config.Flag & PRINT_RET))
		printResult(start, end, primes, ts);

	if (Config.Flag & SAVE_RESUTL)
		freopen(CONSOLE, "w", stdout);

	return primes;
}

static int dumpPrime(double ts, uchar prime[], const char* file)
{
	int primes = 8 + PRIME_PRODUCT / 9699690;
	uint p = Prime[0];
	//assert(prime[1] + prime[2] == 3);
	for (; prime[primes] > 0; primes++) {
		p += prime[primes];
		//PADD_DIFF(p, primes);
	}

	//		assert(Pi2(0, p + 1) == primes);
	printf("\nPrime[%u] = %u, ", primes, p);
	printf("save primes gap use %.2lf sec\n", ts / 1000.0);

	return primes;
}

static int sievePrime(uchar prime[], uint n)
{
	static uint maxp = 100000;
	if (n <= maxp) {
		return 0;
	}
	maxp = n;

	double ts = getTime( );
	int primes = checkSmall(0, 7, NULL);

#if 0
	prime[primes] = 23 - WHEEL30;
	Cmd cmd = {SAVE_BYTEGAP, primes, prime};
#else
	prime[primes] = 11 - 7;
	Cmd cmd = {SAVE_BYTE, primes, prime};
#endif

	const int sieve_size = Config.SieveSize;
	setSieveSize(CpuCache.L2Size);

#if 0
	initMediumWheel(Config.SieveSize, sqrt((double)n) + 10, 0);
	primes += pi(0, n + 4 * WHEEL30, Config.SieveSize, &cmd) - primes;
#else
	primes += pi2(0, (uint64)n + 4 * WHEEL30, &cmd) - primes;
#endif

	setSieveSize(sieve_size); //reset sieve

	prime[primes + 0] = 0;
	prime[primes + 1] = prime[primes + 2] = (uchar)(-1u);

//	if (Config.Flag & PRINT_LOG)
//		dumpPrime(getTime() - ts, Prime + 1, NULL);

	return primes;
}

static uint64 initPiCache(uint64 starti, uint64 endi, int threads, Cmd* cmd)
{
	uint64 pi = 0;
	const uint sieve_size = Config.SieveSize;

	for (uint64 bi = starti; bi < endi; bi += threads) {
		if (bi < sizeof(PiCache) / sizeof(PiCache[0])) {
			if (PiCache[(uint)bi] == 0 || cmd) {
				PiCache[(uint)bi] = segmentedSieve((bi - 1) * sieve_size, sieve_size, cmd);
			}
			pi += PiCache[(uint)bi];
		} else {
			pi += segmentedSieve((bi - 1) * sieve_size, sieve_size, cmd);
		}

		if (Config.Progress && (bi & 511) == 0) {
			printf("%2d%%\r", 100 * bi / endi);
		}
	}

	return pi;
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
		if (0 == (bitarray[p >> 4] & (1 << (p / 2 & 7)))) {

			Prime[++primes] = p - lastp;
			lastp = p;
			if (p > (1 << 16)) {
				continue;
			}

			for (uint j = p * p / 2; j <= maxp / 2; j += p)
				bitarray[j >> 3] |= 1 << (j & 7);
		}
	}

	//pack the last two byte for safety
	Prime[primes + 2] = Prime[primes + 1] = 255;

	return primes;
}

//The first presieved template
//sieve the first 8th prime multiples
static void initPreSieved( )
{
	const uchar smallprimes[ ] = {7, 11, 13, 17, 19, 23, 29};

	for (int i = 0; PRIME_PRODUCT % smallprimes[i] == 0; i++) {
		int start = smallprimes[i], p2 = 2 * smallprimes[i];
		for (; start < PRIME_PRODUCT; start += p2) {
			PreSieved[start / WHEEL30] |= WheelInit30[start % WHEEL30].UnsetBit;
		}
	}
}

static void initBitTable( )
{
	//1. init WordNumBit1 table in 0-2^16
	int i;
	int nbitsize = sizeof(WordNumBit1) / sizeof(WordNumBit1[0]);
#if	POPCNT == 0 && TREE2 == 0
	WordNumBit1[0] = 0;
	for (i = 1; i < nbitsize; i++)
		WordNumBit1[i] = WordNumBit1[i >> 1] + (i & 1);
#endif

	//2. init Left most bit table
	nbitsize = sizeof(Lsb) / sizeof(Lsb[0]);
	for (i = 0; i < nbitsize; i += 2) {
		Lsb[i + 0] = Lsb[i >> 1] + 1;
	}

	for (i = 0; i < nbitsize; i += 2) {
		Lsb[i + 0] = Pattern30[Lsb[i]];
		Lsb[i + 1] = Pattern30[0];
	}
}

static void initWheelGap()
{
	if (WheelGap == NULL)
		WheelGap = (PrimeGap*) (malloc(sizeof(PrimeGap) * (1 << 16) ));
	if (WheelGap[3].Bits)
		return;

	for (uint k = 1; k < (1 << 16) - 1; k++) {
		int pattern = 0, bits = 0;
		for (int j = 0; j < 16; j++) {
			if (k & (1 << j)) {
				if (pattern == 0)
					WheelGap[k].Beg = Pattern30[j];
				else
					WheelGap[k].Gap[bits - 1] = Pattern30[j] - pattern;
				bits++;
				pattern = Pattern30[j];
			}
		}
		WheelGap[k].End = 2 * WHEEL30 - pattern;
		WheelGap[k].Bits = bits;
	}
	//	WheelGap[0].Beg = WHEEL30 * 2; WheelGap[0].Bits = 0;
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

	for (int i = 0; i < WHEEL30; i += 1) {
		for (int pi = 0; pi < 8; pi++) {
			int multiples = 0, sieve_index = i;
			if (i % 2 == 0) {
				multiples += 1;
				sieve_index += Pattern30[pi];
			}
			int wi = WheelInit30[sieve_index % WHEEL30].WheelIndex;
			while (wi < 0) {
				sieve_index += Pattern30[pi] * 2;
				wi = WheelInit30[sieve_index % WHEEL30].WheelIndex;
				multiples += 2;
			}
			WheelElement& fw = WheelFirst30[i][pi];
			fw.NextMultiple = multipleIndex[wi][pi];
			fw.WheelIndex = wi;
			fw.Correct = multiples;
			fw.UnsetBit = 1 << wi;
		}
	}

	for (int si = 0; si < 8; si++) {
		for (int n = 0; n < WHEEL30; n++) {
			int ski = si;
			int multiples = MultipleFactor30[ski++] - n;
			while (multiples <= 0) {
				multiples += MultipleFactor30[ski++];
			}
			Skip30[n][si].WheelIndex = ski & 7;
			Skip30[n][si].NextMultiple = multiples;
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
	const int psize = sizeof(Pattern210) / sizeof(Pattern210[0]);

	for (int j = 0; j < WHEEL210; j += 1) {
		WheelInit210[j].WheelIndex = -1;
		if (WheelInit30[j % WHEEL30].WheelIndex >= 0 && j % 7 != 0)
			WheelInit210[j].WheelIndex = wi++;
		WheelInit210[j].UnsetBit = WheelInit30[j % WHEEL30].UnsetBit;
	}

	for (wi = 0; wi < psize; wi++) {
		for (int pi = 0; pi < psize; pi++) {
			int next = Pattern210[wi] + 2 * Pattern210[pi];
			int multiples = 2;
			while (WheelInit210[next % WHEEL210].WheelIndex < 0) {
				next += Pattern210[pi] * 2;
				multiples += 2;
			}
			WheelElement& wdata = Wheel210[wi][pi];
			wdata.Correct = next / WHEEL30 - Pattern210[wi] / WHEEL30;
			wdata.UnsetBit = WheelInit210[Pattern210[wi]].UnsetBit;
			wdata.WheelIndex = WheelInit210[next % WHEEL210].WheelIndex;
			wdata.NextMultiple = multiples * 7;
		}
	}

	for (i = 0; i < WHEEL210; i += 1) {
		for (int pi = 0; pi < psize; pi++) {
			int multiples = 0, next = i;
			if (i % 2 == 0) {
				multiples += 1;
				next += Pattern210[pi];
			}
			int wi = WheelInit210[next % WHEEL210].WheelIndex;
			while (wi < 0) {
				next += Pattern210[pi] * 2;
				wi = WheelInit210[next % WHEEL210].WheelIndex;
				multiples += 2;
			}
			WheelFirst210[i][pi].WheelIndex = wi;
			WheelFirst210[i][pi].Correct = multiples;
		}
	}

#if 0
	uchar MultipleFactor210[psize];
	for (wi = 0; wi < sizeof(MultipleFactor210); wi++) {
		if (wi % psize != 47)
			MultipleFactor210[wi] = Pattern210[(wi + 1) % psize] - Pattern210[wi % psize];
		else
			MultipleFactor210[wi] = WHEEL210 + Pattern210[(wi + 1) % psize] - Pattern210[wi % psize];
	}

	for (wi = 0; wi < psize; wi++) {
		for (int pi = 0; pi < psize; pi++) {
			uchar skip[psize], cwi = wi;
			for (i = 0; i < sizeof(skip); i++) {
				skip[i] = Wheel210[cwi][pi].NextMultiple / 7;
				cwi = Wheel210[cwi][pi].WheelIndex;
			}

			for (i = 0; i < sizeof(skip) / sizeof(skip[0]); i++) {
				if (memcmp(skip, MultipleFactor210 + i, sizeof(skip)) == 0) {
					WheelFirst210[Pattern210[wi]][pi].NextMultiple = i;
					break;
				}
			}
		}
	}
#endif
}

static void startBenchmark( )
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
		4118054813U,// pi(10^11)
		203280221, // pi(2^32)
		155428406, // prime count of [10^12, 10^12 + 2^32]
		143482916, // prime count of [10^13, 10^13 + 2^32]
		133235063, // prime count of [10^14, 10^14 + 2^32]
		124350420, // prime count of [10^15, 10^15 + 2^32]
		116578809, // prime count of [10^16, 10^16 + 2^32]
		109726486, // prime count of [10^17, 10^17 + 2^32]
		103626726, // prime count of [10^18, 10^18 + 2^32]
		98169972,  // prime count of [10^19, 10^19 + 2^32]
		2895317534U// pi[10^15, 10^15+10^11]
	};

	double ts = getTime();

	Config.Flag &= ~(PRINT_LOG | PRINT_RET);

	uint primes = 0;
	setSieveSize(CpuCache.L2Size);
	Config.Progress = 0;
	for (int i = 1; i <= 10; i++) {
		primes = doSievePrime(0, ipow(10, i));
		printf("pi(10^%02d) = %d\n", i, primes);
		assert(primes == primeCounts[i - 1]);
	}
	printf("pi(2^32)  = %d\n", primeCounts[11]);
	assert(primeCounts[11] == doSievePrime(0, ipow(2, 32)));

	Config.Progress = 7;
	setSieveSize(MAX_SIEVE);
	for (int j = 12; j <= 19; j++) {
		uint64 start = ipow(10, j);
		primes = doSievePrime(start, start + ipow(2, 32));
		printf("pi(10^%d, 10^%d+2^32) = %d               \n", j, j, primes);
		assert(primes == primeCounts[j]);
	}

	Config.Progress = 0;
	printf("Time elapsed %.f sec\n", (getTime() - ts) / 1000);
	puts("All Big tests passed SUCCESSFULLY!");

	printf("Sieving the primes within [10^15, 10^15+10^11] randomly\n");

	uint64 lowerBound = ipow(10, 15);
	uint64 upperBound = lowerBound + ipow(10, 11);
	primes = 0;

	while (lowerBound < upperBound) {
		uint sieve_size = ((uint64)rand() * rand() * rand()) % 1000000000 + 1e8;
		if (sieve_size + lowerBound > upperBound)
			sieve_size = upperBound - lowerBound;

		setSieveSize(rand() % 1024);
		primes += doSievePrime(lowerBound, lowerBound + sieve_size - 1);
		lowerBound += sieve_size;
		printf("Remaining chunk: %lld\r", upperBound - lowerBound);
	}
	printf("\npi(10^15+1e11) = %lld\n", primes);
	assert(primes == primeCounts[20]);
	Config.Flag |= (PRINT_LOG | PRINT_RET);
}

//test case code, test data from third party
static int startRandTest(int tcases, int sieve_size, int powbase)
{
	Config.Flag &= ~(PRINT_RET | PRINT_LOG);

	printf("number case = %d, sieve_size = %d, test file %s flag = %d\n",
			tcases, sieve_size, TEST_FILE, powbase);

	if (powbase == 0) {
		if (!freopen(TEST_FILE, "r", stdin)) {
			puts("can not read test data file");
			Config.Flag |= PRINT_RET;
			freopen(CONSOLE, "r", stdin);
			return -1;
		}
	} else {
		Config.Progress = 0;
		if (!freopen(TEST_FILE, "w", stdout)) {
			puts("can not write test data file");
			Config.Flag |= PRINT_RET;
			freopen(CONSOLE, "w", stdout);
			return -2;
		}
	}

	srand(time(NULL));
	if (sieve_size == 0) {
		setSieveSize(rand() % 257);
	} else if (sieve_size < 6024) {
		setSieveSize(sieve_size);
	} else {
		setSieveSize(32);
	}

	printf("sieve_size = %d\n", Config.SieveSize);

	uint64 maxn = (uint64)pow(10.0, powbase);
	if (maxn < 1E6) {
		maxn = 2e9;
	}

	pi2(0, 1000000000);

	double ts = getTime();
	const char* sformat1 = "%u PI[%u, %u] = %u\n";
	const char* sformat2 = "PI[%u, %u] = %u\n";
	const char* sformat3 = "PI(%u) = %u\n";

	if (sizeof(uint64) != sizeof(uint)) {
		sformat1 = "%u PI[%llu, %llu] = %u\n";
		sformat2 = "PI[%llu, %llu] = %u\n";
		sformat3 = "PI(%llu) = %u\n";
	}

	int failedcases = 0;
	for (int i = 1; i <= tcases; i++) {
		if (powbase) {
			uint64 start = 0 + ((uint64)rand()) * rand() * rand() % maxn;
			uint64 end = 1E6 + ((uint)rand() * rand()) % (1000000000);
			if (start > end)
				start ^= (end ^= (start ^= end));

			if (start < 10) {
				i--;
				continue;
			}

			if (sieve_size == 0) {
				printf(sformat3, end, pi2(0, end));
			} else {
				printf(sformat1, i, start, end, pi2(start, end));
			}
		} else {
			char linebuf[420] = {0};
			uint64 index, start = 0, end;
			int retInFile;
			gets(linebuf);
			if (sscanf(linebuf, sformat1, &index, &start, &end, &retInFile) != 4 &&
				sscanf(linebuf, sformat2, &start, &end, &retInFile) != 3 &&
				sscanf(linebuf, sformat3, &end, &retInFile) != 2) {
				printf("case %d with wrong data %s\n", i, linebuf);
				if (failedcases++ > 10)
					break;
			}
#if 0
			int primes = pi2(start, end);
#else
			int primes = doSievePrime(start, end);
#endif

			if (primes != retInFile) {
				printf(sformat1, i, start, end, primes);
				printf("ecprime %llu -n%llu\n", (uint64)end, (uint64)start);
			}
//			printf("case %d\r", i);
			if ((i & 255) == 0)
				printf("case pass %d%%\r", i * 100 / tcases);
		}
	}

	freopen(CONSOLE, "r", stdin);
	freopen(CONSOLE, "w", stdout);
	Config.Flag |= PRINT_RET;
	printf("test case time use %.2lf ms\n", getTime() - ts);

	return 0;
}

//make sure the input value is valid
static void listPrime(char params[][60])
{
	uint64 start = 0, end = 1e9, step = 1e6;
	uint64 buf[8] = {0, start, end, step, 0};
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
	printf("[%llu:%llu:%llu]\n", (uint64)start, (uint64)end, (uint64)step);

	if (Config.Flag & SAVE_RESUTL) {
		Config.Flag &= ~PRINT_TIME;
		freopen(TEST_FILE, "w", stdout);
	} else if (step < 1000) {
		Config.Flag &= ~PRINT_TIME;
	}

	if (start <= end && end - start >= step - 1) {
		int flag = (int)buf[4];
		for (uint64 i = start; i <= end - step + 1 && i >= start; i += step) {
			uint64 tend = i + step - 1;
			if (flag == 0) {
				doSievePrime(i, tend);
			} else if (flag == 1) {
				doSievePrime(start, tend);
			} else {
				doSievePrime(0, tend);
			}
		}
	}

	if (Config.Flag & SAVE_RESUTL) {
		Config.Flag &= ~SAVE_RESUTL;
		Config.Flag |= PRINT_TIME;
		freopen(CONSOLE, "w", stdout);
	}

	printf("time use %.2lf ms\n", getTime() - ts);
}

static void cpuid(int cpuinfo[4], int id)
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

// http://msdn.microsoft.com/en-us/library/hskdteyh%28v=vs.80%29.aspx
static int getCpuInfo()
{
	char cpuName[255] = {-1};
	int (*pTmp)[4] = (int(*)[4])cpuName;
	cpuid(*pTmp++, 0x80000002);
	cpuid(*pTmp++, 0x80000003);
	cpuid(*pTmp++, 0x80000004);

	for (int i = 0; cpuName[i]; i++) {
		if (cpuName[i] != ' ' || cpuName[i + 1] != ' ')
			putchar(cpuName[i]);
	}

	int cpuinfo[4];
	cpuid(cpuinfo, 0x80000006);
	printf(", L2 = %d kb\n\n", cpuinfo[2] >> 16);

	if (cpuName[0] == 'A') { //amd cpu
		CpuCache.L1Size = 64 * (WHEEL30 << 10);
	} else {
		CpuCache.L1Size = 32 * (WHEEL30 << 10);
	}

	CpuCache.L2Size = L2_DCACHE_SIZE * (WHEEL30 << 10);
	CpuCache.L1Maxp = CpuCache.L1Size / (WHEEL30 * L1_SIEVE_SEG);
	CpuCache.L2Maxp = CpuCache.L2Size / (WHEEL30 * L2_SIEVE_SEG);

	return cpuinfo[2] >> 16;
}

static int getSystemInfo( )
{
#ifdef _WIN32
	SYSTEM_INFO si;
	GetSystemInfo(&si);

	if (si.wProcessorArchitecture == PROCESSOR_ARCHITECTURE_INTEL)
		printf("Cpu arch = x86, ");
#if PROCESSOR_ARCHITECTURE_AMD64
	else if (si.wProcessorArchitecture == PROCESSOR_ARCHITECTURE_AMD64)
		printf("Cpu arch = x86-64, ");
#endif
	return si.dwNumberOfProcessors;
#else
	return sysconf(_SC_NPROCESSORS_CONF);
#endif
}

static void printInfo(int argc)
{
	const char* sepator =
		"-------------------------------------------------------------------------";
	puts(sepator);
	printf("Count/Sieve number of primes in (0, 2^64-1E11), version %s\n", VERSION);
	puts("Implemented by the segmented sieve of eratosthenes [wheel = 30/210]");
	puts("Copyright @ by Huang Yuanbing 2011 - 2012 bailuzhou@163.com");

	puts(sepator);
	puts(sepator);

#ifdef _MSC_VER
	printf("Compiled by MS/vc++ %d", _MSC_VER);
#elif __GNUC__
	printf("Compiled by MinGW/g++ %d.%d.%d",
			__GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#endif

#ifdef X86_64
	printf(" x86-64");
#endif

	printf(" %s %s\n", __TIME__, __DATE__);

	if (argc > 0) {
		puts(sepator);
		printf("[MARCO] : ASM_X86 = %d, SEGS = %d\n", ASM_X86, SEGS);
		printf("[MARCO] : L1_DCACHE_SIZE = %d, L2_DCACHE_SIZE = %d\n",
				L1_DCACHE_SIZE, L2_DCACHE_SIZE);
		printf("[MARCO] : MAX_SIEVE = %dk, BLOCK_SIZE = %dk\n",
			MAX_SIEVE, BLOCK_SIZE >> 7);
	}
	puts(sepator);
	puts(sepator);
}

static void doCompile(const char* flag)
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
		"cl /O2 /Oi /Ot /Oy /GT /GL %s %s %s";
#elif defined X86_64
		"g++ -m64 -msse3 -mpopcnt %s -O3 -funroll-loops -s -pipe %s -o %s";
#else
		"g++ -m32 -msse3 -mpopcnt %s -O3 -funroll-loops -s -pipe %s -o %s";
#endif

	char commpile[255] = {0};
	if (flag == NULL)
		flag = "";
	sprintf(commpile, cxxflag, flag, __FILE__, programming);

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
			case 'H' :
				puts(Help);
				break;
			case 'S':
				setSieveSize(cdata);
				break;
			case 'C':
				setCpuSize(cdata);
				break;
			case 'G':
				Config.Progress = (1 << cdata) - 1;
				if (Config.Progress == 0)
					Config.Progress = (1 << 30) - 1;
				break;
			case 'D':
				Config.Flag ^= PRINT_LOG;
				break;
			case 'F':
				Config.Flag ^= SAVE_RESUTL;
				break;
			case 'A':
				Config.Flag ^= CHCECK_RESUTL;
				break;
			case 'I':
				printf("L1 = %dk, L2 = %dk, SieveSize = %dk\n",
						CpuCache.L1Size / WHEEL30 >> 10,
						CpuCache.L2Size / WHEEL30 >> 10,
						Config.SieveSize / WHEEL30 >> 10);
				break;
			case 'M':
				doCompile(params[i + 1]);
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

	for (int i = 0; i < 64; i++) {
		while (isspace(*ccmd) || ',' == *ccmd) {
			ccmd++;
		}
		char* pc = cmdparams[i];
		if (*ccmd == 0 || *ccmd == ';') {
			break;
		}
		char c = *ccmd;
		bool isvalid = false;
		while (isalnum(c) || c == '^' ||
				c == '+' || c == '-' || c == '*' || c == '=') {
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

bool excuteCmd(const char* cmd)
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
		uint64 start = atoint64(params[cmdi]);
		uint64 end = atoint64(params[cmdi + 1]);
		if (!isdigit(cmdc) && cmdc != 'E') {
			start = atoint64(params[cmdi + 1]);
			end = atoint64(params[cmdi + 2]);
		}

		if (end == 0)
			end = start, start = 0;
		else if (end < start)
			end += start;

		if (cmdc == 'B') {
			puts("-------------------start benchmark------------------------");
			if (isdigit(params[cmdi + 1][0])) {
				for (int i = 11; i < 20; i++) {
					uint64 start = ipow(10, i);
					uint64 size = ipow(10, 9);
					doSievePrime(start, start + size);
				}
			} else
				startBenchmark();
		} else if (cmdc == 'U') {
			puts("-------------------start unit test------------------------");
			int sieve_size = 0, powbase = 0;
			int testcase = atoi(params[cmdi + 1]);
			if (isdigit(params[cmdi + 2][0])) {
				sieve_size = (int)atoint64(params[cmdi + 2]);
				if (isdigit(params[cmdi + 3][0]))
					powbase = (int)atoint64(params[cmdi + 3]);
			}
			startRandTest(testcase, sieve_size, powbase);
		} else if (cmdc == 'L') {
			puts("-------------------start list multi result---------------");
			listPrime(params);
		} else if (cmdc == 'P') {
			puts("-------------------start print prime---------------------");
			Cmd cmd = {PCALL_BACK, 0, (uchar*)printPrime};
			doSievePrime(start, end, &cmd);
		} else if (cmdc == 'Y') {
			//y 1425172824437699411 1476
			//http://www.ieeta.pt/~tos/gaps.html
			puts("-------------------start find max gap -------------------");
			uint64 data[4] = {0}; double ts = getTime();
			Cmd cmd = {FIND_MAXGAP, 0, (uchar*)data};
			doSievePrime(start, end, &cmd);
			printf("maxp prime gap = %d on %llu, time use = %.2f sec\n",
					data[1], data[2] - data[1], (getTime() - ts) / 1000);
		} else if (cmdi >= 0 && end > 0) {
			puts("-------------------start count primes -------------------");
			doSievePrime(start, end);
			if (Config.Flag & CHCECK_RESUTL)
				doSievePrime2(start, end);
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
		getSystemInfo();
		getCpuInfo();

		eratoSieve(123456);
		initPreSieved( );
		initBitTable( );
//		initWheelGap( );
		initWheel30( );
		initWheel210( );
		setSieveSize(sieve_size);
		initOnce = false;
	}
}

int main(int argc, char* argv[])
{
	if (argc < 2) {
		printInfo(argc);
//		puts(Help);
	}

	initPrime(MAX_SIEVE);

	if (argc > 1)
		excuteCmd(argv[1]);

//	uint64 start = (uint64)(Config.SieveSize / 2) * (Config.SieveSize / 2);
//	doSievePrime(start - 1e9, start + 2e9);
//	excuteCmd("b 0");
	excuteCmd("1e14, 10^14+1e9");
	excuteCmd("1e18, 10^9");

	while (true) {
		char ccmd[1023];
		printf("\n[command or number] >> ");
		if (!gets(ccmd) || !excuteCmd(ccmd))
			break;
	}

	return 0;
}

/***
bugs:
	doSievePrime(start - 1e9, start + 2e9);

OS: windows 7 32 bit
MINGW: gcc 4.6.3
CPU: Intel core i5 560m 2.66G (L1 32k, L2 256k, L3 3M)
CXXFLAGS: -Ofast -msse4 -s -pipe  -march=corei7 -funroll-loops

range                            primesieve   primenumber     Oliveira(not adding init time)
[1E10, 1E10+1E10] = 427154205    3.65         3.09
[1E11, 1E11+1E10] = 394050419    5.04         3.85
[1E12, 1E12+1E10] = 361840208    6.15         4.92             6.0
[1E13, 1E13+1E10] = 334067230    7.78         6.41             7.4
[1E14, 1E14+1E10] = 310208140    9.68         7.93             8.6
[1E15, 1E15+1E10] = 289531946    10.9         9.42             9.9
[1E16, 1E16+1E10] = 271425366    13.0         11.2/10.8        11.2
[1E17, 1E17+1E10] = 255481287    15.0         13.3/12.3        12.6
[1E18, 1E18+1E10] = 241272176    20.1         16.0/13.6        14.0
[1E19, 1E19+1E10] = 228568014    31.4         21.6/15.0
[1E19, 1E19+1E11] = 2285693139   184          156.

[1E18, 1E18+1E7 ] = 241295       1.64         0.74
[1E18, 1E18+1E8 ] = 2414886      2.31         1.67
[1E18, 1E18+1E9 ] = 24127085     4.78         3.12
[1E19, 1E19+1E9 ] = 22854258     10.7         6.42
[1E18, 1E18+1E11] = 2412731214   168.         138.             140

windows 7 64 bit, AMD Althon 2 X4 641 2.8G  / Intel i3 350M 2.26G
                                 primesieve      primenumber
PI[1E11, 1E11+1E9] = 39475591    0.41/0.61       0.54/0.56
PI[1E12, 1E12+1E9] = 36190991    0.55/0.83       0.68/0.68
PI[1E13, 1E13+1E9] = 33405006    0.72/1.04       0.87/0.84
PI[1E14, 1E14+1E9] = 31019409    0.89/1.28       1.07/1.06
PI[1E15, 1E15+1E9] = 28946421    1.04/1.50       1.28/1.31
PI[1E16, 1E16+1E9] = 27153205    1.20/1.72       1.46/1.56
PI[1E17, 1E17+1E9] = 25549226    1.55/2.06       1.61/1.74
PI[1E18, 1E18+1E9] = 24127085    1.72/2.15       1.70/1.84
PI[1E19, 1E19+1E9] = 22854258    2.02/2.35       1.80/1.92
PI[1E18, 1E18+1E10]= 241272176   21.1/21.5       18.7/20.3
PI[1E18, 1E18+1E11]= 2412731214  228./212        193./208
PI[1E18, 1E18+1E12]= 24127637783 2270            1937/2080

vc++
	cl /O2 /Oi /Ot /Oy /GT /GL PrimeNumber.cpp
mingw/gcc
	g++ -march=native -O3 -funroll-loops -s -pipe PrimeNumber.cpp -o PrimeNumber
	-fprofile-use -flto -fprofile-generate
else
	g++ -O3 -s -pipe PrimeNumber.cpp -o PrimeNumber -fprofile-generate/use -flto
doc:
	http://www.ieeta.pt/~tos/software/prime_sieve.html
	http://code.google.com/p/primesieve/wiki/Links
	Cache optimized linear sieve
***/
