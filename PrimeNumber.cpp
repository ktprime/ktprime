//TODO: make it simple
//
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>

# define WHEEL           30
# define WHEEL210        210
# define PRIME_PRODUCT   (WHEEL * 7 * 11 * 13 * 17 * 19)
# define TEST_FILE       "prime.pi"
# define VERSION         "2.3"

# define L1_SIEVE_SEG    4 //2 - 6
# define L2_SIEVE_SEG    2 //2 - 6

//SSE4.2/SSE4a POPCNT instruction for fast bit counting.
#if _MSC_VER > 1400
	# define POPCNT      1
	# include <intrin.h>
#elif (__GNUC__ * 10 + __GNUC_MINOR__ > 44)
	# define POPCNT      0
	# include <popcntintrin.h>
#else
	# define POPCNT      0
#endif

#if __x86_64__ || _M_AMD64
	# define TREE2       0
	# define X86_64      1
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
# define L2_DCACHE_SIZE   256

static const int MAX_CACHE = (MAX_SIEVE + 2 * L1_DCACHE_SIZE) << 10;

typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;

#ifdef _WIN32
	typedef __int64 int64;
	typedef unsigned __int64 uint64;
	#define CONSOLE "CON"
	#include <windows.h>
#else //linux/unix
	typedef long long int64;
	typedef unsigned long long uint64;
	#define CONSOLE "/dev/tty"
	#include <sys/time.h>
	#include <unistd.h>
#endif

#if _MSC_VER == 1200
	typedef int64 ltype;
#else
	typedef uint64 ltype;
#endif

#define PADD_DIFF(p, j)   if (p % 2 == 0) { p += -1 + ((uint)Prime[++j] << 8); }
#define NEXT_PRIME(p, j)  p += Prime[++j]
#define MIN(a, b)         (a < b ? a : b)

static struct
{
	uint L1Size;
	uint L1Maxp;
	uint L2Size;
	uint L2Maxp;
}
CpuCache =
{
	32 * (WHEEL << 10),
	(32 << 10) / L1_SIEVE_SEG,
	L2_DCACHE_SIZE * (WHEEL << 10),
	(L2_DCACHE_SIZE << 10) / L2_SIEVE_SEG,
};

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

enum FLAG
{
	PRINT_RET = 1 << 1,
	PRINT_TIME = 1 << 2,
	PRINT_LOG = 1 << 3,
	SAVE_RESUTL = 1 << 4,
	CHCECK_RESUTL = 1 << 5
};

//config
static struct
{
	uint Flag;
	//print calculating time
	uint Progress;
	//number of threads
	uint Threads;
	//sieve size
	uint SieveSize;
}
Config =
{
	PRINT_RET | PRINT_TIME | PRINT_LOG,
	64 - 1, 4, MAX_SIEVE * (WHEEL << 10)
};

enum ECMD
{
	COUNT_BITS = 1,
	COPY_BITS,
	SAVE_PRIME,
	SAVE_BYTEGAP,
	FIND_MAXGAP,
	PCALL_BACK
};

struct CmdData
{
	ECMD Cmd;
	ltype Primes;
	uchar* Data;
};

//primes < pi(sqrt(1e19 + 1e13))
//pi(2^32) = 203280221,
static const uint SQRT_PRIMES = (203280221 + 1000);

struct WheelPrime
{
	//sieve index
	uint SieveIndex; //[0 - 5]: wheel index, [6 - 31]: sieve offset / WHEEL
	//sieve prime
	uint Wp; //[0 - 5]: pattern index, [6 - 31]: p / sieve_size
};

//the bucket data
static struct _BucketInfo
{
	uint CurIndex;
	uint MaxIndex;
	uint BucketSize;
	uint Log2Size;
	uint SieveSize;
} BucketInfo;

const uint BUCKET_SIZE = 1 << 13; //depend maxp/sieve_size 10*2^32 / 2^18 * 30
const uint BLOCK_SIZE = 1 << 11; //12: 32k, best in [10 - 13]
const uint STOCK_SIZE = SQRT_PRIMES / BLOCK_SIZE + BUCKET_SIZE;

struct Stock
{
	WheelPrime* Wheel;
	Stock* NextStock;
};

//each bucket contains multi Stock(by list)
//each stock contains multi block(by array)
static struct _Bucket
{
	WheelPrime* CurWheel;
	Stock* StockHead;
	uint WheelSize;
} Bucket[BUCKET_SIZE];

//free stock list
static Stock* StockArray[STOCK_SIZE];
static int StockSize = 0;
static int AllStocks = 0;

const uint MaxPrime = (280738 / L2_DCACHE_SIZE) * (MAX_CACHE >> 10);
//each segment sieve_index: (SegIndex[i] + start) % p[j] == 0
static WheelPrime MediumWheel[MaxPrime];

//Prime (i + 1) th - (i)th difference
static uchar Prime[MaxPrime];

//presieved small prime number <=17 bit array.
//the crossing out bit module WHEEL, the first
//16 bit of PreSieved map to
//----------------------------------------
//|01/1|07/1|11/1|13/1|17/1|19/1|23/1|29/1| = 0x1111 1111 = PreSieved[0]
//----------------------------------------
//|31/1|37/1|41/1|43/1|47/1|49/0|53/1|59/1| = 0x1101 1111 = PreSieved[1]
//----------------------------------------

static uchar PreSieved[PRIME_PRODUCT / WHEEL];

//position of least significant 1 bit of an integer
static uchar Lsb[1 << 16];

//number of bits 1 binary representation table in Range[0-2^16)
static uchar WordNumBit1[1 << 16];

//cache small segment primes
static int64 PiCache[10001];

//accelerates to save diff prime < 2^31 into Prime[]
//range = [n, n + 2 * WHEEL]
static struct _WheelGap
{
	uchar Gap[13];
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

static struct _Skip
{
	uchar SkipIndex;
	uchar Multiple;
} Skip[WHEEL][8];

struct WheelData
{
	char Index;
	uchar Mask;
	short Leng;
};

static WheelData WheelData210[WHEEL210];
static uchar MultipleFactor210[48];

static const uchar Pattern[ ] =
{
	1,  7,  11, 13, 17, 19, 23, 29,
	31, 37, 41, 43, 47, 49, 53, 59
};

static const WheelData WheelData30[ ] =
{
	{-1, 0, 0}, { 0, 1, 0}, {-1, 0, 1}, {-1, 0, 1}, {-1, 0, 1},
	{-1, 0, 1}, {-1, 0, 1}, { 1, 2, 1}, {-1, 0, 2}, {-1, 0, 2},
	{-1, 0, 2}, { 2, 4, 2}, {-1, 0, 3}, { 3, 8, 3}, {-1, 0, 4},
	{-1, 0, 4}, {-1, 0, 4}, { 4, 16,4}, {-1, 0, 5}, { 5, 32,5},
	{-1, 0, 6}, {-1, 0, 6}, {-1, 0, 6}, { 6, 64,6}, {-1, 0, 7},
	{-1, 0, 7}, {-1, 0, 7}, {-1, 0, 7}, {-1, 0, 7}, { 7, 128,7}
};

//adjacent element difference of pattern,
//MultipleFactor30[i] = Pattern[j] - Pattern[j - 1]
static const uchar MultipleFactor30[ ] =
{
	6, 4, 2, 4, 2, 4, 6, 2,
	6, 4, 2, 4, 2, 4, 6, 2,
	6, 4, 2, 4, 2, 4, 6, 2,
};

typedef void (*call_back)(ltype, ltype);
static int sievePrime(uchar[], uint);
static ltype initPiCache(ltype, ltype, int, CmdData*);
static int segmentedSieve(ltype start, int sieve_size, CmdData*);
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

static int ilog10(ltype n)
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
	ltype rem = 0, root = 0, divisor = 0;

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

//return n % p == 0
static uint fastMod(const ltype n, uint p)
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

static inline void
crossOffWheelFactor(uchar* ppbeg[], const uchar* pend, const uint step)
{
	uchar* ps0 = ppbeg[0], *ps1 = ppbeg[1];
	uchar* ps2 = ppbeg[2], *ps3 = ppbeg[3];

	while (ps3 <= pend) {
		*ps0 |= 1 << 0, ps0 += step;
		*ps1 |= 1 << 1, ps1 += step;
		*ps2 |= 1 << 2, ps2 += step;
		*ps3 |= 1 << 3, ps3 += step;
	}
	//not safe
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
static uchar* crossOffWheelFactor2(uchar *pbeg, const uchar *pend, const uint step)
{
	const uint o = step / WHEEL;
	const uchar windex = WheelData30[step % WHEEL].Index;
	uchar* p = pbeg;

	switch (windex) {
		case 0 :
			while (p < pend) {
				p[o *  0 +  0] |= 1 << 0, p[o *  6 +  0] |= 1 << 1;
				p[o * 10 +  0] |= 1 << 2, p[o * 12 +  0] |= 1 << 3;
				p[o * 16 +  0] |= 1 << 4, p[o * 18 +  0] |= 1 << 5;
				p[o * 22 +  0] |= 1 << 6, p[o * 28 +  0] |= 1 << 7;
				p += step;
			}
			break;
		case 1 :
			while (p < pend) {
				p[o *  0 +  0] |= 1 << 0, p[o *  4 +  0] |= 1 << 7;
				p[o *  6 +  1] |= 1 << 3, p[o * 10 +  2] |= 1 << 2;
				p[o * 16 +  3] |= 1 << 6, p[o * 18 +  4] |= 1 << 1;
				p[o * 24 +  5] |= 1 << 5, p[o * 28 +  6] |= 1 << 4;
				p += step;
			}
			break;
		case 2 :
			while (p < pend) {
				p[o *  0 +  0] |= 1 << 0, p[o *  2 +  0] |= 1 << 6;
				p[o *  6 +  2] |= 1 << 1, p[o *  8 +  2] |= 1 << 7;
				p[o * 12 +  4] |= 1 << 3, p[o * 18 +  6] |= 1 << 5;
				p[o * 20 +  7] |= 1 << 2, p[o * 26 +  9] |= 1 << 4;
				p += step;
			}
			break;
		case 3 :
			while (p < pend) {
				p[o *  0 +  0] |= 1 << 0, p[o *  4 +  1] |= 1 << 6;
				p[o *  6 +  2] |= 1 << 5, p[o * 10 +  4] |= 1 << 2;
				p[o * 12 +  5] |= 1 << 1, p[o * 16 +  6] |= 1 << 7;
				p[o * 22 +  9] |= 1 << 4, p[o * 24 + 10] |= 1 << 3;
				p += step;
			}
			break;
		case 4:
			while (p < pend) {
				p[o *  0 +  0] |= 1 << 0, p[o *  6 +  3] |= 1 << 3;
				p[o *  8 +  4] |= 1 << 4, p[o * 14 +  7] |= 1 << 7;
				p[o * 18 + 10] |= 1 << 1, p[o * 20 + 11] |= 1 << 2;
				p[o * 24 + 13] |= 1 << 5, p[o * 26 + 14] |= 1 << 6;
				p += step;
			}
			break;
		case 5 :
			while (p < pend) {
				p[o *  0 +  0] |= 1 << 0, p[o *  4 +  2] |= 1 << 4;
				p[o * 10 +  6] |= 1 << 2, p[o * 12 +  7] |= 1 << 5;
				p[o * 18 + 11] |= 1 << 3, p[o * 22 + 13] |= 1 << 7;
				p[o * 24 + 15] |= 1 << 1, p[o * 28 + 17] |= 1 << 6;
				p += step;
			}
			break;
		case 6 :
			while (p < pend) {
				p[o *  0 +  0] |= 1 << 0, p[o *  2 +  1] |= 1 << 4;
				p[o *  6 +  4] |= 1 << 5, p[o * 12 +  9] |= 1 << 1;
				p[o * 14 + 10] |= 1 << 6, p[o * 20 + 15] |= 1 << 2;
				p[o * 24 + 18] |= 1 << 3, p[o * 26 + 19] |= 1 << 7;
				p += step;
			}
			break;
		case 7 :
			while (p < pend) {
				p[o *  0 +  0] |= 1 << 0, p[o *  2 +  1] |= 1 << 7;
				p[o *  8 +  7] |= 1 << 6, p[o * 12 + 11] |= 1 << 5;
				p[o * 14 + 13] |= 1 << 4, p[o * 18 + 17] |= 1 << 3;
				p[o * 20 + 19] |= 1 << 2, p[o * 24 + 23] |= 1 << 1;
				p += step;
			}
			break;
	}

/**	uint sieve_index = 1;
	while (p <= pend + step) {
		p[sieve_index / WHEEL] |= WheelData30[sieve_index % WHEEL].Mask;
		sieve_index += MultipleFactor30[skip_index++] * p;
	}
*/
	return p;
}

#if ASM_X86 == 1
	#define ESP_OFFSET4 0x20
_declspec(naked)
#endif
static inline void
crossOff4Factor(uchar* ppbeg[], const uchar* pend, const uint mask, const uint step)
{
#if ASM_X86 && _MSC_VER
	__asm
	{
		pushad //store all register into stack, esp will decrease 0x20
		//ppbeg  esp + 32 + 04
		//step   esp + 32 + 16
		//define ebx, ebp, esi, edi as ppbeg[0], ppbeg[1], ppbeg[2], ppbeg[3]
		mov eax, dword ptr [esp + 4 + ESP_OFFSET4]
		mov ebx, dword ptr [eax + 0]
		mov ebp, dword ptr [eax + 4]
		mov esi, dword ptr [eax + 8]
		mov edi, dword ptr [eax + 12]

		//define eax, edx, ecx as pend, mask step
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
		*ps0 |= umask.bmask[0]; ps0 += step;
		*ps1 |= umask.bmask[1]; ps1 += step;
		*ps2 |= umask.bmask[2]; ps2 += step;
		*ps3 |= umask.bmask[3]; ps3 += step;
#else
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
crossOff2Factor(uchar* ps0, uchar* ps1, const uchar* pend,
		const ushort smask, const int step)
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

	const uchar masks1 = smask >> 8;
	const uchar masks0 = (uchar)smask;

	while (ps1 <= pend) {
		*ps1 |= masks1; ps1 += step;
		*ps0 |= masks0; ps0 += step;
	}
	if (ps0 <= pend)
		*ps0 |= masks0;
#endif
}

//Pre-sieve multiples of small primes <= limit (e.g. 19).
//Resets the sieve array (resets bits to 1) of SieveOfEratosthenes
//objects after each sieved segment and removes the multiples of
//small primes without sieving.
static int preSieve(uchar bitarray[], const ltype start, const int sieve_size)
{
	assert(start % WHEEL == 0);
	const int offset = (int)(start % PRIME_PRODUCT) / WHEEL;
	const int bits = sieve_size / WHEEL * 8 + WheelData30[sieve_size % WHEEL].Leng;
	const int bytes = (bits + 7) / 8;

	if (offset + bytes < sizeof(PreSieved)) {
		memcpy(bitarray, PreSieved + offset, bytes);
	} else {
		memcpy(bitarray, PreSieved + offset, sizeof(PreSieved) - offset);
		memcpy(bitarray + sizeof(PreSieved) - offset,
				PreSieved, bytes + offset - sizeof(PreSieved));
	}

	if (start < WHEEL) {
		//set bit position 1, 2 for prime 3, 5
		bitarray[0] = 0x3;
	}

	//pack the last 2 dword with bit 1
	*((uint*)bitarray + (bits >> 5) + 0) |= ~((1u << (bits & 31)) - 1);
	*((uint*)bitarray + (bits >> 5) + 1) = ~0;
	*((uint*)bitarray + (bits >> 5) + 2) = ~0;

	return bytes;
}

static int sieveBytes(const ltype start, const int sieve_size)
{
	int bits = sieve_size + (int)(start % WHEEL);
	bits = (bits / WHEEL * 8 + WheelData30[bits % WHEEL].Leng);
	return (bits + 7) / 8;
}

//popcnt instruction : INTEL ix/SSE4.2, AMD Phonem/SSE4A
//Use of the SSE 4.2 POPCNT instruction for fast bit counting
#if POPCNT
static inline int countBitsPopcnt(const uint64 n)
{
#if X86_64
	return _mm_popcnt_u64(n);
#else
	return _mm_popcnt_u32(n) + _mm_popcnt_u32((uint)(n >> 32));
#endif
}
#endif

static inline int countBitsTable(const uint64 n)
{
	uint hig = (uint)(n >> 32), low = (uint)n;
	return WordNumBit1[(ushort)hig] + WordNumBit1[(ushort)low] +
		WordNumBit1[hig >> 16] + WordNumBit1[low >> 16];
}

static inline int countBitsTree2(uint64 n)
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

//get prime from bit buffer
static int savePrime(const ushort bitarray[], ltype sieve_index, const int word_size, ltype* prime)
{
	int primes = 0;
	sieve_index -= sieve_index % WHEEL;

	for (int bi = 0; bi <= word_size; bi++) {
		ushort mask = ~bitarray[bi];
		while (mask > 0) {
			prime[primes++] = sieve_index + Lsb[mask];
			mask &= mask - 1;
		}
		sieve_index += WHEEL * 2;
	}

	return primes;
}

//print prime from bit buffer
static void printPrime(ltype primes, ltype prime)
{
	printf("%llu %llu\n", primes, prime);
}

static int callBack(const ushort bitarray[], ltype sieve_index,
		const int word_size, ltype sum_prime, call_back func)
{
	int primes = 0;
	sieve_index -= sieve_index % WHEEL;

	for (int bi = 0; bi <= word_size; bi++) {
		ushort mask = ~bitarray[bi];
		while (mask > 0) {
			func(++primes + sum_prime, sieve_index + Lsb[mask]);
			mask &= mask - 1;
		}
		sieve_index += WHEEL * 2;
	}

	return primes;
}

//get prime from bit buffer
//word_size > 1476 / 60
static int findPrimeGap(const ushort bitarray[], ltype sieve_index, const int word_size, ltype *result)
{
	ushort pdiff = (ushort)result[0];
	ushort max_gap = (ushort)result[1];
	int i = 0, skip_words = max_gap / (WHEEL * 2) - 1;
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
			pdiff += WHEEL * 2;
			continue;
		}

		const ushort gap = pdiff + WheelGap[msk].Beg;
		if (gap > max_gap) { //
			max_gap = gap;
			result[2] = sieve_index + i * WHEEL * 2 + Lsb[msk];
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

static int savePrimeGap(ushort bitarray[], const int word_size, uchar* & prime)
{
	ushort pdiff = *(char*)prime;
	int primes = 0;
	initWheelGap();

	for (int i = 0; i < word_size; i++) {
		const uint msk = ~bitarray[i] & 0xffff;
		if (msk == 0) {
			pdiff += WHEEL * 2;
			continue;
		}

		ushort gap = pdiff + WheelGap[msk].Beg;
		*prime = gap;
		if (gap > 255) {
			*(ushort*)(prime++) = gap + 1;
		}
#if 1
		*(ltype*)(prime + 1) = *(ltype*)(WheelGap[msk].Gap + 0);
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

static int initBucketInfo(const uint sieve_size, const uint sqrtp, ltype range)
{
	BucketInfo.CurIndex = 0;
	BucketInfo.MaxIndex = (uint)(range / sieve_size) + 1;//BUGS
	assert(range / sieve_size < -1u);

	//maxp wheel 210 pattern 10
	BucketInfo.BucketSize = ((ltype)sqrtp * 10) / sieve_size + 2;
	if (BucketInfo.MaxIndex >= BUCKET_SIZE)
		assert(BucketInfo.BucketSize <= BUCKET_SIZE);

	BucketInfo.BucketSize = (2 << ilog2(BucketInfo.BucketSize)) - 1;

	assert(BucketInfo.BucketSize <= BUCKET_SIZE);
	BucketInfo.Log2Size = ilog2(sieve_size / WHEEL);
	BucketInfo.SieveSize = (1 << BucketInfo.Log2Size) - 1;

	const uint pix = (uint)(sqrtp / log((double)sqrtp) * (1 + 1.200 / log((double)sqrtp)));

	//203280221 >> 12
	StockSize = pix / BLOCK_SIZE + MIN(BucketInfo.BucketSize, BucketInfo.MaxIndex) / 2;
	if (range < sqrtp / 4)
		StockSize >>= 2;
	else if (range < sqrtp)
		StockSize >>= 1;

	assert(StockSize < sizeof(StockArray) / sizeof(StockArray[0]));

	//dynamic memory use
	WheelPrime *pwheel =
		(WheelPrime*) malloc(StockSize * BLOCK_SIZE * sizeof(WheelPrime));
	assert(pwheel);

	Stock* pstock = (Stock*) malloc(StockSize * sizeof(Stock));
	for (int j = 0; j < StockSize; j++) {
		StockArray[j] = pstock + j;
	}

	const uint aligned = (ltype)pwheel % (8 * BLOCK_SIZE);
	StockArray[0]->Wheel = pwheel;
	StockArray[1]->Wheel = (WheelPrime*)((ltype)pwheel + (8 * BLOCK_SIZE) - aligned);
	for (int i = 2; i < StockSize; i++) {
		StockArray[i]->Wheel = StockArray[i - 1]->Wheel + BLOCK_SIZE;
	}

	AllStocks = StockSize;
//	printf("new memory %d MB\n", StockSize * BLOCK_SIZE * sizeof(WheelPrime) >> 20);
	return pix;
}

static void inline
pushBucket(const uint sieve_index, const uint wp, const uchar wheel_index)
{
	const uint next_bucket = (sieve_index >> BucketInfo.Log2Size) + BucketInfo.CurIndex;
	if (next_bucket < BucketInfo.MaxIndex) {
		_Bucket* pbucket = &Bucket[next_bucket & BucketInfo.BucketSize];
		if (pbucket->WheelSize++ % BLOCK_SIZE == 0) {
			Stock* pstock = StockArray[-- StockSize];
			pbucket->CurWheel = pstock->Wheel;
			pstock->NextStock = pbucket->StockHead;
			pbucket->StockHead = pstock;
		}

		WheelPrime* nextwheel = pbucket->CurWheel++;
		nextwheel->SieveIndex = (sieve_index & BucketInfo.SieveSize) << 6 | wheel_index;
		nextwheel->Wp = wp;
	}
}

static void initMediumWheel(const uint sieve_size, const uint maxp, const ltype start)
{
	sievePrime(Prime + 1, maxp);

	uint j = 5, p = 11;
	for (; p < maxp; NEXT_PRIME(p, j)) {
		uint medsieve = sieve_size;
		if (p < CpuCache.L2Maxp && sieve_size > CpuCache.L2Size)
			medsieve = CpuCache.L2Size;

		const ltype pp = (ltype)p * p;
		uint sieve_index;
		if (start <= pp) {
			sieve_index = (pp - start) % medsieve;
		} else {
			sieve_index = p - (uint)(start % p);
		}

		const FirstWheel cwn = FirstWheel30[sieve_index % WHEEL][WheelData30[p % WHEEL].Index];
		sieve_index += cwn.Correct * p;

		if (pp > start + medsieve) {
			sieve_index %= medsieve;
		}

		MediumWheel[j + 0].SieveIndex = (sieve_index << 3) | cwn.MultipleIndex & 7;
		MediumWheel[j + 0].Wp = p;
	}
	MediumWheel[j].Wp = -1u;
	assert(j < sizeof(MediumWheel) / sizeof(MediumWheel[0]));
}

static void initBucketStock(const uint sieve_size, const uint sqrtp, const ltype start, const ltype range)
{
	initBucketInfo(sieve_size, sqrtp, range);

	//static ushort bitarray[MAX_CACHE / 4];//prime = 252515
	uint pmax = (start >> 32) + 1, nextp = 0;
	ltype remp = 0;

	for (uint offset = sieve_size / SEGS, segsize = CpuCache.L2Size / 4; offset < sqrtp; offset += segsize) {

		if (segsize > sqrtp - offset)
			segsize = sqrtp - offset;

		//CmdData cmd = {COPY_BITS, 0, (uchar*)bitarray};
		CmdData cmd = {COPY_BITS, 0, 0}; //no copy but no save
		segmentedSieve(offset, segsize, &cmd);

		ushort* wbitarray = (ushort*)cmd.Data;
		ushort mask = 0;
		uint p2 = offset - offset % WHEEL - sizeof(mask) * WHEEL;

		for (int j = 0; j <= 1 + ((int)cmd.Primes) / sizeof(mask); ) {
			if (mask == 0) {
				mask = ~wbitarray[j++];
				p2 += WHEEL * sizeof(mask);
				continue;
			}

			const uint p = p2 + Lsb[mask];
			mask &= mask - 1;

			if (p > nextp) {
				remp = start / (nextp = p + ((ltype)p * p) / pmax);
//				remp = start / (nextp = p + 10000);
				//for large prime overflow
				if (p > nextp) { remp = start / (nextp = -1u); }
			}

			uint sieve_index = p - fastMod(start - remp * p, p);
//			uint sieve_index = p - start % p;
			if (sieve_index >= range) // || (sieve_index % 2 == 0 && p > range))
				continue;

			const uint wp = p / WHEEL210 * 64 + WheelData210[p % WHEEL210].Index;
			const uchar module_wheel = sieve_index % WHEEL210;
			const FirstWheel& cwn = FirstWheel210[module_wheel][wp & 63];

			sieve_index = sieve_index / WHEEL210 * 7 + (wp >> 6) * cwn.Correct * 7;
			sieve_index += (cwn.Correct * (p % WHEEL210) + module_wheel) / WHEEL;
			pushBucket(sieve_index, wp, cwn.WheelIndex);
		}
	}

	assert(StockSize > 0);
}

static inline void
sieveSmall0(uchar bitarray[], const uchar* pend, const uint p, uint sieve_index, int skip_index)
{
#if 0
	uchar* ppbeg[8];
	for (int i = 0; i < 8; i++) {
		ppbeg[WheelData30[sieve_index % WHEEL].Index] = bitarray + sieve_index / WHEEL;
		sieve_index += MultipleFactor30[skip_index++] * p;
	}
	crossOffWheelFactor(ppbeg, pend, p);
#else //p4 fast
	for (int i = 0; i < 8; i++) {
		const uchar mask = WheelData30[sieve_index % WHEEL].Mask;
		bitarray[sieve_index / WHEEL] |= mask;
		if (mask == 1)
			break;
		sieve_index += MultipleFactor30[skip_index++] * p;
	}
	uchar* ps = crossOffWheelFactor2(bitarray + sieve_index / WHEEL, pend + 1 - 0, p);
#if 0
	sieve_index = 1;
	const int sieve_size = (1 + pend - ps) * WHEEL;
	while ((int)sieve_index <= sieve_size) {
		ps[sieve_index / WHEEL] |= WheelData30[sieve_index % WHEEL].Mask;
		sieve_index += MultipleFactor30[skip_index++] * p;
	}
#endif
#endif
}

static inline void
sieveSmall1(uchar bitarray[], const uchar* pend, const uint p, uint sieve_index, int skip_index)
{
	for (int i = 0; i < 4; i++) {
		uchar* ps0 = bitarray + sieve_index / WHEEL;
		ushort smask = WheelData30[sieve_index % WHEEL].Mask;
		sieve_index += MultipleFactor30[skip_index++] * p;

		uchar* ps1 = bitarray + sieve_index / WHEEL;
		smask |= WheelData30[sieve_index % WHEEL].Mask << 8;
		sieve_index += MultipleFactor30[skip_index++] * p;
		crossOff2Factor(ps0, ps1, pend, smask, p);
	}
}

static inline void
sieveSmall2(uchar bitarray[], const uchar* pend, const uint p, uint sieve_index, int skip_index)
{
#if 1
	uchar* ppbeg[8], mask[8];
	for (int i = 0; i < 8; i++) {
		ppbeg[i] = bitarray + sieve_index / WHEEL;
		mask[i] = WheelData30[sieve_index % WHEEL].Mask;
		sieve_index += MultipleFactor30[skip_index++] * p;
	}
	crossOff4Factor(ppbeg + 0, pend, *(uint*)(mask + 0), p);
	crossOff4Factor(ppbeg + 4, pend, *(uint*)(mask + 4), p);
#else
	for (int i = 0; i < 8; i++) {
		uchar* pbeg = bitarray + sieve_index / WHEEL;
		uchar mask = WheelData30[sieve_index % WHEEL].Mask;
		sieve_index += MultipleFactor30[skip_index++] * p;
		for (; pbeg < pend; pbeg += p)
			*pbeg |= mask;
	}
#endif
}

static inline void
eratSieveSmall(uchar bitarray[], const ltype start, const int segsize, uint maxp)
{
	const uchar* pend = bitarray + segsize / WHEEL;

	if ((start + segsize) < ((ltype)maxp) * maxp) {
		maxp = (uint)sqrt((double)start + segsize) + 2;
	}

	for (uint p = Prime[0], j = 8 + PRIME_PRODUCT / 9699690; p < maxp; NEXT_PRIME(p, j)) {
		uint sieve_index = p - (uint)(start % p);
		if (start <= p) {
			sieve_index = p * p - (uint)start;
		}

		const FirstWheel cwn = FirstWheel30[sieve_index % WHEEL][WheelData30[p % WHEEL].Index];
#if 1
		sieveSmall0(bitarray, pend, p, sieve_index + cwn.Correct * p, cwn.MultipleIndex);
#else
		sieveSmall1(bitarray, pend, p, sieve_index + cwn.Correct * p, cwn.MultipleIndex);
#endif
	}
}

static void eratSieveSmall2(uchar bitarray[], const ltype start, const uint sieve_size)
{
	const uint maxp = CpuCache.L1Maxp;
	uint segsize = CpuCache.L1Size;

	//preSieve(bitarray, start, sieve_size);
	//#pragma omp parallel for num_threads(2)
	for (uint sieve_index = 0; sieve_index < sieve_size; sieve_index += segsize) {
		if (segsize + sieve_index > sieve_size)
			segsize = sieve_size - sieve_index;
		preSieve(bitarray + sieve_index / WHEEL, start + sieve_index, segsize);
		eratSieveSmall(bitarray + sieve_index / WHEEL, start + sieve_index, segsize, maxp);
	}
}

static void eratSieveMedium2(uchar bitarray[], const ltype start,
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
		uint sieve_index = p - (uint)(start % p);
		const FirstWheel cwn = FirstWheel30[sieve_index % WHEEL][WheelData30[p % WHEEL].Index];
		sieve_index += cwn.Correct * p;
		sieveSmall1(bitarray, bitarray + sieve_size / WHEEL, p, sieve_index, cwn.MultipleIndex);
	}

	const uint pmax = (start >> 32) + 1;
	ltype remp = 0;
	for (uint nextp = 0; p < maxp; NEXT_PRIME(p, j)) {

		PADD_DIFF(p, j);
		if (p > nextp) {
			remp = start / (nextp = p + ((ltype)p * p) / pmax);
			if (p > nextp) remp = start / (nextp = -1u);
			//remp = start / (nextp = p + 2000);
		}
		uint sieve_index = p - fastMod(start - remp * p, p);
		if (sieve_index <= sieve_size) {
		//	if (sieve_index % 2)
			bitarray[sieve_index / WHEEL] |= WheelData30[sieve_index % WHEEL].Mask;
		}
	}
}

//core code of this algorithm for large range
//sieve prime multiples in [start, start + sieve_size)
static void eratSieveMedium(uchar bitarray[], const ltype start,
		const uint sieve_size, const uint minp, uint maxp)
{
	if ((start + sieve_size) < (ltype)maxp * maxp) {
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
	const uint bytes = sieve_size / WHEEL + 1;
	const uint pmin = MIN(maxp, bytes / 2);
	const uchar* pend = bitarray + bytes;

	for (; p < pmin; p = pwheel->Wp) {
		const uint sieve_index = pwheel->SieveIndex >> 3;
		const uchar skip_index = pwheel->SieveIndex & 7;
		const uint offset = sieve_size - sieve_index;
		const uint rem = offset / p;

		_Skip skipd = Skip[rem % WHEEL][skip_index];
		pwheel++->SieveIndex = ((skipd.Multiple + rem) * p - offset) * 8 | skipd.SkipIndex;
		//300/1840
		sieveSmall1(bitarray, pend, p, sieve_index, skip_index);
	}

	for (; p < maxp; p = pwheel->Wp) {
		uint sieve_index = pwheel->SieveIndex >> 3;
#if SEGS < 6
		if (sieve_index > sieve_size) {
			pwheel++->SieveIndex -= sieve_size << 3; //15%
			continue;
		}
#endif

		uint skip_index = pwheel->SieveIndex & 7;
		//200/1840//
		while (sieve_index < sieve_size) {
			bitarray[sieve_index / WHEEL] |= WheelData30[sieve_index % WHEEL].Mask;
			sieve_index += MultipleFactor30[skip_index++] * p;
		}

		pwheel++->SieveIndex = (sieve_index - sieve_size) << 3 | (skip_index & 7);
	}
}

//This implementation uses a sieve array with WHEEL numbers per byte and
//a modulo 210 wheel that skips multiples of 2, 3, 5 and 7.
static void eratSieveBucket(uchar bitarray[], const int cur_index, const uint sieve_size)
{
	_Bucket* pbucket = &Bucket[cur_index];
	uint wheel_prime = pbucket->WheelSize;
	pbucket->WheelSize = 0;

	int loops = wheel_prime % BLOCK_SIZE;
	for (; wheel_prime > 0; loops = BLOCK_SIZE) {
		Stock* phead = StockArray[StockSize++] = pbucket->StockHead;
		pbucket->StockHead = phead->NextStock;
		WheelPrime* pwheel = phead->Wheel;

		if (loops == 0) loops = BLOCK_SIZE;//why slow moveto loop
		wheel_prime -= loops;

		while (loops--) {

			uint sieve_index = pwheel->SieveIndex;
			const uint wp = pwheel++->Wp;//int more fast on amd

#if CPU == 1
			uint wheel = *(uint*)&NextWheel210[sieve_index & 63][wp & 63];
			bitarray[sieve_index >>= 6] |= wheel >> 24;
			sieve_index += (uchar)(wheel >> 16) + (wheel & 255) * (wp >> 6);
	#if SEGS > 2
			if (sieve_index <= sieve_size) {
				wheel = *(uint*)&NextWheel210[(wheel) >> 8 & 63][wp & 63];
				bitarray[sieve_index] |= wheel >> 24;
				sieve_index += (uchar)(wheel >> 16) + (wheel & 255) * (wp >> 6);
			}
	#endif
			pushBucket(sieve_index, wp, wheel >> 8);
#else
			const WheelFactorization& wheel = NextWheel210[sieve_index & 63][wp & 63];
			/// 135/184
			bitarray[sieve_index >>= 6] |= wheel.UnsetBit;
			sieve_index += wheel.Correct + wheel.MultipleIndex * (wp >> 6);
			pushBucket(sieve_index, wp, wheel.WheelIndex);
#endif
		}
	}
}

static int doSieveResult(ushort bitarray[], const ltype start, const int bytes, CmdData* cmd)
{
	int primes = 0;
	if (cmd == NULL || cmd->Cmd == COUNT_BITS) {
		primes = countBit0sArray((uint64*)bitarray, bytes * 8);
	} else if (cmd->Cmd == SAVE_BYTEGAP) {
		primes = savePrimeGap(bitarray, (bytes + 1) / 2, cmd->Data);
	} else if (cmd->Cmd == SAVE_PRIME) {
		primes = savePrime(bitarray, start, bytes / 2, (ltype*)cmd->Data + cmd->Primes);
	} else if (cmd->Cmd == COPY_BITS) {
		primes = bytes;
		if (cmd->Data) memcpy(cmd->Data, bitarray, bytes + 8);
		else cmd->Data = (uchar*)bitarray;
	} else if (cmd->Cmd == FIND_MAXGAP) {
		primes = findPrimeGap(bitarray, start, (bytes + 1) / 2, (ltype*)cmd->Data);
	} else if (cmd->Cmd == PCALL_BACK) {
		primes = callBack(bitarray, start, bytes / 2, cmd->Primes, (call_back)cmd->Data);
	}

	if (cmd)
		cmd->Primes += primes;
	return primes;
}

//core code of this algorithm
//sieve prime multiples in [start, start + sieve_size)
static int segmentedSieve(const ltype start, const uint sieve_size, const uint wheel_offset, CmdData* cmd = NULL)
{
	//1.sieve small/medium prime factor
	uchar bitarray[MAX_CACHE];
	for (uint sieve_index = 0, segsize = CpuCache.L2Size; sieve_index < sieve_size; sieve_index += segsize) {
		if (segsize + sieve_index > sieve_size)
			segsize = sieve_size - sieve_index;

		eratSieveSmall2(bitarray + sieve_index / WHEEL, start + sieve_index, segsize);
		eratSieveMedium(bitarray + sieve_index / WHEEL, start + sieve_index, segsize, CpuCache.L1Maxp, CpuCache.L2Maxp);
	}

	const uint sqrtp = (uint)sqrt((double)start + sieve_size) + 1;
	const uint minp = MIN(sqrtp, Config.SieveSize / SEGS);
	//600/1840
	eratSieveMedium(bitarray, start, sieve_size, CpuCache.L2Maxp, minp);

	if (minp != sqrtp) {
		//static double time_use = 0; double ts = getTime();
		eratSieveBucket(bitarray, BucketInfo.CurIndex & BucketInfo.BucketSize, BucketInfo.SieveSize);
		BucketInfo.CurIndex++;
		//time_use += getTime() - ts;
		//if (sieve_size != Config.SieveSize) { printf("eratSieveBucket time %.f ms\n", time_use); time_use = 0; }
	}

	if (wheel_offset > 0) {
		memset(bitarray, -1u, wheel_offset / WHEEL);
		bitarray[wheel_offset / WHEEL] |= (1 << WheelData30[wheel_offset % WHEEL].Leng) - 1;
	}

	const uint bytes = sieveBytes(start, sieve_size);
	return doSieveResult((ushort*)bitarray, start, bytes, cmd);
}

static int segmentedSieve(ltype start, int sieve_size, CmdData* cmd = NULL)
{
	const uint sqrtp = (uint)sqrt((double)start + sieve_size) + 1;
	const uint wheel_offset = start % WHEEL;
	start -= wheel_offset;
	sieve_size += wheel_offset;

	uchar bitarray[MAX_CACHE];
	eratSieveSmall2(bitarray, start, sieve_size);
	eratSieveMedium2(bitarray, start, sieve_size, CpuCache.L1Maxp, sqrtp);

	bitarray[0] |= (1 << WheelData30[wheel_offset].Leng) - 1;

	const uint bytes = sieveBytes(start, sieve_size);
	return doSieveResult((ushort*)bitarray, start, bytes, cmd);
}

static void setCpuSize(uint cdata)
{
	if ((cdata & (cdata - 1)) != 0) {
		cdata = 1 << ilog2(cdata);
	}

	if (cdata >= 16 && cdata < 256) { //L1
		CpuCache.L1Size = cdata * (WHEEL << 10);
		CpuCache.L1Maxp = CpuCache.L1Size / (WHEEL * L1_SIEVE_SEG);
	} else if (cdata >= 256 && cdata <= 1024) { //L2
		CpuCache.L2Size = cdata * (WHEEL << 10);
		CpuCache.L2Maxp = CpuCache.L2Size / (WHEEL * L2_SIEVE_SEG);
	}
}

static int setSieveSize(int sieve_size)
{
	if (sieve_size <= 0) {
		memset(PiCache, 0, sizeof(PiCache));
		return 0;
	}

	if (sieve_size < 12) {
		sieve_size = (WHEEL << 10) << sieve_size;
	} if (sieve_size < 2048) {
		sieve_size *= (WHEEL << 10);
	}

	if (sieve_size > MAX_CACHE * WHEEL)
		sieve_size = MAX_CACHE * WHEEL - 8 * CpuCache.L1Size;

//	sieve_size -= sieve_size % (WHEEL210 * 8);
	sieve_size = (1 << ilog2(sieve_size / WHEEL + 1)) * WHEEL;

	if (sieve_size != Config.SieveSize) {
		memset(PiCache, 0, sizeof(PiCache));
	}

	return Config.SieveSize = sieve_size;
}

static int checkSmall(const ltype start, const ltype end, ltype prime[], bool print = false)
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
static ltype Pi2(const ltype start, const ltype end, CmdData* cmd = NULL)
{
	assert (start <= end && start >= 0);

	ltype primes = 0;
	if (Config.Flag & SAVE_RESUTL) {
		freopen("prime.txt", "w", stdout);
		CmdData cmdbuf = {PCALL_BACK, 0, (uchar*)printPrime};
		cmd = &cmdbuf;
	}
	if (start <= 7) {
		if (cmd) {
			primes = checkSmall(start, end, (ltype*)cmd->Data, cmd->Cmd == PCALL_BACK);
		} else
			primes = checkSmall(start, end, NULL);
		if (cmd && (cmd->Cmd == SAVE_PRIME || cmd->Cmd == PCALL_BACK))
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

static ltype Pi(ltype start, const ltype end, uint sieve_size, CmdData* cmd)
{
	double ts = getTime();

	int module_start = (int)(start % WHEEL210);
	start -= module_start;
	ltype primes = 0;
	double segs = (end - start) * 1.0 / sieve_size;

	//overflow ci
	for (uint ci = 0, segsize = sieve_size; start < end; start += segsize) {

		if (start + segsize > end)
			segsize = int(end - start) + (end & 1);

		const int seg2 = segmentedSieve(start, segsize, module_start, cmd);
		primes += seg2;

#if _DEBUG && 0
		const int seg = segmentedSieve(start + module_start, segsize - module_start + 1);
		if (seg != seg2) {
			printf("start = %I64d, seg2 = %d != %d = seg\n", start, seg, seg2);
		}
#endif
		if ((ci++ & Config.Progress) == Config.Progress - 1) {
			double ratio = ci * 1000.0 / segs;
			printf("%02d%%, sieve time ~ %.2f sec, ret ~= %lld\r",
					((int)ratio) / 10, (getTime() - ts) / ratio, (ltype)(1000 * primes / ratio));
		}
		module_start = 0;
	}

	return primes;
}

static void printPiResult(const ltype start, const ltype end, ltype primes, double ts)
{
	int sta10 = ilog10(start);
	int end10 = ilog10(end);
	int dif10 = ilog10(end - start + 1);

	if (start > 0) {
		if (start % ipow(10, sta10) == 0 && sta10 > 2)
			printf("PI[%de%d, ", (int)(start / ipow(10, sta10)), sta10);
		else
			printf("PI[%llu, ", (ltype)start);

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
		printf("PI(%llu)", (ltype)end);
	}

	printf(" = %llu", (ltype)primes);
	if (Config.Flag & PRINT_TIME)
		printf(", time use %.2f sec\t", (getTime() - ts) / 1000.0);
	putchar('\n');
}

ltype doSievePrime2(const ltype start, const ltype end, CmdData* cmd = NULL)
{
	double ts = getTime();
	const uint sqrtp = (uint)sqrt((double)end);

	sievePrime(Prime + 1, sqrtp);

	const ltype primes = Pi2(start, end, cmd);

	if (Config.Flag & PRINT_RET)
		printPiResult(start, end, primes, ts);

	return primes;
}

ltype doSievePrime(const ltype start, const ltype end, CmdData* cmd = NULL)
{
	assert (start <= end && end - start < (ltype)Config.SieveSize << 32);

	if (Config.Flag & SAVE_RESUTL)
		freopen("prime.txt", "w", stdout);

	double ts = getTime();

	ltype primes = checkSmall(start, end, NULL);

	uint sieve_size = Config.SieveSize;
	const uint sqrtp = (uint)sqrt((double)end) + 1;
	if (sqrtp > sieve_size / SEGS) {
		if (sieve_size < CpuCache.L2Size && sqrtp > 1000000)
			sieve_size = setSieveSize(CpuCache.L2Size);
	}

	const ltype segment_low = start - (int)(start % WHEEL210);
	const uint minp = MIN(sqrtp, sieve_size / SEGS);
	initMediumWheel(sieve_size, minp, segment_low);

	bool bucket_sieve = sqrtp > sieve_size / SEGS && sqrtp >= CpuCache.L2Maxp;
	if (bucket_sieve) {
		initBucketStock(sieve_size, sqrtp, segment_low, end - segment_low);
		if (Config.Flag & PRINT_LOG) {
			printf("init bucket time use %.2f sec and sieve size = %d k\n",
					(getTime() - ts) / 1000.0, sieve_size / (WHEEL << 10));
		}
		ts = getTime();
	}

	primes += Pi(start, end, sieve_size, cmd);

	if ((!cmd || cmd->Cmd != FIND_MAXGAP) && (Config.Flag & PRINT_RET))
		printPiResult(start, end, primes, ts);

	if (bucket_sieve) {
		assert(AllStocks == StockSize);
		free(StockArray[0]->Wheel);
		free(StockArray[0]);
	}

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
	printf("save primes gap use %.2f sec\n", ts / 1000.0);

	return primes;
}

static int sievePrime(uchar prime[], uint n)
{
	static uint maxp = 10000;
	if (n <= maxp) {
		return 0;
	}
	maxp = n;

	double ts = getTime( );
	int primes = checkSmall(0, 7, NULL);
	prime[primes] = 23 - WHEEL;
	CmdData cmd = {SAVE_BYTEGAP, primes, prime + primes};

	const int sieve_size = Config.SieveSize;
	Config.SieveSize = CpuCache.L2Size;

#if 0
	initMediumWheel(Config.SieveSize, sqrt((double)n) + 10, 0);
	primes += Pi(0, n + 4 * WHEEL, Config.SieveSize, &cmd) - primes;
#else
	primes += Pi2(0, (ltype)n + 4 * WHEEL, &cmd) - primes;
#endif

	setSieveSize(sieve_size);

	prime[primes + 0] = 0;
	prime[primes + 1] = prime[primes + 2] = (uchar)(-1u);

//	if (Config.Printlog)
//		dumpPrime(getTime() - ts, Prime + 1, NULL);

	return primes;
}

static ltype initPiCache(ltype starti, ltype endi, int threads, CmdData* cmd)
{
	ltype pi = 0;
	const uint sieve_size = Config.SieveSize;

	for (ltype bi = starti; bi < endi; bi += threads) {
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
			PreSieved[start / WHEEL] |= WheelData30[start % WHEEL].Mask;
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
		Lsb[i + 0] = Pattern[Lsb[i]];
		Lsb[i + 1] = Pattern[0];
	}
}

static void initWheelGap()
{
	if (WheelGap[3].Bits)
		return;

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
			int multiples = 0, sieve_index = i;
			if (i % 2 == 0) {
				multiples += 1;
				sieve_index += Pattern[pi];
			}
			int wi = WheelData30[sieve_index % WHEEL].Index;
			while (wi < 0) {
				sieve_index += Pattern[pi] * 2;
				wi = WheelData30[sieve_index % WHEEL].Index;
				multiples += 2;
			}
			FirstWheel first_wheel = {multipleIndex[wi][pi], wi, multiples, 1 << wi};
			FirstWheel30[i][pi] = first_wheel;
		}
	}

	for (int si = 0; si < 8; si++) {
		for (int mulp = 0; mulp < WHEEL; mulp++) {
			int ski = si;
			int sum = MultipleFactor30[ski++] - mulp;
			while (sum <= 0) {
				sum += MultipleFactor30[ski++];
			}
			Skip[mulp][si].SkipIndex = ski & 7;
			Skip[mulp][si].Multiple = sum;
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

#if 0
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
#endif
}

void initPrime(int sieve_size)
{
	static bool initOnce = true;
	if (initOnce) {
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
		ltype start = ipow(10, j);
		primes = doSievePrime(start, start + ipow(2, 32));
		printf("pi(10^%d, 10^%d+2^32) = %d               \n", j, j, primes);
		assert(primes == primeCounts[j]);
	}

	Config.Progress = 0;
	printf("Time elapsed %.f sec\n", (getTime() - ts) / 1000);
	puts("All Big tests passed SUCCESSFULLY!");

	printf("Sieving the primes within [10^15, 10^15+10^11] randomly\n");

	ltype lowerBound = ipow(10, 15);
	ltype upperBound = lowerBound + ipow(10, 11);
	primes = 0;

	while (lowerBound < upperBound) {
		uint sieve_size = ((ltype)rand() * rand() * rand()) % 1000000000 + 1e8;
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

	ltype maxn = (ltype)pow(10.0, powbase);
	if (maxn < 1E6) {
		maxn = 2e9;
	}

	Pi2(0, 1000000000);

	double ts = getTime();
	const char* sformat1 = "%u PI[%u, %u] = %u\n";
	const char* sformat2 = "PI[%u, %u] = %u\n";
	const char* sformat3 = "Pi(%u) = %u\n";

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

			if (sieve_size == 0) {
				printf(sformat3, end, Pi2(0, end));
			} else {
				printf(sformat1, i, start, end, Pi2(start, end));
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
#if 0
			int primes = Pi2(start, end);
#else
			int primes = doSievePrime(start, end);
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
	Config.Flag |= PRINT_RET;
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

	if (Config.Flag & SAVE_RESUTL) {
		Config.Flag &= ~PRINT_TIME;
		freopen(TEST_FILE, "w", stdout);
	} else if (step < 1000) {
		Config.Flag &= ~PRINT_TIME;
	}

	if (start <= end && end - start >= step - 1) {
		int flag = (int)buf[4];
		for (ltype i = start; i <= end - step + 1 && i >= start; i += step) {
			ltype tend = i + step - 1;
			if (flag == 0) {
				doSievePrime2(i, tend);
			} else if (flag == 1) {
				doSievePrime2(start, tend);
			} else {
				doSievePrime2(0, tend);
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
	printf(", L2 cache = %d kb\n", cpuinfo[2] >> 16);

	//amd cpu
	if (cpuName[0] == 'A') {
		CpuCache.L1Size = 64 * (WHEEL << 10);
		CpuCache.L2Size = 4 * CpuCache.L1Size;
	} else {
		CpuCache.L1Size = 32 * (WHEEL << 10);
		CpuCache.L2Size = 8 * CpuCache.L1Size;
	}
	CpuCache.L1Maxp = CpuCache.L1Size / (WHEEL * L1_SIEVE_SEG);
	CpuCache.L2Maxp = CpuCache.L2Size / (WHEEL * L2_SIEVE_SEG);

	return cpuinfo[2] >> 16 ;
}

static int getSystemInfo( )
{
#ifdef _WIN32
	SYSTEM_INFO si;
	GetSystemInfo(&si);

	Config.Threads = si.dwNumberOfProcessors;

	if (si.wProcessorArchitecture == PROCESSOR_ARCHITECTURE_INTEL)
		printf("Cpu arch = x86, ");
#if PROCESSOR_ARCHITECTURE_AMD64
	else if (si.wProcessorArchitecture == PROCESSOR_ARCHITECTURE_AMD64)
		printf("Cpu arch = x86-64, ");
#endif
	return Config.Threads;
#else
	return sysconf(_SC_NPROCESSORS_CONF);
#endif
}

static void printInfo(int argc)
{
	const char* sepator =
		"--------------------------------------------------------------------";
	puts(sepator);
	printf("Count/Sieve number of primes in (0, 2^64-1E11), version %s\n", VERSION);
	puts("Implemented by the segmented sieve of eratosthenes [wheel = 30/210]");
	puts("Copyright @ by Huang Yuanbing 2011 - 2012 bailuzhou@163.com");

	puts(sepator);
	puts(sepator);

#ifdef _MSC_VER
	printf("Compiled by MS/vc++ %d", _MSC_VER);
#elif __GNUC__
	printf("Compiled by Mingw/g++ %d.%d.%d",
			__GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#endif

#if X86_64
	printf(" on x86-64 bit");
#endif

	printf(" on %s %s\n", __TIME__, __DATE__);

	Config.Threads = getSystemInfo();

	getCpuInfo();

	if (argc > 0) {
		puts(sepator);
		printf("[MARCO] : ASM_X86 = %d, LIANGBCH = %d\n", ASM_X86, 0);
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
						CpuCache.L1Size / WHEEL >> 10,
						CpuCache.L2Size / WHEEL >> 10,
						Config.SieveSize / WHEEL >> 10);
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

		if (cmdc == 'B') {
			puts("-------------------start benchmark------------------------");
			if (isdigit(params[cmdi + 1][0])) {
				for (int i = 11; i < 20; i++) {
					ltype start = ipow(10, i);
					ltype size = ipow(10, 9);
					doSievePrime(start, start + size);
				}
			}
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
			CmdData cmd = {PCALL_BACK, 0, (uchar*)printPrime};
			doSievePrime(start, end, &cmd);
		} else if (cmdc == 'Y') {
			//y 1425172824437699411 1476
			//http://www.ieeta.pt/~tos/gaps.html
			puts("-------------------start find max gap -------------------");
			ltype data[4] = {0}; double ts = getTime();
			CmdData cmd = {FIND_MAXGAP, 0, (uchar*)data};
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

void static testCode()
{
	printf("\tswitch (windex) {\n");
	for (int i = 0; i < 8; i++) {
		int w = 1;
		int p = Pattern[i];
		printf("\t\tcase %d :\n", Pattern[i]);
		printf("\t\t\twhile (p < pend) {\n");
		for (int j = 0; j < 30; j += 2) {
			if (WheelData30[w % WHEEL].Index >= 0) {
				printf("\t\t\t\tp[o * %2d + %2d] |= 1 << %d;\n", j, w / WHEEL, WheelData30[w % WHEEL].Index);
			}
			w += 2 * p;
		}
		printf("\t\t\t\tp += step;\n");
		printf("\t\t\t}\n");
		printf("\t\t\tbreak;\n");
	}
}


int main(int argc, char* argv[])
{
	if (argc < 2) {
		printInfo(argc);
		puts(Help);
	}

	initPrime(MAX_SIEVE);

	if (argc > 1)
		excuteCmd(argv[1]);

//	testCode();

//	excuteCmd("1e16 1e10 s9");
//	excuteCmd("1e19 1e10 s8 ");

	excuteCmd("1e18 1e9");
//	excuteCmd("1e15 1e8; e10 a ");

	while (true) {
		char ccmd[1023];
		printf("\n[command or number] : ");
		if (!gets(ccmd) || !excuteCmd(ccmd))
			break;
	}

	return 0;
}

/***
bugs:

OS: windows 7 32 bit
MINGW: gcc 4.6.3
CPU: Intel core i5 560m 2.66G (L1 32k, L2 256k, L3 3M)
CXXFLAGS: -Ofast -msse4 -s -pipe  -march=corei7 -funroll-loops

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
[1E19, 1E19+1E10] = 228568014    31.4         31.7       21.6/15.0
[1E19, 1E19+1E11] = 2285693139   184          197        163./156.

[1E18, 1E18+1E7 ] = 241295       1.64         2.91       0.74
[1E18, 1E18+1E8 ] = 2414886      2.31         3.10       1.67
[1E18, 1E18+1E9 ] = 24127085     4.78         4.71       3.12
[1E19, 1E19+1E9 ] = 22854258     10.7         12.5       6.42
[1E18, 1E18+1E11] = 2412731214   168.         169        138.             140

windows 7 64 bit, AMD Althon 2 X4 641 2.8G  / Intel i3 350M 2.26G
                                 primesieve      ktprime
PI[1E11, 1E11+1E9] = 39475591    0.41/0.61       0.54/0.56
PI[1E12, 1E12+1E9] = 36190991    0.55/0.85       0.68/0.68
PI[1E13, 1E13+1E9] = 33405006    0.72/1.06       0.87/0.84
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
mingw
	g++ -mpopcnt -march=native -O3 -funroll-loops -s -pipe PrimeNumber.cpp -o PrimeNumber
gc++
if (cpu spupport popcnt instruction and gcc version > 4.3)
	g++ -mpopcnt -march=native -O3 -s -pipe -lpthread PrimeNumber.cpp -o PrimeNumber
else
	g++ -O3 -s -pipe -lpthread PrimeNumber.cpp -o PrimeNumber -fprofile-generate/use -flto
doc:
	http://code.google.com/p/primesieve/
	http://code.google.com/p/primesieve/wiki/Links
***/
