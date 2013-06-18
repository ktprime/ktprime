//Copyright @ by Huang Yuanbing 2011 - 2013 bailuzhou AT 163.com
//
#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <memory.h>
#include <stdio.h>
#include <time.h>
#include <assert.h>

# define WHEEL30         30
# define WHEEL210        210
# define PRIME_PRODUCT   (210 * 11 * 13 * 17 * 19)
# define TEST_FILE       "prime.pi"
# define KVERSION        "3.0"
# define WHEEL_SKIP      0x26424246
# define L2_DCACHE_SIZE  256

//performance marco
# define MPART           8 //8 - 10
# define WPART           6 //6 - 7
# define ERATBIG         2 //2 - 6
# define ERATSMALL       4 //2 - 6
# define ERATMEDIUM      4 //2 - 6

#if 0
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

//for test
# define CHECK            0
# define RELEASE          0
//x86 cpu only
# define BIT_SCANF        1
# define POPCNT           0
# define TREE2            0

//SSE4.2/SSE4a POPCNT instruction for fast bit counting.
#if _MSC_VER > 1400
	# include <intrin.h>
#elif (__GNUC__ * 10 + __GNUC_MINOR__ > 44)
//	# include <popcntintrin.h>
#endif

#if __x86_64__ || _M_AMD64
	# define X86_64      1
#endif

#if X86_64 && _MSC_VER > 1500
	# define ASM_X86     0 //can not on x64 ms vc
#else
	# define ASM_X86     1
#endif

#define  INTEL  0
#define  AMD    1

//#define CPU 1
#if CPU == INTEL //Intel core ix
	# define L1_DCACHE_SIZE  32
	# define MAX_SIEVE       1024
#elif CPU == AMD
	# define L1_DCACHE_SIZE  64
	# define MAX_SIEVE       512
#else  //intel old pentium4
	# define L1_DCACHE_SIZE  32
	# define MAX_SIEVE       512
#endif

#ifdef _WIN32
//	# pragma warning(disable: 4616 4996 4244 4018 6328 6031)
	#if _MSC_VER == 1200
		typedef __int64 uint64;
	#else
		typedef unsigned __int64 uint64;
	#endif

	#define CONSOLE "CON"
	#include <windows.h>
#else //linux/unix
	typedef unsigned long long uint64;
	#define CONSOLE "/dev/tty"
	#include <sys/time.h>
	#include <unistd.h>
#endif
typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;

#define MIN(a, b)         (a < b ? a : b)
#define PRIME_DIFF(p, j)  p += Prime[++j]
#define PRIME_ADD_DIFF(p)    if (p % 2 == 0) { p += 255; }

#if BIT_SCANF == 0
	#define PRIME_OFFSET(mask)   Lsb[mask]
	typedef ushort stype;
#elif X86_64
	typedef uint64 stype;
	#define PRIME_OFFSET(mask)   Pattern30[bitScanForward(mask)]
#else
	typedef uint stype;
	#define PRIME_OFFSET(mask)   Pattern30[bitScanForward(mask)]
#endif

static const char* const Help = "\
	[B: Benchmark]\n\
	[D: Debug log]\n\
	[M: Progress of calculating]\n\
	[C: Cpu L1/L2 data cache size (L1:16-128, L2:256-1024)]\n\
	[I: Info of sieve]\n\
	[F: Save result to prime.pi]\n\
	[A: Result compare by two algorithm]\n\
	[S: Set sieve segment size (16 - 2048)]\n\
	[U: Unit test prime.pi (cases) (cache) (rw flag)]\n\
	[G: Find maxp adjacent prime gap [start, end]]\n\
	[P: Print prime in [start, end]]\n\
	[L: List (start) (end/count) (step) (type 0 1 2)]\n\n\
Example:\n\
	1e16 1e16+1e10\n\
	g 1e19 2^32\n\
	p 1e10 100";

enum EFLAG
{
	PRINT_RET = 1 << ('B' - 'A'),
	PRINT_TIME = 1 << ('C' - 'A'),
	PRINT_LOG = 1 << ('D' - 'A'),
	SAVE_RESUTL = 1 << ('F' - 'A'),
	CHCECK_RESUTL = 1 << ('A' - 'A')
};

enum ECMD
{
//	COUNT_BITS = 0,
	COPY_BITS = 1,
	SAVE_PRIME,
	SAVE_BYTE,
	SAVE_BYTEGAP,
	FIND_MAXGAP,
	PCALL_BACK
};

enum EBITMASK
{
	BITZ = 0,
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
	uint L1Index;
	uint L2Size;
	uint L2Maxp;
	uint L2Index;
}
CpuCache =
{
	L1_DCACHE_SIZE * (WHEEL30 << 10),
	(L1_DCACHE_SIZE << 10) / ERATSMALL,
	0,
	L2_DCACHE_SIZE * (WHEEL30 << 10),
	(L2_DCACHE_SIZE << 10) / ERATMEDIUM,
};

//config
static struct
{
	uint Flag;
	//print calculating time
	uint Progress;
	//sieve size
	uint SieveSize;
	//sieve_size / ERATBIG
	uint MinPrime;
}
Config =
{
	PRINT_RET | PRINT_TIME | PRINT_LOG,
	16 - 1,
	MAX_SIEVE * (WHEEL30 << 10),
	MAX_SIEVE * (WHEEL30 << 10) / ERATBIG,
};

struct Cmd
{
	ECMD Oper;
	uchar* Data;
	uint64 Primes;
};

struct WheelPrime
{
	//[0 - 7]: p % wheel, [8 - 31]: p / sieve_size
	uint Wp;
	//[0 - 7]: index % sieve_size, [8 - 31]: index / 30
	uint SieveIndex;
};

struct Stock
{
	WheelPrime* Wheel; //read only
	Stock* Next;
//	uint WheelSize;
};

//each bucket contains multi Stock(by list)
//each stock contains multi wheelprime(by array)
//Head->Stock1->Stock2->....->Stockn
//
struct _Bucket
{
	WheelPrime* Wheel;//write only
	Stock* Head;
};

struct _BucketInfo
{
	uint CurIndex;
	uint MaxIndex;
	uint LoopSize;
	uint SieveSize;
	uint Log2Size;
	bool UseSieve;
};

static const uint MAX_BUCKET = 1 << 13; //> 10*2^32 / 2^18 * 30 = 10^14 / 3
static const uint WHEEL_BLOCK = 1 << 11; //12: 32k, best in [10 - 13]
static const uint MEM_BLOCK = (1 << 22) / WHEEL_BLOCK; //8M
static const uint MEM_SIZE  = WHEEL_BLOCK * sizeof(WheelPrime);
static const uint MAX_CACHE = (MAX_SIEVE + 2 * L1_DCACHE_SIZE) << 10;
static const uint MAX_STOCK = 203280221 / WHEEL_BLOCK + MAX_BUCKET; //pi(2^32) = 203280221
static const uint MAX_MEDIUM = (280738 / L2_DCACHE_SIZE) * (MAX_CACHE >> 10);

//thread data
static uint StockSize, AllStocks, MemSize;
static _BucketInfo BucketInfo;
static Stock MemStock [MAX_STOCK];
static Stock* StockPool [MAX_STOCK];
static _Bucket Bucket [MAX_BUCKET + 2];
static WheelPrime* StockMem [(1 << 17) / MEM_BLOCK]; //2G vm

//each segment sieve_index: (SegIndex[i] + start) % p[j] == 0
static WheelPrime MediumWheel[MAX_MEDIUM];

//presieved small prime number <=19 bit array.
//the crossing out bit module WHEEL30, the first
//16 bit of PreSieved map to
//----------------------------------------
//|01/1|07/1|11/1|13/1|17/1|19/1|23/1|29/1| = 0x1111 1111 = PreSieved[0]
//----------------------------------------
//|31/1|37/1|41/1|43/1|47/1|49/0|53/1|59/1| = 0x1101 1111 = PreSieved[1]
//----------------------------------------
static uchar PreSieved[PRIME_PRODUCT / WHEEL30];
//Prime (i + 1) th - (i)th difference
#if CHECK
static uchar Prime[MAX_MEDIUM * 50];
#else
static uchar Prime[MAX_MEDIUM * 1];
#endif

//position of least significant 1 bit of an integer
#if BIT_SCANF == 0
static uchar Lsb[1 << 16];
#endif

//number of bits 1 binary representation table in Range[0-2^16)
#if TREE2 == 0 && POPCNT == 0
static uchar WordNumBit1[1 << 16];
#endif

//cache small segment primes
static uint64 PiCache[1001];

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

static PrimeGap* WheelGap = NULL;

struct WheelElement
{
	uchar WheelIndex;
	uchar UnsetBit;
	uchar Correct;
	uchar NextMultiple;
};

struct WheelInit
{
	uchar WheelIndex;
	uchar UnsetBit;
};

typedef WheelElement WheelFirst;
static WheelElement Wheel210[48][64];
static WheelFirst WheelFirst210[WHEEL210][64];
static WheelInit WheelInit210[WHEEL210];

static WheelElement Wheel30[8][8];
static WheelFirst WheelFirst30[WHEEL30][8];
static WheelInit WheelInit30[WHEEL30];

static const uchar Pattern30[] =
{
	1  ,7  ,11 ,13 ,17 ,19 ,23 ,29 ,
	31 ,37 ,41 ,43 ,47 ,49 ,53 ,59 ,
	61 ,67 ,71 ,73 ,77 ,79 ,83 ,89 ,
	91 ,97 ,101,103,107,109,113,119,
	121,127,131,133,137,139,143,149,
	151,157,161,163,167,169,173,179,
	181,187,191,193,197,199,203,209,
	211,217,221,223,227,229,233,239,
};

typedef void (*call_back)(uint64, uint64);
static int segmentedSieve(uint64 start, uint sieve_size, Cmd*);
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

	for (int i = 0; i < 32; i ++) {
		root <<= 1;
		rem = ((rem << 2) + (uint)(n >> 62));
		n <<= 2;
		divisor = (root << 1) + 1;
		if (divisor <= rem) {
			rem -= divisor;
			root ++;
		}
	}

	return (uint)root;
}

#if BIT_SCANF
inline static uint bitScanForward(stype n)
{
#if _MSC_VER
	unsigned long index;
	#if X86_64
	_BitScanForward64(&index, n);
	#else
	_BitScanForward(&index, n);
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

//return n % p < 2^32
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
		str ++;
	}

	while (isdigit(*str)) {
		ret = ret * 10 + *str ++ - '0';
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

#define CHECK_SET(pbit, orbit)  if (pbit <= pend) *pbit |= orbit
static void crossOffWheelFactor(uchar* ppbeg[], const uchar* pend, const uint p)
{
	uchar* ps0 = ppbeg[0], *ps1 = ppbeg[1];
	uchar* ps2 = ppbeg[2], *ps3 = ppbeg[3];
	while (ps3 <= pend) {
		*ps0 |= BIT0, ps0 += p;
		*ps1 |= BIT1, ps1 += p;
		*ps2 |= BIT2, ps2 += p;
		*ps3 |= BIT3, ps3 += p;
	}
	CHECK_SET(ps0, BIT0); CHECK_SET(ps1, BIT1); CHECK_SET(ps2, BIT2);

	ps0 = ppbeg[4], ps1 = ppbeg[5];
	ps2 = ppbeg[6], ps3 = ppbeg[7];
	while (ps3 <= pend) {
		*ps0 |= BIT4, ps0 += p;
		*ps1 |= BIT5, ps1 += p;
		*ps2 |= BIT6, ps2 += p;
		*ps3 |= BIT7, ps3 += p;
	}
	CHECK_SET(ps0, BIT4); CHECK_SET(ps1, BIT5); CHECK_SET(ps2, BIT6);
}

//30% fast on old cpu pentium 4
static uchar* crossOffWheelFactor2(uchar* pbeg, const uchar* pend, const uint step)
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

static void crossOff4Factor(uchar* ppbeg[], const uchar* pend, const uint mask, const uint p)
{
#if 0
	typedef struct umask {
		uchar mask0;
		uchar mask1;
		uchar mask2;
		uchar mask3;
	} smask;
	smask umask = {mask, mask >> 8, mask >> 16, mask >> 24};
#endif

	uchar* ps0 = ppbeg[0], *ps1 = ppbeg[1];
	uchar* ps2 = ppbeg[2], *ps3 = ppbeg[3];

	while (ps3 <= pend) {
#if 0
		*ps0 |= umask.mask0; ps0 += p;
		*ps1 |= umask.mask1; ps1 += p;
		*ps2 |= umask.mask2; ps2 += p;
		*ps3 |= umask.mask3; ps3 += p;
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
}

inline static void crossOff2Factor(uchar* ps0, uchar* ps1, const uchar* pend, const ushort smask, const int p)
{
	const uchar masks1 = smask >> 8;
	const uchar masks0 = (uchar)smask;

	while (ps1 <= pend) {
		*ps1 |= masks1; ps1 += p;
		*ps0 |= masks0; ps0 += p;
	}

	if (ps0 <= pend)
		*ps0 |= masks0;
}

//popcnt instruction : INTEL ix/SSE4.2, AMD Phonem/SSE4A
//Use of the SSE 4.2 POPCNT instruction for fast bit counting
#if POPCNT
inline static int countBitsPopcnt(const uint64 n)
{
#ifdef X86_64
	return (int)_mm_popcnt_u64(n);
#else
	return _mm_popcnt_u32(n) + _mm_popcnt_u32((uint)(n >> 32));
#endif
}
#endif

#if POPCNT == 0 && TREE2 == 0
inline static int countBitsTable(const uint64 n)
{
	uint hig = (uint)(n >> 32), low = (uint)n;
//	return __builtin_popcount(hig) + __builtin_popcount(low)
//	return __builtin_popcountll(n);
	return WordNumBit1[(ushort)hig] + WordNumBit1[(ushort)low] +
		WordNumBit1[hig >> 16] + WordNumBit1[low >> 16];
}
#endif

inline static int countBitsTree(uint64 n)
{
	n -= (n >> 1) & 0x5555555555555555;
	n = (n & 0x3333333333333333) + ((n >> 2) & 0x3333333333333333);
	n = (n + (n >> 4)) & 0x0F0F0F0F0F0F0F0F;
	n += n >> 8;
	n += n >> 16;
	n += n >> 32;
	return ((uint)n) & 0x00000000FF;
}

//count number of bit 0 in binary representation of array
//bit 0 mean it's a prime position
static int countBit0sArray(const uint64 bitarray[], const int bitsize)
{
	int loops = bitsize >> 6;
	int bit1s = (1 + loops) << 6;

	while (loops -- >= 0) {
#if POPCNT
		bit1s -= countBitsPopcnt(*bitarray ++);
#elif TREE2
		bit1s -= countBitsTree(*bitarray ++);
#else
		bit1s -= countBitsTable(*bitarray ++);
#endif
	}

	return bit1s;
}

//Pre-sieve multiples of small primes <= limit (e.g. 19).
//Resets the sieve array (resets bits to 1) of SieveOfEratosthenes
//objects after each sieved segment and removes the multiples of
//small primes without sieving.
static int preSieve(uchar bitarray[], const uint64 start, const int sieve_size)
{
	assert(start % WHEEL30 == 0);
	const int offset = (int)(start % PRIME_PRODUCT) / WHEEL30;
	const int bits = sieve_size / WHEEL30 * 8 + WheelInit30[sieve_size % WHEEL30].WheelIndex;
	const int bytes = (bits + 7) / 8;

	if (offset + bytes < sizeof(PreSieved)) {
		memcpy(bitarray, PreSieved + offset, bytes);
	} else {
		memcpy(bitarray, PreSieved + offset, sizeof(PreSieved) - offset);
		memcpy(bitarray + sizeof(PreSieved) - offset, PreSieved, bytes + offset - sizeof(PreSieved));
	}

	if (start < WHEEL30) {
		//set bit position 1, 2 for prime 1, 7
		bitarray[0] = BIT0 | BIT1;
	}

	//pack the last 2 dword with bit 1
	uint* pdata = (uint*)bitarray + (bits >> 5);
	pdata[0] |= ~((1u << (bits & 31)) - 1);
	pdata[1] = pdata[2] = ~BITZ;

	return bytes;
}

static int sieveBytes(const uint64 start, const int sieve_size)
{
	int bits = sieve_size + (int)(start % WHEEL30);
	bits = bits / WHEEL30 * 8 + WheelInit30[bits % WHEEL30].WheelIndex;
	return (bits + 7) / 8;
}


//print prime from bit buffer
static void printPrime(uint64 primes, uint64 prime)
{
	printf("%llu %llu\n", primes, prime);
}

static int callBack(const stype bitarray[], uint64 sieve_index, const int size, uint64 sum_prime, call_back func)
{
	int primes = 0;

	for (int bi = 0; bi <= size; bi ++) {
		stype mask = ~bitarray[bi];
		while (mask > 0) {
			func(++primes + sum_prime, sieve_index + PRIME_OFFSET(mask));
			mask &= mask - 1;
		}
		sieve_index += WHEEL30 * sizeof(mask);
	}

	return primes;
}

//get prime from bit buffer
static int savePrime(const stype bitarray[], uint64 sieve_index, const int size, uint64* prime)
{
	int primes = 0;

	for (int bi = 0; bi <= size; bi ++) {
		stype mask = ~bitarray[bi];
		while (mask > 0) {
			prime[primes ++] = sieve_index + PRIME_OFFSET(mask);
			mask &= mask - 1;
		}
		sieve_index += WHEEL30 * sizeof(mask);
	}

	return primes;
}

static int savePrimeByte(const stype bitarray[], uint sieve_index, const int size, uchar* prime)
{
	int primes = 0;
	int lastp = sieve_index + (char)prime[0];

	for (int bi = 0; bi < size; bi ++) {
		stype mask = ~bitarray[bi];
		while (mask > 0) {
			int curp = sieve_index + PRIME_OFFSET(mask);
			mask &= mask - 1;
#if CHECK
			prime[primes ++] = curp - lastp + (curp - lastp) / 256;
#else
			prime[primes ++] = curp - lastp;
#endif
			lastp = curp;
		}
		sieve_index += WHEEL30 * sizeof(mask);
	}

	prime[primes] = lastp - sieve_index;

	return primes;
}

//get prime from bit buffer
//word_size > 1476 / 60
static int findPrimeGap(const ushort bitarray[], uint64 sieve_index, const int word_size, uchar* prime)
{
	uint64* result = (uint64*)prime;
	ushort pdiff = (ushort)result[0];
	ushort max_gap = (ushort)result[1];
	int i = 0, skip_words = max_gap / (WHEEL30 * 2) - 1;
	if (skip_words < 1)
		skip_words = 1;

	//find first
	for (i = 0; pdiff == 0; i ++) {
		initWheelGap();
		pdiff = WheelGap[~bitarray[i]].End;//bug for small range
	}

	for (; i < word_size; i ++) {
		const ushort mask = ~bitarray[i];
		if (mask == 0) {
			pdiff += WHEEL30 * sizeof(mask);
			continue ;
		}

		const ushort gap = pdiff + WheelGap[mask].Beg;
		if (gap > max_gap) { //
			max_gap = gap;
			result[2] = sieve_index + i * WHEEL30 * sizeof(mask) + PRIME_OFFSET(mask);
		}
		uint next = i + skip_words;
		if (bitarray[next] != 0xffff && next < word_size) {
			i = next - 1;
			pdiff = 0;
		} else {
			pdiff = WheelGap[mask].End;
		}
	}

	result[0] = pdiff;
	result[1] = max_gap;

	return 0;
}

static int savePrimeGap(ushort bitarray[], const uint64 start, const int word_size, uchar* prime)
{
	ushort pdiff = *(char*)prime;
	int primes = 0;
	initWheelGap();

	for (int i = 0; i < word_size; i ++) {
		const short mask = ~bitarray[i];
		if (mask == 0) {
			pdiff += WHEEL30 * sizeof(mask);
			continue ;
		}

		ushort gap = pdiff + WheelGap[mask].Beg;
		*prime = gap;
		if (gap > 255) {
			*(ushort*)(prime ++) = gap + 1;
		}
		*(uint64*)(prime + 1) = *(uint64*)(WheelGap[mask].Gap + 0);
		pdiff = WheelGap[mask].End;
		uchar bits = WheelGap[mask].Bits;
		if (bits > 8)
			*(uint*)(prime + 9) = *(uint*)(WheelGap[mask].Gap + 8);
		prime += bits;
		primes += bits;
	}

	prime[0] = pdiff;
	prime[1] = 0;

	return primes;
}

static void newWheelBlock(uint blocks)
{
	//dynamic memory use
	WheelPrime *pwheel = (WheelPrime*) malloc((blocks + 1) * MEM_SIZE);
	StockMem[MemSize++] = pwheel;
	assert(pwheel && MemSize < sizeof(StockMem) / sizeof(StockMem[0]));

	pwheel = (WheelPrime*)((size_t)pwheel + MEM_SIZE - (size_t)pwheel % MEM_SIZE);
	for (int i = 0; i < blocks; i ++) {
		StockPool[i + StockSize]->Wheel = pwheel + WHEEL_BLOCK * i;
		StockPool[i + StockSize]->Next = NULL;
	}

	StockSize += blocks;
	AllStocks += blocks;
}

static int initBucketInfo(const uint sieve_size, const uint sqrtp, uint64 range)
{
	BucketInfo.CurIndex = 0;
	assert (range / sieve_size < -1u);
	BucketInfo.MaxIndex = range / sieve_size + 1;//overflow for large range

	//maxp wheel 210 pattern, max pattern difference is 10
	BucketInfo.LoopSize = (2 << ilog2((uint64)sqrtp * 10 / sieve_size + 2)) - 1;
	assert (BucketInfo.LoopSize < MAX_BUCKET);

	BucketInfo.Log2Size = ilog2(sieve_size / WHEEL30);
	BucketInfo.SieveSize = (1 << BucketInfo.Log2Size) - 1;

	const uint pix = (uint)(sqrtp / log((double)sqrtp) * (1 + 1.200 / log((double)sqrtp)));

	//203280221 >> 12
	int blocks = pix / WHEEL_BLOCK + MIN(BucketInfo.LoopSize, BucketInfo.MaxIndex) * 2 / 3;
	if (range < sqrtp / 4)
		blocks >>= 2;
	else if (range < sqrtp / 2)
		blocks >>= 1;

	assert(blocks < sizeof(MemStock) / sizeof(MemStock[0]));
	for (int j = 0; j < blocks; j ++) {
		StockPool[j] = MemStock + j;
	}

#if 1
	for (int i = 0; i < blocks; i += MEM_BLOCK) {
		if (i + MEM_BLOCK > blocks)
			newWheelBlock(blocks - i);
		else
			newWheelBlock(MEM_BLOCK);
	}
#else
	//dynamic memory use
	WheelPrime *pwheel = (WheelPrime*) malloc((blocks + 1) * MEM_SIZE);
	assert(pwheel);
	StockPool[0]->Wheel = pwheel;

	pwheel = (WheelPrime*)((size_t)pwheel + MEM_SIZE - (size_t)pwheel % MEM_SIZE);
	for (int i = 1; i < blocks; i ++) {
		StockPool[i]->Wheel = pwheel + WHEEL_BLOCK * (i - 1);
		StockPool[i]->Next = NULL;
	}
	AllStocks = StockSize = blocks - 0;
#endif

	return pix;
}

static void pushBucket(const uint sieve_index, const uint wp, const uchar wheel_index)
{
	const uint next_bucket = (sieve_index >> BucketInfo.Log2Size) + BucketInfo.CurIndex;
	if (next_bucket < BucketInfo.MaxIndex) {
		_Bucket* pbucket = Bucket + (next_bucket & BucketInfo.LoopSize);
		WheelPrime* wheel = pbucket->Wheel ++;
		if ((size_t)wheel % MEM_SIZE == 0) {
//			if (StockSize == 0) newWheelBlock(MEM_BLOCK);
			Stock* pstock = StockPool[-- StockSize];
			pbucket->Wheel = pstock->Wheel;
			wheel = pbucket->Wheel ++;
			pstock->Next = pbucket->Head;
			pbucket->Head = pstock;
		}

		wheel->Wp = wp;
		wheel->SieveIndex = ((sieve_index & BucketInfo.SieveSize) << MPART) | wheel_index;
	}
}

static void initMediumWheel(const uint sieve_size, const uint maxp, const uint64 start, const bool check)
{
	uint j = CpuCache.L1Index, p = CpuCache.L1Maxp;
	uint minp = MIN(maxp, CpuCache.L2Maxp);
	for (; p < minp; PRIME_DIFF(p, j)) {
		uint medsieve = CpuCache.L2Size;
		const uint64 p2 = (uint64)p * p;
		uint sieve_index = p - (uint)(start % p);
		if (check)
			medsieve = sieve_size;

		if (p2 >= start) {
			sieve_index = (p2 - start) % medsieve;
		}

		const uint wi = WheelInit30[p % WHEEL30].WheelIndex;
		WheelFirst wf = WheelFirst30[sieve_index % WHEEL30][wi];
		sieve_index = (sieve_index + wf.Correct * p) / WHEEL30;

		if (p2 > start + medsieve) {
			sieve_index %= (medsieve / WHEEL30);
		}

		MediumWheel[j].Wp = (p / WHEEL30 << MPART) + wi;
		MediumWheel[j].SieveIndex = (sieve_index << MPART) | wf.WheelIndex;
	}

	for (; p < maxp; PRIME_DIFF(p, j)) {
		const uint medsieve = sieve_size;
		const uint64 p2 = (uint64)p * p;
#if WHEEL == WHEEL30
		uint sieve_index = p - (uint)(start % p);
#else
		uint64 sieve_index = p - (uint)(start % p);
#endif

		if (p2 >= start) {
#if WHEEL == WHEEL30
			sieve_index = (p2 - start) % medsieve;
#else
			sieve_index = p2 - start;
#endif
		}

		const uint wi = WHEEL_INIT[p % WHEEL].WheelIndex;
		WheelFirst wf = WHEEL_FIRST[sieve_index % WHEEL][wi];
		sieve_index = (sieve_index + wf.Correct * p) / WHEEL30;

		if (p2 > start + medsieve) {
			sieve_index %= (medsieve / WHEEL30);
		}

		MediumWheel[j].Wp = (p / WHEEL << MPART) + wi;
		MediumWheel[j].SieveIndex = (sieve_index << MPART) | wf.WheelIndex;
	}

//	MediumWheel[j].Wp = -1u;
	assert(j < sizeof(MediumWheel) / sizeof(MediumWheel[0]));
}

static void initWheelBucket(uint offset, const uint sqrtp, const uint64 start, const uint64 range)
{
	uint nextp = 0;
	uint64 remp = 0;
	if ((uint64)offset * offset >= start) {
		offset = sqrt((double)start) + 0;
	}

	for (uint segsize = CpuCache.L2Size / 1; offset < sqrtp; offset += segsize) {
		if (segsize > sqrtp - offset)
			segsize = sqrtp - offset;

		//Cmd cmd = {COPY_BITS, 0, (uchar*)bitarray};
		Cmd cmd = {COPY_BITS, 0, 0}; //no copy but no save
		segmentedSieve(offset, segsize, &cmd);

		stype mask = 0, *bitarray = (stype*)cmd.Data;
		uint p2 = offset - offset % WHEEL30 - sizeof(mask) * WHEEL30;

#if CPU == INTEL
		const uint pmax = 1024 << 1;
#else
		const uint pmax = ((uint64)offset * offset / (uint)(start >> 32)) - 0;
#endif

		for (int j = 0; j <= 1 + ((int)cmd.Primes) / sizeof(mask); ) {
			if (mask == 0) {
				mask = ~bitarray[j ++];
				p2 += sizeof(mask) * WHEEL30;
				continue ;
			}

			const uint p = p2 + PRIME_OFFSET(mask);
			mask &= mask - 1;

#if ASM_X86
			if (p > nextp) {
				remp = start / (nextp = p + pmax);
				if (p > nextp) { remp = start / (nextp = -1u); }
			}

			uint64 reminder = start - remp * p;
//			if ((uint)(reminder >> 32) > p) { reminder -= (uint64)p << 32; }
			uint sieve_index = p - fastMod(reminder, p);
#else
			uint sieve_index = p - (uint)(start % p);
#endif
//			if (sieve_index % 2 == 0) sieve_index += p;
#if CPU == INTEL
			if (sieve_index >= range)
				continue ;
#endif

			const uint wp = (p / WHEEL210 << WPART) + WheelInit210[p % WHEEL210].WheelIndex;
			const uint module_wheel = sieve_index % WHEEL210;
			const WheelFirst wf = WheelFirst210[module_wheel][wp % (1 << WPART)];
#if CPU == INTEL || X86_64 == 0
			sieve_index = (sieve_index / WHEEL210 + (wp >> WPART) * wf.Correct) * (WHEEL210 / WHEEL30);
			sieve_index += (wf.Correct * (p % WHEEL210) + module_wheel) / WHEEL30;
#else
			sieve_index = (sieve_index + wf.Correct * (uint64)p) / WHEEL30;
#endif
			pushBucket(sieve_index, wp, wf.WheelIndex);
		}
	}

	assert(StockSize > 0);
}

inline static void
sieveSmall0(uchar bitarray[], const uchar* pend, const uint p, uint sieve_index, uint multiples)
{
	uchar* ppbeg[8];
	for (int i = 0; i < 8; i ++) {
		ppbeg[WheelInit30[sieve_index % WHEEL30].WheelIndex] = bitarray + sieve_index / WHEEL30;
		sieve_index += (multiples % 16) * p; multiples /= 16;
	}
	crossOffWheelFactor(ppbeg, pend, p);
}

inline static void
sieveSmall1(uchar bitarray[], const uchar* pend, const uint p, uint sieve_index, uint multiples)
{
	for (int i = 0; i < 4; i ++) {
		uchar* ps0 = bitarray + sieve_index / WHEEL30;
		ushort smask = WheelInit30[sieve_index % WHEEL30].UnsetBit;
		sieve_index += (multiples % 16) * p; multiples /= 16;

		uchar* ps1 = bitarray + sieve_index / WHEEL30;
		smask |= WheelInit30[sieve_index % WHEEL30].UnsetBit << 8;
		sieve_index += (multiples % 16) * p; multiples /= 16;
		crossOff2Factor(ps0, ps1, pend, smask, p);
	}
}

inline static void
sieveSmall2(uchar bitarray[], const uchar* pend, const uint p, uint sieve_index, uint multiples)
{
	uchar* ppbeg[8], mask[8];
	for (int i = 0; i < 8; i ++) {
		ppbeg[i] = bitarray + sieve_index / WHEEL30;
		mask[i] = WheelInit30[sieve_index % WHEEL30].UnsetBit;
		sieve_index += (multiples % 16) * p; multiples /= 16;
	}
	crossOff4Factor(ppbeg + 0, pend, *(uint*)(mask + 0), p);
	crossOff4Factor(ppbeg + 4, pend, *(uint*)(mask + 4), p);
}

inline static void
sieveSmall3(uchar bitarray[], const uchar* pend, const uint p, uint sieve_index, uint multiples)
{
	for (int i = 0; i < 8; i ++) {
		const uchar mask = WheelInit30[sieve_index % WHEEL30].UnsetBit;
		bitarray[sieve_index / WHEEL30] |= mask;
		if (mask == BIT0)
			break;
		sieve_index += (multiples % 16) * p; multiples /= 16;
	}
	crossOffWheelFactor2(bitarray + sieve_index / WHEEL30, pend + 1, p);
}

inline static WheelElement*
sieveSmall5(uchar bitarray[], const uchar* pend, const uint wi, const uint pi, WheelElement* wheel)
{
	const uint p = wi * WHEEL30 + Pattern30[pi];
	WheelElement* wdata = Wheel30[pi];

	for (int i = 0; i < 4; i ++) {
		uchar* ps0 = bitarray;
		ushort smask = wheel->UnsetBit;
		bitarray += wheel->Correct + wheel->NextMultiple * wi;
		wheel = wdata + wheel->WheelIndex;

		uchar* ps1 = bitarray;
		smask |= ((ushort)wheel->UnsetBit) << 8;
		bitarray += wheel->Correct + wheel->NextMultiple * wi;
		if (i < 3)
			wheel = wdata + wheel->WheelIndex;

		crossOff2Factor(ps0, ps1, pend, smask, p);
	}

	return wheel;
}

static void eratSieveL1(uchar bitarray[], const uint64 start, const int segsize, uint maxp)
{
	const uchar* pend = bitarray + segsize / WHEEL30;

	if ((start + segsize) < ((uint64)maxp) * maxp) {
		maxp = (uint)sqrt((double)start + segsize) + 1;
	}

	for (uint p = Prime[0], j = 8 + PRIME_PRODUCT / 9699690; p < maxp; PRIME_DIFF(p, j)) {
		uint sieve_index = p - (uint)(start % p);
		if (start <= p) {
			sieve_index = p * p - (uint)start;
		}

		const WheelFirst wf = WheelFirst30[sieve_index % WHEEL30][WheelInit30[p % WHEEL30].WheelIndex];
		const uint multiples = (WHEEL_SKIP >> wf.NextMultiple * 4) | (WHEEL_SKIP << (32 - 4 * wf.NextMultiple));
		sieve_index += wf.Correct * p;
#if __GNUC__ && CPU == INTEL
//		sieveSmall3(bitarray, pend, p, sieve_index, multiples);
		sieveSmall0(bitarray, pend, p, sieve_index, multiples);
#elif _MSC_VER
		sieveSmall2(bitarray, pend, p, sieve_index, multiples);
#else
		sieveSmall1(bitarray, pend, p, sieve_index, multiples);
#endif
	}
}

static void eratSieveSmall(uchar bitarray[], const uint64 start, const uint sieve_size)
{
	const uint maxp = CpuCache.L1Maxp;
	uint segsize = CpuCache.L1Size;

	//preSieve(bitarray, start, sieve_size);
	//#pragma omp parallel for num_threads(2)
	for (uint sieve_index = 0; sieve_index < sieve_size; sieve_index += segsize) {
		if (segsize + sieve_index > sieve_size)
			segsize = sieve_size - sieve_index;
		preSieve(bitarray + sieve_index / WHEEL30, start + sieve_index, segsize);
		eratSieveL1(bitarray + sieve_index / WHEEL30, start + sieve_index, segsize, maxp);
	}
}

static void eratSieveMedium2(uchar bitarray[], const uint64 start, const uint sieve_size, uint sqrtp)
{
	const uint l1maxp = CpuCache.L1Maxp;
	if (sqrtp < l1maxp)
		return ;

	uint j = CpuCache.L1Index, p = l1maxp;
	const uint mins = MIN(sqrtp, sieve_size);
	for (const uchar* pend = bitarray + sieve_size / WHEEL30; p < mins; PRIME_DIFF(p, j)) {
		uint sieve_index = p - (uint)(start % p);
		const WheelFirst wf = WheelFirst30[sieve_index % WHEEL30][WheelInit30[p % WHEEL30].WheelIndex];
		const uint multiples = (WHEEL_SKIP >> wf.NextMultiple * 4) | (WHEEL_SKIP << (32 - 4 * wf.NextMultiple));
		sieve_index += wf.Correct * p;

//		sieveSmall0(bitarray, pend, p, sieve_index, multiples);
		sieveSmall1(bitarray, pend, p, sieve_index, multiples);
//		sieveSmall2(bitarray, pend, p, sieve_index, multiples);
//		sieveSmall3(bitarray, pend, p, sieve_index, multiples);
	}

#if CHECK
	const uint pmax = (uint)(start >> 32) + 1;
	uint64 remp = 0;
	for (uint nextp = 0; p < sqrtp; PRIME_DIFF(p, j)) {

		PRIME_ADD_DIFF(p);
		if (p > nextp) {
			remp = start / (nextp = p + ((uint64)p * p) / pmax);
			if (p > nextp) remp = start / (nextp = -1u);
		}
		uint sieve_index = p - fastMod(start - remp * p, p);
//		uint sieve_index = p - start % p;
		if (sieve_index <= sieve_size) {
			bitarray[sieve_index / WHEEL30] |= WheelInit30[sieve_index % WHEEL30].UnsetBit;
		}
	}
#endif
}

static void eratSieveL2(uchar bitarray[], const uint64 start, uint sieve_size)
{
	uint maxp = CpuCache.L2Maxp;
	if ((start + sieve_size) < (uint64)maxp * maxp) {
		maxp = (uint)sqrt((double)start + sieve_size) + 1;
	}

	WheelPrime* pwheel = MediumWheel + CpuCache.L1Index;

	maxp = (maxp / WHEEL30 << MPART) + WheelInit30[maxp % WHEEL30].WheelIndex;
	sieve_size = sieve_size / WHEEL30 + WheelInit30[sieve_size % WHEEL30].WheelIndex;
	const uchar* pend = bitarray + sieve_size + 1;

	for (uint wp1 = pwheel->Wp; wp1 < maxp; wp1 = pwheel->Wp) {

		const uint wi = wp1 >> MPART;
		const uchar pi = wp1 % (1 << MPART);
		const uint p = wi * WHEEL30 + Pattern30[pi];
		uint sieve_index = pwheel->SieveIndex >> MPART;
		WheelElement* wdata = Wheel30[pi];
		WheelElement* wheel = wdata + pwheel->SieveIndex % (1 << MPART);

		wheel = sieveSmall5(bitarray + sieve_index, pend, wi, pi, wheel);

		sieve_index += (sieve_size - sieve_index) / p * p;
		while (sieve_index < sieve_size) {
			wheel = wdata + wheel->WheelIndex;
			sieve_index += wheel->Correct + wheel->NextMultiple * wi;
		}
		pwheel++ ->SieveIndex = (sieve_index - sieve_size) << MPART | wheel->WheelIndex;
		if (wp1 > pwheel->Wp)
			break;
	}
}

//core code of this algorithm for large range
//sieve prime multiples in [start, start + sieve_size)
static void eratSieveMedium(uchar bitarray[], const uint64 start, uint sieve_size, uint maxp)
{
	if ((start + sieve_size) < (uint64)maxp * maxp) {
		maxp = (uint)sqrt((double)start + sieve_size) + 1;
	}

	WheelPrime* pwheel = MediumWheel + CpuCache.L2Index;
	maxp = (maxp / WHEEL << MPART) + WHEEL_INIT[maxp % WHEEL].WheelIndex;
	sieve_size = sieve_size / WHEEL30 + WheelInit30[sieve_size % WHEEL30].WheelIndex;

	for (uint wp = pwheel->Wp; wp < maxp; wp = pwheel->Wp) {
		uint sieve_index = pwheel->SieveIndex >> MPART;
		if (sieve_index >= sieve_size) {
			pwheel++ ->SieveIndex -= sieve_size << MPART;
			continue;
		}

		const uint wi = wp >> MPART;
		WheelElement* wdata = WHEEL_MAP[wp % (1 << MPART)];
#if 0
		WheelElement wheel = wdata[pwheel->SieveIndex % (1 << MPART)];
		bitarray[sieve_index] |= wheel.UnsetBit;
		sieve_index += wheel.Correct + wheel.NextMultiple * wi;
		while (sieve_index < sieve_size) {
			wheel = wdata[wheel.WheelIndex];
			bitarray[sieve_index] |= wheel.UnsetBit;
			sieve_index += wheel.Correct + wheel.NextMultiple * wi;
		}
		pwheel++ ->SieveIndex = (sieve_index - sieve_size) << MPART | wheel.WheelIndex;
#else
		uint wheel = *(uint*)(wdata + pwheel->SieveIndex % (1 << MPART));
		bitarray[sieve_index] |= wheel >> 8;
		sieve_index += (uchar)(wheel >> 16) + (wheel >> 24) * wi;
		while (sieve_index < sieve_size) {
			wheel = *(uint*)(wdata + ((uchar)wheel));
			bitarray[sieve_index] |= wheel >> 8;
			sieve_index += (uchar)(wheel >> 16) + (wheel >> 24) * wi;
		}
		pwheel++ ->SieveIndex = (sieve_index - sieve_size) << MPART | ((uchar)wheel);
#endif
	}
}

//This implementation uses a sieve array with WHEEL30 numbers per byte and
//a modulo 210 wheel that skips multiples of 2, 3, 5 and 7.
static void eratSieveBucket(uchar bitarray[], const uint sieve_size)
{
	_Bucket* pbucket = Bucket + (BucketInfo.CurIndex & BucketInfo.LoopSize);

	int loops = (uint64)pbucket->Wheel % MEM_SIZE / sizeof(WheelPrime);
	if (loops == 0) { loops = WHEEL_BLOCK; }

//	for (Stock* phead = (Stock*)pbucket; phead = phead->Next; loops = WHEEL_BLOCK) {
	for (Stock* phead = pbucket->Head, *pnext = phead; phead = pnext; loops = WHEEL_BLOCK) {
		pnext = phead->Next;

		WheelPrime* cur_wheel = phead->Wheel;

		while (loops --) {
			const uint wp = cur_wheel->Wp;
			uint sieve_index = cur_wheel ++ ->SieveIndex;

#if CPU == INTEL
			WheelElement* wheel = &Wheel210[wp % (1 << WPART)][sieve_index % (1 << MPART)];
			/// 135/1760
			bitarray[sieve_index >>= MPART] |= wheel->UnsetBit;
			sieve_index += wheel->Correct + wheel->NextMultiple * (wp >> WPART);
#if ERATBIG > 2
			if (sieve_index < sieve_size) {
				wheel = &Wheel210[wp % (1 << WPART)][wheel->WheelIndex];
				bitarray[sieve_index] |= wheel->UnsetBit;
				sieve_index += wheel->Correct + wheel->NextMultiple * (wp >> WPART);
			}
#endif
			pushBucket(sieve_index, wp, wheel->WheelIndex);
#else
			uint wheel = *(uint*)&Wheel210[wp % (1 << WPART)][sieve_index % (1 << MPART)];
			bitarray[sieve_index >>= MPART] |= wheel >> 8;
			sieve_index += (uchar)(wheel >> 16) + (wheel >> 24) * (wp >> WPART);
#if ERATBIG > 2
			if (sieve_index < sieve_size) {
				wheel = *(uint*)&Wheel210[wp % (1 << WPART)][(uchar)wheel];
				bitarray[sieve_index] |= wheel >> 8;
				sieve_index += (uchar)(wheel >> 16) + (wheel >> 24) * (wp >> WPART);
			}
#endif
			pushBucket(sieve_index, wp, wheel);
#endif
		}
		StockPool[StockSize ++] = phead;
	}

	pbucket->Wheel = NULL;
	pbucket->Head = NULL;
	BucketInfo.CurIndex ++;
}

static int segmentProcessed(uchar bitarray[], const uint64 start, const int sieve_size, Cmd* cmd)
{
	assert(start % WHEEL30 == 0);
	const uint bytes = sieveBytes(start, sieve_size);
	if (cmd == NULL)
		return countBit0sArray((uint64*)bitarray, bytes * 8);

	const int oper = cmd->Oper;
	int primes = 0, dwords = (1 + bytes) / sizeof(stype);

	if (COPY_BITS == oper) {
		cmd->Data = bitarray;
		primes = bytes; //if (cmd->Data) memcpy(cmd->Data, bitarray, bytes + 8); else
	} else if (SAVE_BYTE == oper) {
		primes = savePrimeByte((stype*)bitarray, start, dwords, cmd->Data + cmd->Primes);
//	} else if (SAVE_BYTEGAP == oper) {
//		primes = savePrimeGap(bitarray, start, dwords / 2, cmd->Data + cmd->Primes);
	} else if (FIND_MAXGAP == oper) {
		primes = findPrimeGap((ushort*)bitarray, start, (bytes + 1) / 2, cmd->Data + cmd->Primes);
	} else if (PCALL_BACK == oper) {
		primes = callBack((stype*)bitarray, start, dwords, cmd->Primes, (call_back)cmd->Data);
//	} else if (SAVE_PRIME == oper) {
//		primes = savePrime((stype*)bitarray, start, dwords, (uint64*)cmd->Data + cmd->Primes);
	}

	cmd->Primes += primes;
	return primes;
}

//sieve prime multiples in [start, start + sieve_size)
static int segmentedSieve(const uint64 start, const uint sieve_size, const uint wheel_offset, Cmd* cmd = NULL)
{
	//1.sieve small/medium prime factor
	uchar bitarray[MAX_CACHE];
	for (uint sieve_index = 0, segsize = CpuCache.L2Size; sieve_index < sieve_size; sieve_index += segsize) {
		if (segsize + sieve_index > sieve_size)
			segsize = sieve_size - sieve_index;

		//sieve p in [23, L1_DCACHE_SIZE / 30] 240/1560
		eratSieveSmall(bitarray + sieve_index / WHEEL30, start + sieve_index, segsize);
		//sieve p in [L1_DCACHE_SIZE / 30, L2_DCACHE_SIZE / 30] 147/1760 180/1560
		static double time_use = 0; double ts = getTime();
		eratSieveL2(bitarray + sieve_index / WHEEL30, start + sieve_index, segsize);
//		eratSieveMedium(bitarray + sieve_index / WHEEL30, start + sieve_index, segsize, CpuCache.L1Index, CpuCache.L2Maxp);
//		time_use += getTime() - ts; if (sieve_size != Config.SieveSize) { printf("eratSieveL2 time %.f ms\n", time_use); time_use = 0; }
	}

	const uint sqrtp = (uint)sqrt((double)start + sieve_size) + 1;
	const uint medium = MIN(sqrtp, Config.MinPrime);

	//2.sieve p in [L2_DCACHE_SIZE / 30, sieve_size / 2]
	//654/1760 560/1560
	static double time_use1 = 0; double ts = getTime();
	if (medium >= CpuCache.L2Maxp)
		eratSieveMedium(bitarray, start, sieve_size, medium);
//	time_use1 += getTime() - ts; if (sieve_size != Config.SieveSize) { printf("eratSieveMedium time %.f ms\n", time_use1); time_use1 = 0; }

	if (BucketInfo.UseSieve) {
		static double time_use2 = 0; double ts2 = getTime();
		//3.sieve p [sieve_size / 2, sqrtp]
		//630/1760, 610/1560
		eratSieveBucket(bitarray, Config.SieveSize / WHEEL30);
//		time_use2 += getTime() - ts2; if (sieve_size != Config.SieveSize) { printf("eratSieveBucket time %.f ms\n", time_use2); time_use2 = 0; }
	}

	if (wheel_offset > 0) {
		memset(bitarray, ~BITZ, wheel_offset / WHEEL30);
		bitarray[wheel_offset / WHEEL30] |= (1 << WheelInit30[wheel_offset % WHEEL30].WheelIndex) - 1;
	}

	return segmentProcessed(bitarray, start, sieve_size, cmd);
}

static int segmentedSieve(uint64 start, uint sieve_size, Cmd* cmd = NULL)
{
	const uint sqrtp = (uint)sqrt((double)start + sieve_size) + 1;
	const uint wheel_offset = start % WHEEL30;
	start -= wheel_offset;
	sieve_size += wheel_offset;

	static double time_use = 0; double ts = getTime();
	uchar bitarray[MAX_CACHE / 1];
	eratSieveSmall(bitarray, start, sieve_size);
	eratSieveMedium2(bitarray, start, sieve_size, sqrtp);

	bitarray[0] |= (1 << WheelInit30[wheel_offset].WheelIndex) - 1;
//	time_use += getTime() - ts; if (sieve_size % WHEEL30) { printf("eratSieveSqrt time %.f ms\n", time_use); time_use = 0; }
	return segmentProcessed(bitarray, start, sieve_size, cmd);
}

static void initMaxp(uint maxp)
{
	uint p = 11, k = 5;
	for (; p < maxp; PRIME_DIFF(p, k)) {
		if (p >= CpuCache.L1Maxp) {
			CpuCache.L1Index = k;
			CpuCache.L1Maxp = p;
			break;
		}
	}
	for (; p < maxp; PRIME_DIFF(p, k)) {
		if (p >= CpuCache.L2Maxp) {
			CpuCache.L2Index = k;
			CpuCache.L2Maxp = p;
			break;
		}
	}
}

static void setCpuSize(uint cdata)
{
	if ((cdata & (cdata - 1)) != 0) {
		cdata = 1 << ilog2(cdata);
	}

	if (cdata >= 16 && cdata < L2_DCACHE_SIZE) { //L1
		CpuCache.L1Size = cdata * (WHEEL30 << 10);
		CpuCache.L1Maxp = CpuCache.L1Size / WHEEL30 / ERATSMALL;
	} else if (cdata >= L2_DCACHE_SIZE && cdata <= 1024) { //L2
		CpuCache.L2Size = cdata * (WHEEL30 << 10);
		CpuCache.L2Maxp = CpuCache.L2Size / WHEEL30 / ERATMEDIUM;
	}

	initMaxp(CpuCache.L2Maxp + 512);
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

	sieve_size = (1 << ilog2(sieve_size / WHEEL30 + 1)) * WHEEL30;

	if (sieve_size != Config.SieveSize) {
		memset(PiCache, 0, sizeof(PiCache));
	}

	Config.MinPrime = sieve_size / ERATBIG + 1;

	return Config.SieveSize = sieve_size;
}

static int checkSmall(const uint64 start, const uint64 end, Cmd* cmd)
{
	int primes = 0;
	const uchar smallp[] = {2, 3, 5, 7};
	for (int i = 0; i < sizeof(smallp) / sizeof(smallp[0]); i ++) {
		if (start <= smallp[i] && smallp[i] <= end) {
			primes ++;
			if (cmd && cmd->Oper == PCALL_BACK) {
				(*(call_back)cmd->Data)(primes, smallp[i]);
				cmd->Primes += 1;
			} else if (cmd && cmd->Oper == SAVE_PRIME) {
				((uint64*)cmd->Data)[primes - 1] = smallp[i];
				cmd->Primes += 1;
			}
		}
	}

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

//calculate number of prime in Range[start, end] with start <= end
static uint64 pi2(const uint64 start, const uint64 end, Cmd* cmd = NULL)
{
	assert (start <= end);

	if (Config.Flag & SAVE_RESUTL) {
		freopen("prime.txt", "w", stdout);
		Cmd cmdbuf = {PCALL_BACK, (uchar*)printPrime, 0};
		cmd = &cmdbuf;
	}

	uint64 primes = checkSmall(start, end, cmd);

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
		if (sieve_size > end - start)
			sieve_size = uint(end - start) + (end & 1);

		const int seg = segmentedSieve(start, sieve_size, module_start, cmd);
		primes += seg;

#if CHECK
		const int ok = segmentedSieve(start + module_start, sieve_size - module_start + 1, cmd);
		if (ok != seg) {
			printf("[%u]: sqrt = %d, start = %I64d, seg = %d != %d = ok\n", (uint)si, (int)sqrt((double)start), start, seg, ok);
		}
#endif
		if ((si ++ & Config.Progress) == Config.Progress - 1) {
			double ratio = si * segs;
			printf(">> %02d%%, total time ~ %.2f sec, primes ~= %lld\r",
					((int)ratio) / 10, (getTime() - ts) / ratio, (uint64)(1000 * primes / ratio));
		}
		module_start = 0;
	}

	return primes;
}

static void printPiResult(const uint64 start, const uint64 end, uint64 primes, double ts)
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
		printf(" (%.2lf sec)\t\t", (getTime() - ts) / 1000.0);
	putchar('\n');
}

static int sievePrime(uchar prime[], uint n)
{
	static uint maxp = 100;
	if (n <= maxp) {
		return 0;
	}
	maxp = n;

	double ts = getTime( );
	int primes = 4;

#if 0
	prime[primes] = 23 - WHEEL30;
	Cmd cmd = {SAVE_BYTEGAP, prime, primes};
#else
	prime[primes] = 11 - 7;
	Cmd cmd = {SAVE_BYTE, prime, primes};
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

	initMaxp(maxp);

//	if (Config.Flag & PRINT_LOG)
//		dumpPrime(getTime() - ts, prime);

	return primes;
}

uint64 doSievePrime2(const uint64 start, const uint64 end, Cmd* cmd = NULL)
{
	double ts = getTime();
	const uint sqrtp = (uint)sqrt((double)end) + 1;

	sievePrime(Prime + 1, sqrtp);

	const uint64 primes = pi2(start, end, cmd);

	if ((!cmd || cmd->Oper != FIND_MAXGAP) && (Config.Flag & PRINT_RET))
		printPiResult(start, end, primes, ts);

	return primes;
}

uint64 doSievePrime(const uint64 start, const uint64 end, Cmd* cmd = NULL)
{
	assert (start <= end && (end - start) < 1e16);

	double ts = getTime();
	const uint sqrtp = (uint)sqrt((double)end) + 1;

	if (Config.Flag & SAVE_RESUTL) {
		freopen("prime.txt", "w", stdout);
		Cmd cmdbuf = {PCALL_BACK, (uchar*)printPrime, 0};
		cmd = &cmdbuf;
	}

	uint sieve_size = Config.SieveSize;
	uint64 primes = checkSmall(start, end, cmd);
	if (sieve_size < CpuCache.L2Size && sqrtp > Config.MinPrime && sqrtp > 10000000) {
		sieve_size = setSieveSize(MAX_SIEVE);
	}
	const uint minp = Config.MinPrime;
	const uint medium = MIN(sqrtp, minp);

#if CHECK
	sievePrime(Prime + 1, sqrtp);
#else
	sievePrime(Prime + 1, medium + 256);
#endif

	const uint64 segment_low = start - (start % WHEEL210);
	initMediumWheel(sieve_size, medium + 256, segment_low, sieve_size < CpuCache.L2Size);

	BucketInfo.UseSieve = sqrtp > minp && sqrtp >= CpuCache.L2Maxp;
	if (BucketInfo.UseSieve) {
		initBucketInfo(sieve_size, sqrtp, end - segment_low);
		initWheelBucket(minp, sqrtp, segment_low, end - segment_low);
#if RELEASE == 0
		if (Config.Flag & PRINT_LOG) {
			printf("sieve size = %d k, init time (%.2f sec)\n", sieve_size / (WHEEL30 << 10), (getTime() - ts) / 1000);
			printf("	MaxIndex = %d, LoopSize = %d\n\tStock remain = %2.2f%%, Mem = %d M\n",
					BucketInfo.MaxIndex, BucketInfo.LoopSize, StockSize * 100.0 / AllStocks,
					AllStocks * MEM_SIZE / (1 << 20));
		}
		ts = getTime();
#endif
	}

	primes += pi(start, end, sieve_size, cmd);

	if ((!cmd || cmd->Oper != FIND_MAXGAP) && (Config.Flag & PRINT_RET))
		printPiResult(start, end, primes, ts);
	if (Config.Flag & SAVE_RESUTL)
		freopen(CONSOLE, "w", stdout);
	if (BucketInfo.UseSieve) {
		for (int i = 0; i < MemSize; i++)
			free(StockMem[i]);
		assert(AllStocks == StockSize);
		MemSize = AllStocks = StockSize = 0;
	}

	return primes;
}

static int dumpPrime(double ts, uchar prime[])
{
	int j = 8 + PRIME_PRODUCT / 9699690, p = Prime[0];
	for (; prime[j + 0] > 0; j ++) {
		p += prime[j];
		PRIME_ADD_DIFF(p);
	}

	//assert(pi2(0, p + 1) == j);
	printf("\nPrime[%u] = %u, ", j, p);
	printf("save primes gap use %.2lf sec\n", ts / 1000.0);

	return j;
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

	Prime[0] = (PRIME_PRODUCT % 19 == 0) ? 23 : 19;
	uint lastp = Prime[1] = 2;

	for (uint p = 3; p <= maxp; p += 2) {
		if (0 == (bitarray[p >> 4] & (1 << (p / 2 & 7)))) {

			Prime[ ++primes] = p - lastp;
			lastp = p;
			if (p > (1 << 16)) {
				continue ;
			}

			for (uint j = p * p / 2; j <= maxp / 2; j += p)
				bitarray[j >> 3] |= 1 << (j & 7);
		}
	}

	//pack the last two byte for safety
	Prime[primes + 2] = Prime[primes + 1] = 254;

	return primes;
}

//The first presieved template
//sieve the first 8th prime multiples
static void initPreSieved( )
{
	const uchar smallprime[ ] = {7, 11, 13, 17, 19, 23, 29};

	memset(PreSieved, BITZ, sizeof(PreSieved));
	for (int i = 0; PRIME_PRODUCT % smallprime[i] == 0; i ++) {
		int p2 = 2 * smallprime[i];
		for (int sieve_index = smallprime[i]; sieve_index < PRIME_PRODUCT; sieve_index += p2) {
			PreSieved[sieve_index / WHEEL30] |= WheelInit30[sieve_index % WHEEL30].UnsetBit;
		}
	}
}

static void initBitTable( )
{
	//1. init WordNumBit1 table in 0-2^16
	int i = 0;
	int nbitsize = 1 << 16;
#if	POPCNT == 0 && TREE2 == 0
	WordNumBit1[0] = 0;
	for (i = 1; i < nbitsize; i ++)
		WordNumBit1[i] = WordNumBit1[i >> 1] + (i & 1);
#endif

#if BIT_SCANF == 0
	//2. init Left most bit table
	nbitsize = sizeof(Lsb) / sizeof(Lsb[0]);
	for (i = 0; i < nbitsize; i += 2) {
		Lsb[i + 0] = Lsb[i >> 1] + 1;
	}

	for (i = 0; i < nbitsize; i += 2) {
		Lsb[i + 0] = Pattern30[Lsb[i]];
		Lsb[i + 1] = Pattern30[0];
	}
#endif
}

static void initWheelGap()
{
	if (WheelGap == NULL) {
		WheelGap = (PrimeGap*) (malloc(sizeof(PrimeGap) * (1 << 16)));
		memset(WheelGap, 0, sizeof(PrimeGap) * (1 << 16));
	}
	if (WheelGap[3].Bits)
		return;

	for (uint k = 1; k < (1 << 16) - 1; k ++) {
		int pattern = 0, bits = 0;
		for (int j = 0; j < 16; j ++) {
			if (k & (1 << j)) {
				if (pattern == 0)
					WheelGap[k].Beg = Pattern30[j];
				else
					WheelGap[k].Gap[bits - 1] = Pattern30[j] - pattern;
				bits ++;
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

	for (int j = 0, wj = 0; j < WHEEL30; j += 1) {
		WheelInit30[j].UnsetBit = BITZ;
		WheelInit30[j].WheelIndex = wj;
		if (j == Pattern30[wj]) {
			WheelInit30[j].UnsetBit = 1 << (wj ++);
		}
	}

	for (int i = 0; i < WHEEL30; i += 1) {
		for (int pi = 0; pi < 8; pi ++) {
			int multiples = 0, sieve_index = i;
			if (i % 2 == 0) {
				multiples += 1;
				sieve_index += Pattern30[pi];
			}
			while (WheelInit30[sieve_index % WHEEL30].UnsetBit == BITZ) {
				sieve_index += Pattern30[pi] * 2;
				multiples += 2;
			}
			int wi = WheelInit30[sieve_index % WHEEL30].WheelIndex;
			WheelElement& wf = WheelFirst30[i][pi];
			wf.NextMultiple = multipleIndex[wi][pi];
			wf.WheelIndex = wi;
			wf.Correct = multiples;
			wf.UnsetBit = 1 << wi;
		}
	}

	for (int wi = 0; wi < 8; wi ++) {
		for (int pi = 0; pi < 8; pi ++) {
			int next = Pattern30[wi] + 2 * Pattern30[pi];
			int multiples = 2;
			while (WheelInit30[next % WHEEL30].UnsetBit == BITZ) {
				next += Pattern30[pi] * 2;
				multiples += 2;
			}

			WheelElement& wdata = Wheel30[pi][wi];
			wdata.Correct = next / WHEEL30 - Pattern30[wi] / WHEEL30;
			wdata.UnsetBit = WheelInit30[Pattern30[wi]].UnsetBit;
			wdata.WheelIndex = WheelInit30[next % WHEEL30].WheelIndex;
			wdata.NextMultiple = multiples * (WHEEL30 / WHEEL30);
		}
	}
}

static void initWheel210()
{
	int wi = 0, i = 0;
	const int psize = 48;
	int Pattern210[WHEEL210] = {0};

	for (i = 0; i < WHEEL210; i += 1) {
		uchar mask = WheelInit30[i % WHEEL30].UnsetBit;
		WheelInit210[i].UnsetBit = mask;
		if (mask && i % (WHEEL210 / WHEEL30)) {
			Pattern210[wi] = i;
			WheelInit210[i].WheelIndex = wi++;
		} else {
			WheelInit210[i].WheelIndex = wi;
		}
	}

	for (wi = 0; wi < psize; wi ++) {
		for (int pi = 0; pi < psize; pi ++) {
			int next = Pattern210[wi] + 2 * Pattern210[pi];
			int multiples = 2;
			while (WheelInit210[next % WHEEL210].WheelIndex ==
					WheelInit210[(next + 1) % WHEEL210].WheelIndex) {
				next += Pattern210[pi] * 2;
				multiples += 2;
			}

			WheelElement& wdata = Wheel210[pi][wi];
			wdata.Correct = next / WHEEL30 - Pattern210[wi] / WHEEL30;
			wdata.UnsetBit = WheelInit210[Pattern210[wi]].UnsetBit;
			wdata.WheelIndex = WheelInit210[next % WHEEL210].WheelIndex;
			wdata.NextMultiple = multiples * (WHEEL210 / WHEEL30);
		}
	}

	for (i = 0; i < WHEEL210; i ++) {
		for (int pi = 0; pi < psize; pi ++) {
			int multiples = 0, next = i;
			if (i % 2 == 0) {
				multiples += 1;
				next += Pattern210[pi];
			}

			while (WheelInit210[next % WHEEL210].WheelIndex ==
					WheelInit210[(next + 1) % WHEEL210].WheelIndex) {
				next += Pattern210[pi] * 2;
				multiples += 2;
			}

			WheelFirst210[i][pi].WheelIndex = WheelInit210[next % WHEEL210].WheelIndex;
			WheelFirst210[i][pi].Correct = multiples;
		}
	}
}

static void fixRangeTest(uint64 lowerBound, const uint64 range, uint64 Ret)
{
	uint llog10 = ilog10(lowerBound), rlog10 = ilog10(range);
	printf("Sieving the primes within [10^%d, 10^%d+10^%d] randomly\n", llog10, llog10, rlog10);
	uint64 primes = 0, upperBound = lowerBound + range;

	while (lowerBound < upperBound) {
		uint sieve_size = ((uint64)rand() * rand() * rand()) % 2000000000 + ipow(10, 8);
		if (sieve_size + lowerBound > upperBound)
			sieve_size = upperBound - lowerBound;

		setSieveSize(rand() % MAX_SIEVE + 128);
		primes += doSievePrime(lowerBound, lowerBound + sieve_size - 1);
		lowerBound += sieve_size;
		printf("Remaining chunk: %.2lf%%\r", (upperBound - lowerBound) * 100.0 / range);
	}
	printf("Pi[10^%d, 10^%d+10^%d] = %lld\n", llog10, llog10, rlog10, primes);
	assert(primes == Ret);
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

	uint64 primes = 0;
//	setSieveSize(CpuCache.L2Size);
	Config.Progress = 0;
	for (int i = 1; i <= 10; i ++) {
		primes = doSievePrime(0, ipow(10, i));
		printf("pi(10^%02d) = %d\n", i, primes);
		assert(primes == primeCounts[i - 1]);
	}
	printf("pi(2^32)  = %d\n", primeCounts[11]);
	assert(primeCounts[11] == doSievePrime(0, ipow(2, 32)));

	Config.Progress = 7;
	setSieveSize(MAX_SIEVE);
	for (int j = 12; j <= 19; j ++) {
		uint64 start = ipow(10, j);
		primes = doSievePrime(start, start + ipow(2, 32));
		printf("pi(10^%d, 10^%d+2^32) = %d               \n", j, j, primes);
		assert(primes == primeCounts[j]);
	}

	Config.Progress = 0;
	printf("Time elapsed %.f sec\n", (getTime() - ts) / 1000);
	puts("All Big tests passed SUCCESSFULLY!");

	uint64 RangeData[][3] =
	{
		{ipow(10, 13), ipow(10, 11), 3340141707ul},
		{ipow(10, 15), ipow(10, 11), 2895317534ul},
		{ipow(10, 17), ipow(10, 11), 2554712095ul},
		{ipow(10, 14), ipow(10, 12), 31016203073ul},
		{ipow(10, 16), ipow(10, 12), 27143405794ul},
		{ipow(10, 18), ipow(10, 12), 24127637783ul},
	};

	for (int k = 0; k < sizeof(RangeData) / sizeof(RangeData[0]); k++) {
		fixRangeTest(RangeData[k][0], RangeData[k][1], RangeData[k][2]);
	}

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
	const char* sformat3 = "Pi(%u) = %u\n";

	if (sizeof(uint64) != sizeof(uint)) {
		sformat1 = "%u PI[%llu, %llu] = %u\n";
		sformat2 = "PI[%llu, %llu] = %u\n";
		sformat3 = "PI(%llu) = %u\n";
	}

	int failedcases = 0;
	for (int i = 1; i <= tcases; i ++) {
		if (powbase) {
			uint64 start = 0 + ((uint64)rand()) * rand() * rand() % maxn;
			uint64 end = 1E6 + ((uint)rand() * rand()) % (1000000000);
			if (start > end)
				start ^= (end ^= (start ^= end));

			if (start < 10) {
				i --;
				continue ;
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
				if (failedcases ++ > 10)
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
			if ((i & 511) == 0)
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

	for (int i = 1; params[i][0] && ni < sizeof(buf) / sizeof(buf[0]); i ++) {
		if (isdigit(params[i][0]) || params[i][0] == 'e')
			buf[ni ++] = atoint64(params[i]);
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
	printf(", L2 = %d kb\n\n", cpuinfo[2] >> 16);

	if (cpuName[0] == 'A') { //amd cpu
		CpuCache.L1Size = 64 * (WHEEL30 << 10);
	} else {
		CpuCache.L1Size = 32 * (WHEEL30 << 10);
	}

	CpuCache.L2Size = L2_DCACHE_SIZE * (WHEEL30 << 10);
	CpuCache.L1Maxp = CpuCache.L1Size / WHEEL30 / ERATSMALL;
	CpuCache.L2Maxp = CpuCache.L2Size / WHEEL30 / ERATMEDIUM;

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

static void printInfo( )
{
	const char* sepator =
		"-------------------------------------------------------------------------";
	puts(sepator);
	printf("Count/Sieve primes in (0, 2^64 - 4 * 2^32), KVERSION %s\n", KVERSION);
	puts("Implemented by the segmented sieve of eratosthenes [wheel = 30/210]");
	puts("Copyright @ by Huang Yuanbing 2011 - 2013 bailuzhou@163.com");

	puts(sepator);
	puts(sepator);

#ifdef _MSC_VER
	printf("Compiled by MS/vc ++ %d", _MSC_VER);
#elif __GNUC__
	printf("Compiled by MinGW/g ++ %d.%d.%d", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#endif

#ifdef X86_64
	printf(" x86-64");
#endif

	printf(" %s %s\n", __TIME__, __DATE__);

	puts(sepator);
	printf("[MARCO] : ASM_X86, POPCNT, BIT_SCANF, TREE2 = (%d, %d, %d, %d)\n", ASM_X86, POPCNT, BIT_SCANF, TREE2);
	printf("[MARCO] : MAX_SIEVE = %dk, WHEEL_BLOCK = %dk, ERATBIG = %d, WHEEL = %d\n",
			MAX_SIEVE, WHEEL_BLOCK >> 7, ERATBIG, WHEEL);
	printf("[ARG ]  : L1Size = %dk, L2Size = %dk, SieveSize = %dk\n",
			CpuCache.L1Size / WHEEL30 >> 10, CpuCache.L2Size / WHEEL30 >> 10, Config.SieveSize / WHEEL30 >> 10);
	printf("[ARG ]  : L1Index/L1Maxp/L2Index/L2Maxp/MinPrime = (%d,%d,%d,%d,%d)\n",
		CpuCache.L1Index, CpuCache.L1Maxp, CpuCache.L2Index, CpuCache.L2Maxp, Config.MinPrime);
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
#elif defined X86_64
		"g++ -m64 -msse3 -mpopcnt %s -O3 -funroll-loops -s -pipe %s -o %s";
#else
		"g++ -m32 -msse3 -mpopcnt %s -O3 -funroll-loops -s -pipe %s -o %s";
#endif

	char commpile[257] = {0};
	if (flag == NULL)
		flag = "";
	sprintf(commpile, cxxflag, flag, __FILE__, programming);

	puts(commpile);
	system(commpile);
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
			continue ;
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
			case 'M':
				Config.Progress = (1 << cdata) - 1;
				if (Config.Progress == 0)
					Config.Progress = (1 << 30) - 1;
				break;
			case 'D':
			case 'F':
			case 'A':
				Config.Flag ^= (1 << c - 'A');
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

bool excuteCmd(const char* cmd)
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
			puts("-------------------start benchmark------------------------");
			if (isdigit(params[cmdi + 1][0])) {
				for (int i = 11; i < 20; i ++) {
					uint64 start = ipow(10, i);
					uint64 size = ipow(10, 10);
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
			Cmd cmd = {PCALL_BACK, (uchar*)printPrime, 0};
			doSievePrime(start, end, &cmd);
		} else if (cmdc == 'G') {
			//g 1425172824437699411 1476
			//http://www.ieeta.pt/~tos/gaps.html
			puts("-------------------start find max gap -------------------");
			uint64 data[4] = {0}; double ts = getTime();
			Cmd cmd = {FIND_MAXGAP, (uchar*)data, 0};
			doSievePrime(start, end, &cmd);
			printf("maxp prime gap = %d on %llu, time use = %.2f sec\n",
					(int)data[1], (data[2] - data[1]), (getTime() - ts) / 1000);
		} else if (cmdc == 'Z') {
			puts("-------------------start compile and run -----------------");
			doCompile(params[cmdi + 1]);
		} else if (cmdi >= 0 && end > 0) {
			puts("-------------------start count primes -------------------");
			if (Config.Flag & CHCECK_RESUTL)
				doSievePrime2(start, end);
			doSievePrime(start, end);
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
		initBitTable( );
		initWheel30( );
		initWheel210( );
		initPreSieved( );
		setSieveSize(sieve_size);
		initOnce = false;
	}
}

int main(int argc, char* argv[])
{
#if RELEASE
	printInfo();
//	puts(Help);
#endif

	initPrime(MAX_SIEVE);

	if (argc > 1)
		excuteCmd(argv[1]);

#if RELEASE == 0
	excuteCmd("e12 e9;");
	excuteCmd("1e16 e9");
	excuteCmd("1e19 10^10");
	excuteCmd("0-1e11 1e9");
#endif
//
//	uint64 start = (uint64)Config.MinPrime * Config.MinPrime; uint64 range = atoint64("10e9");
//	doSievePrime(start - range, start + range);

	while (true) {
		char ccmd[257];
		printf("\n>> ");
		if (!gets(ccmd) || !excuteCmd(ccmd))
			break;
	}

	return 0;
}

/***
bugs

MINGW: gcc 4.6.3
CXXFLAGS: -Ofast -msse4 -s -pipe  -march=corei7 -funroll-loops
windows 7 64 bit, AMD Althon 2 X4 641 2.8G  / Intel i3 350M 2.26G
                                 primenumber/primesieve
PI[1e11, 1e11+1e10] = 394050419        4.84 / 6.42
PI[1e12, 1e12+1e10] = 361840208        5.74 / 8.51
PI[1e13, 1e13+1e10] = 334067230        7.37 /10.50
PI[1e14, 1e14+1e10] = 310208140        9.50 /12.81
PI[1e15, 1e15+1e10] = 289531946        11.46/14.94
PI[1e16, 1e16+1e10] = 271425366        13.53/17.10
PI[1e17, 1e17+1e10] = 255481287        18.80/19.78
PI[1e18, 1e18+1e10] = 241272176        19.15/24.27
PI[1e19, 1e19+1e10] = 228568014        25.74/33.16

PI[1e18, 1e18+1e6] = 24280              0.74/1.16
PI[1e18, 1e18+1e8] = 2414886            1.52/1.80
PI[1e18, 1e18+1e9] = 24217085           4.00/4.53

sieve size = 512 k, init time (2.46 sec)
PI[1e18, 1e18+1e12] = 24127637783 (1796.57 sec)
PI[1e16, 1e16+1e12] = 27143405794 (1406.57 sec)
PI[1e14, 1e14+1e12] = 31016203073 (1006.57 sec)

vc++
	cl /O2 /Oi /Ot /Oy /GT /GL PrimeNumber.cpp
mingw/gcc
	g++ -march=native -O3 -funroll-loops -s -pipe PrimeNumber76.cpp -o PrimeNumber
	-fprofile-use -flto -fprofile-generate
else
	g++ -O3 -s -pipe PrimeNumber.cpp -o PrimeNumber -fprofile-generate/use -flto
doc:
	http://www.ieeta.pt/~tos/software/prime_sieve.html
	http://code.google.com/p/primesieve/wiki/Links
	Cache optimized linear sieve
***/
