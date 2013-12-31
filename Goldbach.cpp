/************************************************************
	this programming is the most fast algorithm for count
goldbach partition, performance and throughput improvement by minimizing cache conflicts and misses
in the last level caches of multi-cores processors.

http://www.trnicely.net
http://www.ieeta.pt/~tos/primes.html
http://code.google.com/p/primesieve/
http://numbers.computation.free.fr/Constants/Primes/twin.html.
**************************************************************/

# include <ctype.h>
# include <memory.h>
# include <stdlib.h>
# include <time.h>
# include <stdio.h>
# include <math.h>
# include <assert.h>

# define GVERSION       "2.1"
# define MAXN           "1e16"
# define MINN           10000000

# define MAX_L1SIZE     (64 << 13)
# define SEGMENT_SIZE   (510510 * 19)
# define MAX_THREADS    32

//SSE4a popcnt instruction, make sure cpu support it
#if _MSC_VER > 1400
	# define POPCNT      1
	# include <intrin.h>
#elif (__GNUC__ * 10 + __GNUC_MINOR__ > 44)
	# define POPCNT      1
	# include <popcntintrin.h>
#else
	# define POPCNT      0
#endif

//#pragma pack (16)

# if OMP
	#include <omp.h>
# endif

# define FAST_CROSS      0
# define OPT_L1CACHE     1
# define OPT_L2CACHE     1
# define PRIME_DIFF      1
# define TREE2           0
# define OPT_SEG_SIEVE   1

#if _MSC_VER && _M_AMD64
# define NO_ASM_X86      1
#endif

typedef unsigned char  uchar;
typedef unsigned short ushort;
typedef unsigned int   uint;

#ifdef _WIN32
	typedef __int64 uint64;
	#define CONSOLE "CON"
	#include <windows.h>
#else
	typedef unsigned long long uint64;
	#define CONSOLE "/dev/tty"
	#include <unistd.h>
	#include <sys/time.h>
	#include <pthread.h>
#endif

# define BSHIFT 5
# if BSHIFT == 3
	typedef uchar utype;
	# define MASK 7
# elif BSHIFT ==  4
	typedef ushort utype;
	# define MASK 15
# elif BSHIFT ==  5
	typedef uint utype;
	# define MASK 31
# endif

# define MASK_N(n)          (1 << (n & MASK))
# define SET_BIT(a, n)      a[n >> BSHIFT] |= MASK_N(n)
//# define FLP_BIT (a, n)     a[(n) >> BSHIFT] ^= MASK_N(n)
//# define CLR_BIT(a, n)      a[(n) >> BSHIFT] &= ~MASK_N(n)
//# define TST_BIT(a, n)      (a[(n) >> BSHIFT] & MASK_N(n))
# define TST_BIT2(a, n)     (a[(n / 2) >> BSHIFT] & MASK_N(n / 2))

# define CHECK_FLAG(flag)   (Config.Flag & flag)
# define SET_FLAG(flag)     Config.Flag |= flag
# define CLR_FLAG(flag)     Config.Flag &= ~flag

enum EFLAG
{
	PRINT_RET = 1 << 30,
	PRINT_TIME = 1 << ('P' - 'A'),
	PRINT_LOG = 1 << ('D' - 'A'),
	SAVE_RESUTL = 1 << ('S' - 'A'),
	CHECK_PATTERN = 1 << ('R' - 'A'),
	SAVE_TASK = 1 << ('A' - 'A'),
};

static const char* const HelpConfig = "\
	[P: Print time use]\n\
	[D: Debug log]\n\
	[S: Save result to file]\n\
	[R: Runtime check pattern]\n\
	[A: Save/Continue last task]\n\
	[M: Monitor progress m(0 - 30)]\n\
	[F: Factorial of whell prime factor f(7 - 29)]\n\
	[T: Threads number t(2 - 64)]\n\
	[C: Cpu L1/L2 data cache size c(L1:16-128, L2:128-1024)]\n";

static const char* const HelpCmd = "\
	[B: Benchmark (start) (end)]\n\
	[U: Unit test (1 - 10000)]\n\
	[N: Number of patterns (start) (count)]\n\
	[I: List (powbase) (start index) (end index)]\n\
	[L: List (start) (end/count) (step)]\n";

static const char* const HelpUse = "\
	All command/config as follow:\n\
	B, B e9 e10\n\
	C31, C1024\n\
	U, U 1000+2 2\n\
	N 2e8+20 100\n\
	I 2 10 20\n\
	P S L 1e9 2e9 1e8";

static const char* const TaskFormat =
"[Task]\n\
Wheel = %d\n\
Patterns = %d\n\
Tasks = %d\n\
Pbegi = %d\n\
Pendi = %d\n\
N = %lld\n\
Result = %lld";

static const char* const PrintFormat =
#if _MSC_VER == 1200
	"G(%I64d) = %I64d";
#else
	"G(%lld) = %lld";
#endif

/************************************/
# define PRIME_NUMS (5761455 + 160)

#if PRIME_DIFF
	static uchar Prime[PRIME_NUMS];
	#define PRIME_NEXT(p, j) p += Prime[++j]
#else
	static int Prime[PRIME_NUMS];
	#define PRIME_NEXT(p, j) p = Prime[++j]
#endif

//the smallest Moudle[i] * wheel % Prime[i] = 1
static uint Moudle[PRIME_NUMS];

typedef uchar ptype;
static ptype* Pattern = NULL;//[99532800 + 10];

#if FAST_CROSS
//CrossedTpl cross out prime <= 17
static utype CrossedTpl[(SEGMENT_SIZE >> (BSHIFT + 1)) + 100];
#endif

//bit 1 left most table
static uchar LeftMostBit1[1 << 16];

//number of bits 1 binary representation table
#if POPCNT == 0 && TREE2 == 0
static uchar WordNumBit1[1 << 16];
#endif

#if WORD_REVERSE
//table of binary representation reverse (i < 2^16 = 65536)
static ushort WordReverse[1 << 16];
#endif

//the first 10 even prime numbers
static const uchar SmallPrime[ ] =
{
	2, 3, 5, 7, 11, 13,
	17, 19, 23, 29, 31
};

//config
static struct
{
	uint Flag;
	//cpu L1/L2 size
	uint CpuL1Size;
	uint CpuL2Size;
	//work threads
	uint Threads;
	//print the progress gap
	uint PrintGap;
}
Config =
{
	PRINT_RET | PRINT_TIME, MAX_L1SIZE, 1024, 4, (1 << 9) - 1
};

static struct
{
	int Wheel;
	int Moudel;
	int SqrtN;
	int FirstIndex;

	int Patterns;
	uint64 N;
} Gp;

static struct Task
{
	int Tasks;
	int Pbegi;
	int Pendi;

	uint64 Result;
} LastTask;

static struct ThreadInfo
{
	int Pbegi;
	int Pendi;
	uint64 Result;
} TData[MAX_THREADS];

static uint64 sievePattern(int, int);

#ifdef _WIN32
static DWORD WINAPI threadProc(void* ptinfo)
#else
static void* threadProc(void* ptinfo)
#endif
{
	struct ThreadInfo* ptdata = (struct ThreadInfo*)(ptinfo);
	ptdata->Result = sievePattern(ptdata->Pbegi, ptdata->Pendi);
	return 0;
}

static void devideTaskData(int threads, int pbegi, int pendi)
{
	if (pendi > Gp.Patterns) {
		pendi = Gp.Patterns;
	}

	int tsize = (pendi - pbegi) / threads;
	tsize += tsize & 1;
	TData[0].Pbegi = pbegi;
	for (int i = 1; i < threads; i++) {
		TData[i].Pbegi = TData[i - 1].Pendi =
		TData[i - 1].Pbegi + tsize;
	}
	TData[threads - 1].Pendi = pendi;
}

static uint64 startWorkThread(int threads, int pbegi, int pendi)
{
	int i;
	//assert(pendi >= pbegi && pbegi >= 0);
	if (threads > MAX_THREADS) {
		threads = 4;
	}
	devideTaskData(threads, pbegi, pendi);

#ifdef _WIN32
	HANDLE thandle[MAX_THREADS];
	DWORD tid[MAX_THREADS];
	for (i = 0; i < threads; i++) {
		thandle[i] = CreateThread(NULL, 0, threadProc, (LPVOID)(&TData[i]), 0, &tid[i]);
	}
	for (i = 0; i < threads; i++) {
		WaitForSingleObject(thandle[i], INFINITE);
		CloseHandle(thandle[i]);
	}
#else
	pthread_t tid[MAX_THREADS];
	for (i = 0; i < threads; i++) {
		pthread_create(&tid[i], NULL, threadProc, &TData[i]);
	}
	for (i = 0; i < threads; i++) {
		pthread_join(tid[i], NULL);
	}
#endif

	uint64 gpn = 0;
	for (i = 0; i < threads; i++) {
		gpn += TData[i].Result;
	}

	return gpn;
}

//us
static double getTime()
{
#ifdef WIN32
	LARGE_INTEGER s_freq;
	LARGE_INTEGER performanceCount;
	QueryPerformanceFrequency(&s_freq);
	QueryPerformanceCounter(&performanceCount);
	return 1000 * performanceCount.QuadPart / (double)s_freq.QuadPart;
#else
	struct timeval tmVal;
	gettimeofday(&tmVal, NULL);
	return tmVal.tv_sec * 1000 + tmVal.tv_usec / 1000.;
#endif
}

//min common factor of a, b
static int gcd(int a, int b)
{
	int r = a % b;
	while (r != 0) {
		a = b;
		b = r;
		r = a % b;
	}

	return b;
}

//http://en.wikipedia.org/wiki/Extended_Euclidean_algorithm
/*
function extended_gcd(a, b)
    if a mod b = 0
        return {0, 1}
    else
        {x, y} := extended_gcd(b, a mod b)
        return {y, x-y*(a / b)}
*/
//get min Result y: ay % b = 1
//and param a, b : gcd(a, b) = 1
static int extendedEuclidean(int a, int b, int &y)
{
	int x;
	if (a == 0) {
		y = x = 0;
		if (b == 1 || b == -1) {
			y = b;
		}
	} else {
		y = extendedEuclidean(b % a, a, x);
		x -= b / a * y;
	}

	return x;
}

static uint64 ipow(const uint x, uint n)
{
	uint64 pown = 1;
	while (n --) {
		pown *= x;
	}

	return pown;
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

//convert str to uint64 ((x)*E(y)+-*(z))
//invalid input format: 123456 1234-12 e9 2e7+2^30 2e10-2 10^11-25 2e6*2
static uint64 atoint64(const char* str, uint64 defaultn = 0)
{
	uint64 n = 0;

	while (isspace(*str)) {
		str++;
	}

	if (!isdigit(*str) && !(toupper(*str) == 'E')) {
		return defaultn;
	}

	while (isdigit(*str)) {
		n = n * 10 + *str++ - '0';
	}

	if (*str && isdigit(str[1])) {
		int pows = atoi(str + 1);
		if (*str == '^') {
			n = ipow(n, pows);
		} else if (toupper(*str) == 'E') {
			if (n == 0) {
				n = 1;
			}
			n *= ipow(10, pows);
		}
		if (*str == '^' || toupper(*str) == 'E') {
			str++;
			while (isdigit(*str))
				str++;
		}
	}

	if (*str == '+') {
		n += atoint64(str + 1, 0);
	} else if (*str == '-') {
		n -= atoint64(str + 1, 0);
	} else if (*str == '*') {
		n *= atoint64(str + 1, 1);
	}

	return n;
}

//the last Qword(64 bit integer) which the leng poisition is filled with bit 1
static void packQwordBit1(utype* bitarray, const int leng)
{
	uint64* mem = (uint64*)bitarray + leng / 64;
	mem[0] |= ~(((uint64)1 << (leng % 64)) - 1);
	mem[1] = (uint64)(~0);
}

//use asm to accelerate hot spot
static void setBitArray(utype bitarray[], int s1, int s2, const int step)
{
	if (s2 > s1) {
		s2 ^= (s1 ^= (s2 ^= s1));
	}

	for (; s2 > 0; ) {
		SET_BIT(bitarray, s2); s2 += step;
		SET_BIT(bitarray, s1); s1 += step;
	}

	if (s1 > 0) {
		SET_BIT(bitarray, s1);
	}
}

static void set3BitArray(utype bitarray[], uint64& start, const int step, const int leng)
{
	int s2 = (int)(start >> 32);
	int s1 = (int)start;

	if (s2 > s1) {
		s2 ^= (s1 ^= (s2 ^= s1));
	}

	for (; s2 > 0; ) {
		SET_BIT(bitarray, s2); s2 += step;
		SET_BIT(bitarray, s1); s1 += step;
	}

	if (s1 > 0) {
		SET_BIT(bitarray, s1); s1 += step;
	}

	s1 += leng;
	s2 += leng;

	start = ((uint64)s2 << 32) | s1;
}

static void set2BitArray(utype bitarray[], uint64& start, const int step, const int leng)
{
	int s2 = (int)(start >> 32);
	int s1 = (int)start;

	if (s1 > s2) {
		s1 ^= (s2 ^= (s1 ^= s2));
	}

	for (; s2 < leng; ) {
		SET_BIT(bitarray, s2); s2 += step;
		SET_BIT(bitarray, s1); s1 += step;
	}
	if (s1 < leng) {
		SET_BIT(bitarray, s1); s1 += step;
	}

	s1 -= leng;
	s2 -= leng;

	start = ((uint64)s2 << 32) | s1;
}

//the ith bit of bitarray is map to start + 2 * i + 1
//it's difference with crossOutEvenFactor, only
//6k + 1, 6k + 5 number which is multiple of factor will be crossed out
//and gain performance 1/3 improvement
//TODO:performance issiue
static void crossOutFactor(utype bitarray[], const uint64 start, const int leng, int factor)
{
	int s1 = factor - (int)(start % factor);
	if (s1 % 2 == 0) {
		s1 += factor;
	} else if (start <= factor) {
		s1 += 2 * factor;
	}

	if (s1 > leng)
		return;

	const int bits = leng >> 1;
	if (factor < 7) {
		for (s1 >>= 1; s1 <= bits; s1 += factor) {
			SET_BIT(bitarray, s1);
		}
		return;
	}

#if 1
	const int mrid6 = ((start + s1) / factor) % 6;
	s1 >>= 1;

	int s2 = s1;
	if (mrid6 == 1) {
		s2 += factor * 2;
	} else if (mrid6 == 3) {
		s1 += factor;
		s2 += factor * 2;
	} else {
		s2 += factor;
	}

	//6k + 1, 6k + 5 number will be corssed out
	for (factor *= 3; s2 <= bits;) {
		SET_BIT(bitarray, s1); s1 += factor; //6k + 1
		SET_BIT(bitarray, s2); s2 += factor; //6k + 5
	}
	if (s1 <= bits) {
		SET_BIT(bitarray, s1);
	}
#else
	int rid30 = (s1 + start) / factor % 30;

	while (rid30 < 31 && s1 < leng) {
		SET_BIT(bitarray, s1 / 2); s1 += factor * 2;
		rid30 += 2;
	}

	s1 >>= 1;
	while (s1 + 15 * factor <= bits) {
		SET_BIT(bitarray, s1); s1 += factor * 3; //30k + 1
		SET_BIT(bitarray, s1); s1 += factor * 2; //30k + 7
		SET_BIT(bitarray, s1); s1 += factor * 1; //30k + 11
		SET_BIT(bitarray, s1); s1 += factor * 2; //30k + 13
		SET_BIT(bitarray, s1); s1 += factor * 1; //30k + 17
		SET_BIT(bitarray, s1); s1 += factor * 2; //30k + 19
		SET_BIT(bitarray, s1); s1 += factor * 3; //30k + 23
		SET_BIT(bitarray, s1); s1 += factor * 1; //30k + 29
	}

	uint m = 0x13212123;
	while (s1 <= bits) {
		SET_BIT(bitarray, s1); s1 += factor * (m & 15);
		m = (m >> 4) | (m << 28);
	}
#endif
}

//sieve multiple of each prime factor of wheelsie in [start, start + leng]
static void sieveWheelFactor(utype bitarray[], const uint64 start, const int leng, const uint wheel)
{
	assert(wheel % 6 == 0);

	for (int i = 2, p = 3; wheel % p == 0; PRIME_NEXT(p, i)) {
		crossOutFactor(bitarray, start, leng, p);
		if (start <= p) {
			assert(start == 0);
			SET_BIT(bitarray, p / 2);
		}
	}
}

//sieve prime in [start, start + leng]
static void segmentedSieve(utype bitarray[], const uint64 start, const int leng, int first)
{
	const int sqrtn = (int)(sqrt((double)start + leng)) + 1;
	for (int i = first, p = SmallPrime[first - 1]; p < sqrtn; PRIME_NEXT(p, i)) {
		crossOutFactor(bitarray, start, leng, p);
	}
	if (start == 0) {
		*(ushort*)bitarray = 0x3491;
	}
}

//reverse bit order of a byte with binary representation
//c = ((c * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
static uchar reverseByte(const uchar c)
{
	uchar n =
		(c & 0x55) << 1 | (c & 0xAA) >> 1;
	n = (n & 0x33) << 2 | (n & 0xCC) >> 2;
	n = (n & 0x0F) << 4 | (n & 0xF0) >> 4;
	return n;
}

static int reverseWord(ushort n)
{
#if WORD_REVERSE
	return WordReverse[n];
#else
//	return reverseByte(n >> 8) | (reverseByte(n) << 8);
	n = (n & 0x5555) << 1 | (n & 0xAAAA) >> 1;
	n = (n & 0x3333) << 2 | (n & 0xCCCC) >> 2;
	n = (n & 0x0F0F) << 4 | (n & 0xF0F0) >> 4;
	n = (n & 0x00FF) << 8 | (n & 0xFF00) >> 8;
	return n;
#endif
}

static inline int countBitOnes(uint64 n)
{
#if POPCNT
	//popcnt instruction : INTEL i7/SSE4.2, AMD Phonem/SSE4A
	#if _M_AMD64 || __x86_64__
	return _mm_popcnt_u64(n);
	#else
	return _mm_popcnt_u32(n) + _mm_popcnt_u32(n >> 32);
	#endif
#elif TREE2 == 0
	uint hig = (int)(n >> 32), low = (uint)n;
	return WordNumBit1[(ushort)low] + WordNumBit1[low >> 16] +
		WordNumBit1[(ushort)hig] + WordNumBit1[hig >> 16];
#else
	n -= (n >> 1) & 0x5555555555555555ull;
	n = (n & 0x3333333333333333ull) + ((n >> 2) & 0x3333333333333333ull);
	n = (n + (n >> 4)) & 0x0F0F0F0F0F0F0F0Full;
	n += n >> 8;
	n += n >> 16;
	n += n >> 32;
	return (n & 0x00000000FF);
#endif
}

//count number of bit 0 in binary representation
//!!! buffer of bitarray after position leng packeked with bit 1
static int countBitZeros(utype bitarray[], const int leng)
{
	int loops = leng >> 6;
	int bit0s = (1 + loops) << 6;

	packQwordBit1(bitarray, leng);
	for (uint64* psbuf = (uint64*) bitarray; loops >= 0; loops--) {
		bit0s -= countBitOnes(*psbuf++);
	}

	return bit0s;
}

//reverse word array bitarray with length = leng
static void reverseByteArray(ushort bitarray[], const int leng)
{
	assert(leng % 8 == 0);
	ushort* ps = bitarray;
	ushort* pe = (ushort*)((uchar*)ps + leng / 8 - 2);

	while (ps < pe) {
		ushort tmp = reverseWord(*ps);
		*ps++ = reverseWord(*pe);
		*pe-- = tmp;
	}

	if (ps == pe) {
		*ps = reverseWord(*ps);
	} else if ((uchar*)pe + 1 == (uchar*)ps) {
		*((uchar*)ps) = reverseWord(*ps) >> 8;
	}
}

//shift bitarray to low address and the offset is shiftbit bits
static int shiftBitToLow(uchar bitarray[], const int leng, const int leftbitshift)
{
	uint* plowdword = (uint*)bitarray;
	uint* phigdword = (uint*)(bitarray + 2);

	const int rightbitshift = 16 - leftbitshift;

	//bugs for -O3 and utype = uchar
	for (int i = leng / 32 + 1; i > 0; i--) {
		*plowdword++ = ((ushort)(*plowdword >> leftbitshift)) |
			(*phigdword++ << rightbitshift);
	}

	return leng / 8;
}

//reverse bit array bitarray with length = leng, bugs
static void reverseBitArray(utype bitarray[], const int leng)
{
	const int bitmod = leng % 8;
	if (bitmod == 0) {
		reverseByteArray((ushort*)bitarray, leng);
	} else {
		reverseByteArray((ushort*)bitarray, leng + 8 - bitmod);
		shiftBitToLow((uchar*)bitarray, leng + bitmod, 8 - bitmod);
	}
	packQwordBit1(bitarray, leng);
}

//make sure no divide overflow
//improvement of 100%
static inline int
asmMulDiv(const uint moudle, const uint pattern, uint p)
{
#ifdef NO_ASM_X86
	p = ((uint64)moudle) * pattern % p;
#elif !defined _MSC_VER
	__asm
	(
#if 0
		"imul %%edx\n"
		"divl %%ecx\n"
		: "=d" (p)
		: "d"(moudle), "a"(pattern), "c"(p)
#else
		"movl %1, %%eax\n"
		"imull %2\n"
		"divl %0\n"
		"movl %%edx, %0\n"
		: "+m" (p)
		: "g"(pattern), "g"(moudle)
		: "%eax", "%edx"
#endif
	);
#else
	__asm
	{
		mov eax, pattern
		imul moudle
		div p
		mov p, edx
	}
#endif

	return p;
}

static inline int
asmMulDivSub(const uint moudle, const uint pattern, uint p, const uint leng)
{
#ifdef NO_ASM_X86
	p = (((uint64)moudle) * pattern - leng) % p + leng;
#elif !(defined _MSC_VER)
	__asm
	(
#if 1
		"movl %3, %%ecx\n"
		"movl %1, %%eax\n"

		"imull %2\n"
		"subl %%ecx, %%eax\n"
		"sbbl $0, %%edx\n"

		"idivl %0\n" //no overflow!!!
		"addl %%ecx, %%edx\n"
		"movl %%edx, %0\n"
		: "+m" (p)
		: "g"(pattern), "g"(moudle), "g"(leng)
		: "%eax", "%ecx", "%edx"
#else
		"imull %3\n"

		"subl %%ecx, %%eax\n"
		"sbbl $0, %%edx\n"

		"idivl %4\n" //no overflow!!!
		"addl %%ecx, %%edx\n"
		: "=d" (p)
		: "a"(pattern), "c"(leng), "g"(moudle), "g"(p)
#endif
	);
#else
	__asm
	{
		mov eax, pattern
		mov ecx, leng

		imul moudle //moudle * pattern
		sub eax, ecx // - leng
		sbb edx, 0

		idiv p //p
		add edx, ecx // + leng
		mov p, edx
	}
#endif

	return p;
}

/****
the classic sieve of Eratosthenes implementation by bit packing
all prime less than sqrt(n) will be saved in Prime[] by difference
Prime[1] = 2, Prime[2] = 3 - 2 = 1, Prime[3] = 5 - 3, ....
Prime[i] = (ith Prime) - ((i-1) th Prime);
*/
static int simpleEratoSieve(const uint sqrtn)
{
	int primes = 1;

	utype bitarray[30000 >> (BSHIFT + 1)] = {0};
	assert(sqrtn < sizeof(bitarray) * 16);

	uint lastprime = Prime[primes++] = 2;

	for (uint p = 3; p <= sqrtn; p += 2) {
		//bit position with vlaue 0 is prime number
		if (!TST_BIT2(bitarray, p)) {
#if PRIME_DIFF
			Prime[primes++] = p - lastprime;
			lastprime = p;
#else
			Prime[primes++] = p;
#endif

			if (p > sqrtn / p) {
				continue;
			}
			for (uint j = p * p / 2; j <= sqrtn / 2; j += p) {
				SET_BIT(bitarray, j);
			}
		}
	}

	return primes;
}

//init bit tables
static void initBitTable( )
{
	//1. init WordNumBit1 table in 0-2^16, can use popcnt replace it
	int i;

#if 0 == POPCNT && 0 == TREE2
	WordNumBit1[0] = 0;
	for (i = 1; i < sizeof(WordNumBit1) / sizeof(WordNumBit1[0]); i++) {
		WordNumBit1[i] = WordNumBit1[i >> 1] + (i & 1);
	}
#endif

	//2. init bit WordReverse table
#if WORD_REVERSE
	uchar bytereverse[256] = {0};
	//reverse bit order of byte(with 8 bit) in [0, 2^8)
	for (i = 1; i < (1 << 8); i++) {
		bytereverse[i] = reverseByte((uchar)i);
	}
	//reverse bit order of short(with 16 bit) in [0, 2^16)
	for (i = 1; i < sizeof(WordReverse) / sizeof(WordReverse[0]); i++) {
		WordReverse[i] = bytereverse[i >> 8] | (bytereverse[i & 255] << 8);
	}
#endif

	//4. init LeftMostBit1 table
	for (int m = 2; m < sizeof(LeftMostBit1) / sizeof(LeftMostBit1[0]); m += 2) {
		LeftMostBit1[m + 0] = LeftMostBit1[m >> 1] + 1;
		LeftMostBit1[m + 1] = 0;
	}

#if FAST_CROSS
	//3. init CrossedTpl table, pre sieve the factor in array sievefactor
	sieveWheelFactor(CrossedTpl, 0, sizeof(CrossedTpl) * 16, SEGMENT_SIZE);
#endif
}

//
static int savePrime(const utype bitarray[], const int start, const int leng, int prime[])
{
	int primes = 0;
	for (int p = 1; p < leng; p += 2) {
		if (!TST_BIT2(bitarray, p)) {
			prime[primes++] = start + p;
		}
	}

	return primes;
}

static int savePattern(const utype bitarray[], int leng, ptype pi1pattern[])
{
	int pi1n = 0, diff = pi1pattern[0];
	int lastpattern = 0;

	for (int p = 1; p < leng; p += 2) {
		if (!TST_BIT2(bitarray, p)) {
			pi1pattern[pi1n++] = p - lastpattern;
			lastpattern = p;
//			assert(p - lastpattern < (sizeof(ptype) << 8));
		}
	}
	pi1pattern[0] += diff;
	pi1pattern[pi1n] = leng - lastpattern;

	return pi1n;
}

//min prime in range [start, n / 2]
static int getSegpattern(const int n, int start, int wheel, ptype pattern[])
{
	int gpn = 0;
	const int leng = n / 2 + ((n / 2) & 1);

	if (pattern) {
		pattern[0] = 0;
	}

	for (int sleng = SEGMENT_SIZE; start < leng; start += sleng) {
		utype bitarray[(SEGMENT_SIZE >> (BSHIFT + 1)) + 100];
		if (sleng >= leng - start) {
			sleng = leng - start;
		}

		memset(bitarray, 0, (sleng >> 4) + 1);

		sieveWheelFactor(bitarray, n - start - sleng, sleng + 16, wheel);
		reverseBitArray(bitarray, sleng >> 1);
		sieveWheelFactor(bitarray, start, sleng + 16, wheel);
		if (pattern) {
			gpn += savePattern(bitarray, sleng, pattern + gpn);
		} else {
			gpn += countBitZeros(bitarray, (sleng + 1) / 2);
		}
	}

	return gpn;
}

//optimize for memory
static int getGppattern(const uint64 n, const uint wheel, ptype pattern[])
{
	const uint moudle = (uint)(n % wheel);
	int gpn = getSegpattern(moudle, 0, wheel, pattern);
	if (pattern) {
		pattern += gpn + 1;
	}
	gpn += getSegpattern(moudle + wheel, moudle, wheel, pattern);
	if (pattern) {
		pattern[-1] = 0;
	}

	return gpn;
}

//
static uint64 sieveGpL1(utype bitarray[], const int pattern1, const int pattern2, const int leng)
{
	int k = Gp.FirstIndex;
	int p = SmallPrime[k++];
	int sleng = Config.CpuL1Size;
	int sqrtn = Gp.SqrtN;
	const int minp = sqrtn < sleng ? sqrtn : sleng;
	const int offset = leng - leng % sleng;
	//343k
	uint64 spos[43990 + 1];

	for (; p <= minp; PRIME_NEXT(p, k)) {
		const int moudle = Moudle[k];
#if OPT_SEG_SIEVE
		int s1 = asmMulDiv(moudle, pattern1, p);
		//assert(s1 * wheel + pattern1 % p == 0)
		int s2 = leng - asmMulDivSub(moudle, pattern2, p, leng);
		if (s2 < 0) s2 += p;
#else
		int s1 = leng - asmMulDiv(moudle, pattern2, p);
		int s2 = asmMulDivSub(moudle, pattern1, p, leng);
		if (s2 > leng) s2 -= p;

		s1 -= offset;
		s2 -= offset;
#endif
		spos[k] = ((uint64)s2 << 32) | s1;
	}

#if OPT_SEG_SIEVE
	for (int start = 0; start < leng; start += sleng) {
		if (start + sleng > leng) {
			sleng = leng - start;
		}
#else
	for (int start = offset; start >= 0; start -= sleng) {
#endif
		k = Gp.FirstIndex;
		p = SmallPrime[k++];

		//40%
		for (; p <= minp; PRIME_NEXT(p, k)) {
			utype* pstart = bitarray + (start >> BSHIFT);
#if OPT_SEG_SIEVE
			set2BitArray(pstart, spos[k], p, sleng);
#else
			set3BitArray(pstart, spos[k], -p, sleng);
#endif
		}
	}

	return ((uint64)k << 32) | p;
}

//set goldbach partition (k*wheel + pattern1) with pattern = pattern1
static int sieveGp(utype bitarray[], const int pattern1)
{
	int pattern2 = (Gp.Wheel + Gp.Moudel - pattern1) % Gp.Wheel;
	int leng = (int)((Gp.N - pattern1 - pattern2) / Gp.Wheel);
	int k = Gp.FirstIndex;
	int p = SmallPrime[k++];
	int sqrtn = Gp.SqrtN;

#if OPT_L1CACHE
	//performance improvement from 100 -> 72
	if (leng > Config.CpuL1Size) {
		uint64 kpos = sieveGpL1(bitarray, pattern1, pattern2, leng);
		k = (int)(kpos >> 32); p = (int)kpos;
	}
#endif

#if OPT_L2CACHE
	const uint minp = sqrtn < leng ? sqrtn : leng;
#else
	const uint minp = sqrtn;
#endif

	for (; p <= minp; PRIME_NEXT(p, k)) {
		const uint moudle = Moudle[k];
#if 1
		int s1 = leng - asmMulDiv(moudle, pattern2, p);
		int s2 = asmMulDivSub(moudle, pattern1, p, leng);
		if (s2 > leng)
			s2 -= p;
		setBitArray(bitarray, s1, s2, -p);
#else
		int s1 = asmMulDiv(moudle, pattern1, p);
		int s2 = leng - asmMulDivSub(moudle, pattern2, p, leng);
		if (s2 < 0)
			s2 += p;
		uint64 start = ((uint64)s1 << 32) | s2;
		set2BitArray(bitarray, start, p, leng);
#endif
	}

#if OPT_L2CACHE
	for (; p <= sqrtn; PRIME_NEXT(p, k)) {
		const uint moudle = Moudle[k];
#if 0
		int s1 = leng - asmMulDiv(moudle, pattern2, p);
		int s2 = asmMulDivSub(moudle, pattern1, p, leng);
		if (s2 > leng)
			s2 -= p;
		if (s2 > 0)
			SET_BIT(bitarray, s2);
		if (s1 > 0)
			SET_BIT(bitarray, s1);
#else
		int s1 = asmMulDiv(moudle, pattern1, p);
		int s2 = leng - asmMulDivSub(moudle, pattern2, p, leng);
		if (s2 < 0)
			s2 += p;
		if (s1 < leng)
			SET_BIT(bitarray, s1);
		if (s2 < leng)
			SET_BIT(bitarray, s2);
#endif
	}
#endif

	bitarray[0] |= 1;
	if (pattern1 == pattern2) {
		leng = 1 + leng / 2;
	}

	return countBitZeros(bitarray, leng);
}

//
static ptype* getNextPattern(int& pattern, ptype* pnext)
{
	ptype pdiff = *pnext++;
	pattern += pdiff;
	if (pdiff == 0) {
		pattern = Gp.Moudel + *pnext++;
	}

	return pnext;
}

static void printProgress(int tid, double time_use, int aves, int pcnt)
{
	double currper = 100.0 * pcnt / Gp.Patterns;
	double totaltime = time_use / (10 * currper);

	printf("thread[%2d]: (%.2lf%%), time ~= ", tid, currper);
	if (totaltime < 3600 * 3) {
		printf("%.2lf sec", totaltime);
	} else {
		printf("%.2lf hour", totaltime / 3600);
	}

	printf(", goldbach partition ~= %lld\n", (uint64)aves * Gp.Patterns);
}

//thread call: get result form pattern pbegi to pendi
//calcultate wheel * k + Pattern[pbegi, pendi] <= n
static uint64 sievePattern(const int pbegi, const int pendi)
{
	static int stid = 0;
	static int scnt = 0;

	if (pbegi == 0) {
		stid = scnt = 0;
	}

	int tid = ++stid;
	
	double tstart = getTime( );

	const int sieve_byte = (int)(Gp.N / Gp.Wheel / 8) + 1;
	utype* bitarray = (utype*)malloc(sieve_byte + 63);

	ptype* pnext = Pattern;
	int pattern = 0;
	for (int i = 0; i < pbegi; i++) {
		pnext = getNextPattern(pattern, pnext);
	}

	if (CHECK_FLAG(PRINT_LOG)) {
		printf("maxp / leng = %.2lf\n", Gp.SqrtN / (sieve_byte * 8.0));
	}

	uint64 gpn = 0;

	for (int pcuri = pbegi; pcuri < pendi; pcuri++) {
		pnext = getNextPattern(pattern, pnext);
#ifdef CHECK_PATTERNS
		if (CHECK_FLAG(CHECK_PATTERN)) {
			if (gcd(Gp.Wheel, pattern) != 1)
				printf("error pattern = %d\n", pattern);
			continue;
		}
#endif
		memset(bitarray, 0, sieve_byte);

		gpn += sieveGp(bitarray, pattern);

#if _M_AMD64 && _MSC_VER
		InterlockedAdd((LONG*)(&scnt), 1);
#else
		scnt++;
#endif
		if ((scnt & Config.PrintGap) == 7) {
			printProgress(tid, getTime() - tstart, (int)(gpn / (1 + pcuri - pbegi)), scnt);
		}
	}

	free(bitarray);

	if (CHECK_FLAG(PRINT_LOG)) {
		printf("thread %d: pattern[%3d - %3d] = %lld\n", tid, pbegi, pendi, gpn);
	}

	return gpn;
}

//factorial of prime factor of wheel
static int getFactorial(const uint wheel)
{
	int factorial = 1;
	for (int i = 0, p = SmallPrime[i]; wheel % p == 0; p = SmallPrime[++i]) {
		factorial *= p;
	}
	return factorial;
}

//
static int countPiPattern(uint wheel)
{
	int patterns = 1;
	for (int i = 0, p = SmallPrime[i]; wheel % p == 0; p = SmallPrime[++i]) {
		patterns *= p - 1;
		wheel /= p;
	}
	return patterns * wheel;
}

//get the frist prime index in Prime which is not a factor of wheel
static int getFirstPrime(const uint wheel)
{
	int i = 0;
	for (int p = SmallPrime[i]; wheel % p == 0; p = SmallPrime[++i]) {

	}
	return i;
}

//16 bits number p1, p2 and gcd((p1 + p2), factorial) = 1
// ---------- gp ------------------*/
static int initPattern(const uint64 n, const int wheel, ptype pattern[])
{
	double ts = getTime( );

	int pns = getGppattern(n, wheel, pattern);

	if (CHECK_FLAG(PRINT_LOG)) {
		printf("wheel patterns = %d\n", pns);
		printf("init pattern time use %.2lf ms, %.2lf fast than pi\n",
				getTime( ) - ts, 1.0 * countPiPattern(wheel) / (pns * 2));
	} else if (CHECK_FLAG(CHECK_PATTERN)) {
		printf("wheel patterns = %d\n", pns);
	}
//	assert(pns < sizeof(Pattern));

	return pns;
}

static void initStartp(const int wheel, int maxp, uint moudle[])
{
	double ts = getTime( );

	for (int i = 1, p = 2; p < maxp; PRIME_NEXT(p, i)) {
		int y;
		y = extendedEuclidean(-wheel % p, p, y);
		if (y < 0) {
			y += p;
		}
		moudle[i] = y;
		//assert(moudle[j] < p && moudle[j] >= 0);
	}

	if (CHECK_FLAG(PRINT_LOG)) {
		printf("init moudle time use %.2lf ms\n", getTime( ) - ts);
	}
}

static uint getDefaultWheel(const uint64 n)
{
	uint wheel = Gp.Wheel;
	if (wheel > 30 && wheel % 30 == 0 && n > wheel) {
		return wheel;
	}

	wheel = 1;
	for (int i = 0; SmallPrime[i] <= Gp.Wheel && SmallPrime[i] < 24;) {
		wheel *= SmallPrime[i++];
	}

	if (wheel > 30 && n > wheel) {
		return wheel;
	}

	const int powten = ilog10(n);
	if (powten <= 6) {
		wheel = 30;
	} else if (powten <= 8) {
		wheel = 210;
	} else if (powten <= 9) {
		wheel = 2310;
	} else if (powten <= 11) {
		wheel = 30030;
	} else if (powten <= 13) {
		wheel = 510510;
	} else if (powten <= 15) {
		wheel = 9699690;
	} else if (powten <= 17) {
		wheel = 223092870;
	} else {
		wheel = 223092870ul * 29;
	}

	return wheel;
}

//set sieve buffer size and adjust wheel based
//on cpu L2 cache size and n
static int getWheel(const uint64 n)
{
	uint wheel = getDefaultWheel(n);

	int cachesize = (int)(n / wheel);

	int blocks = cachesize / (Config.CpuL2Size << 13);

	wheel *= (blocks + 1);

	return wheel;
}

//get prime number with diff Result in array Prime
//segmented sieve of EratoSieve to enum prime number
static int getPrime(const int n)
{
	int pi1n = 1;
#if PRIME_DIFF
	if (Prime) {
		Prime[2] = 1 - 3;
	}
#endif

	for (int start = 0, sleng = SEGMENT_SIZE; start < n; start += sleng) {
		utype bitarray[(SEGMENT_SIZE >> (BSHIFT + 1)) + 100];
		if (sleng >= n - start) {
			sleng = n - start;
		}

#if FAST_CROSS
		memcpy(bitarray, CrossedTpl, (sleng >> 4) + 2);
		segmentedSieve(bitarray, start, sleng + 16, 8);
#else
		memset(bitarray, 0, (sleng >> 4) + 2);
		segmentedSieve(bitarray, start, sleng + 16, 2);
#endif

#if PRIME_DIFF
		pi1n += savePattern(bitarray, sleng, Prime + pi1n + 1);
#else
		pi1n += savePrime(bitarray, start, sleng, ((int*)Prime) + pi1n + 1);
#endif
	}

	assert(pi1n < sizeof(Prime) / sizeof(Prime[0]));
	if (CHECK_FLAG(PRINT_LOG)) {
		printf("pi(%d) = %d\n", n, pi1n);
	}

	return pi1n;
}

//[0 - maxn] - [n - maxn, n]
static int getPartition(const uint64 n, int maxn)
{
	int gp = 0;
	if (maxn > n / 2) {
		maxn = n / 2;
	}

	maxn += maxn & 1;

#if (OMP)
	omp_set_num_threads(Config.Threads);
	#pragma omp parallel for if (n >= MINN * 3)
#endif
	for (int start = 0; start < maxn; start += SEGMENT_SIZE) {
		utype bitarray[(SEGMENT_SIZE >> (BSHIFT + 1)) + 100];
		int sleng = SEGMENT_SIZE;
		if (sleng >= maxn - start) {
			sleng = maxn - start;
		}

		memset(bitarray, 0, (sleng >> 4) + 1);
		segmentedSieve(bitarray, start, sleng + 16, 2);
		reverseBitArray(bitarray, sleng >> 1);
		segmentedSieve(bitarray, n - start - sleng, sleng + 16, 2);

		gp += countBitZeros(bitarray, sleng >> 1);
	}

	return gp;
}

//get small partition in range[0 - min(Wheel, sqrt(n))]
//if n is less than a small fix value MIN
static int getSmallGp(const uint64 n)
{
	double ts = getTime( );

	//adjust leng for last few number
	int leng = 0;
	if (n < MINN) {
		leng = (int)(n / 2);
	} else {
		leng = Gp.SqrtN;
		if (leng < Gp.Wheel) {
			leng = Gp.Wheel;
		}
	}

	int ret = getPartition(n, leng);

	if (CHECK_FLAG(PRINT_LOG)) {
		printf("sieve small n = %lld, leng = %d", n, leng);
		printf("\nand small ret = %d, and time use %.lf ms\n", ret, getTime( ) - ts);
	}

	return ret;
}

//init Prime, Pattern, Moudle
static uint initGp(const uint64 n)
{
	if (n > Gp.N) {
		getPrime((int)sqrt((double)n) + 2001);
	}
	const uint wheel = getWheel(n);

	Gp.N = n;
	Gp.Wheel = wheel;
	Gp.Moudel = (uint)(n % wheel);
	Gp.SqrtN = (uint)sqrt((double)n + 0.1);
	Gp.FirstIndex = getFirstPrime(wheel);

	Gp.Patterns = initPattern(n, wheel, 0);
	Pattern = (ptype*)malloc(sizeof(ptype) * Gp.Patterns + 8);
	Gp.Patterns = initPattern(n, wheel, Pattern);

	initStartp(wheel, Gp.SqrtN + 256, Moudle);

	if (CHECK_FLAG(PRINT_LOG)) {
		const int factorial = getFactorial(wheel);
		printf("wheel = %d * %d, sievesize = %dk\n", factorial, wheel /factorial, (n / wheel) >> 13);
	}

	return wheel;
}

static void saveTask(struct Task& curtask)
{
	if (freopen("prime.ta", "rb", stdin)) {
		freopen(CONSOLE, "r", stdin);
		remove("prime.ta.bak");
		if (rename("prime.ta", "prime.ta.bak") != 0) {
			perror("back data fail");
		}
	}

	if (curtask.Tasks == 0) {
		curtask.Pendi += Gp.Patterns / 3;
	} else {
		curtask.Pendi += Gp.Patterns / curtask.Tasks;
	}

	if (curtask.Pendi > Gp.Patterns) {
		curtask.Pendi = Gp.Patterns;
	}

	freopen("prime.ta", "wb", stdout);
	printf(TaskFormat, Gp.Wheel, Gp.Patterns, curtask.Tasks, curtask.Pbegi, curtask.Pendi, Gp.N, curtask.Result);
	freopen(CONSOLE, "r", stdin);
	freopen(CONSOLE, "w", stdout);
}

static int parseTask(struct Task &curtask)
{
	int ret = 0;
	char linebuf[400] = {0}, taskdata[400] = {0};
	freopen("prime.ta", "rb", stdin);

	while (gets(linebuf)) {
		strcat(taskdata, strcat(linebuf, "\n"));
	}

	//read last task data
	int wheel, patterns; uint64 N;
	if (sscanf(taskdata, TaskFormat,
				&wheel, &patterns, &curtask.Tasks, &curtask.Pbegi, &curtask.Pendi, &N, &curtask.Result) != 7) {
		printf("invalid task data format: %s : %s\n", taskdata, TaskFormat);
		ret = -1;
	}

	//check task data
	if (Gp.N != N || Gp.Wheel != wheel || Gp.Patterns != patterns) {
		printf("last task \n%s is not match with current task\n", taskdata);
		ret = -2;
	}

	freopen(CONSOLE, "r", stdin);
	return ret;
}

//
static int loadTask(struct Task &curtask)
{
	int ret = 0;
	if (!freopen("prime.ta", "rb", stdin)) {
		puts("create a default task file prime.ta\n");
		saveTask(curtask);
	}

	if (parseTask(curtask) != 0) {
		curtask.Result = curtask.Pbegi = curtask.Pendi = 0;
		saveTask(curtask);
	}

	if (curtask.Tasks == 0) {
		curtask.Tasks = 4;
	}

	if (curtask.Tasks > 0) {
		curtask.Pendi = curtask.Pbegi + Gp.Patterns / curtask.Tasks + 1;
	}

	if (curtask.Pendi > Gp.Patterns) {
		curtask.Pendi = Gp.Patterns;
	}

	if (CHECK_FLAG(PRINT_LOG)) {
		printf("load last Task Data with pattern[%d - %d] ok\n", curtask.Pbegi, curtask.Pendi);
	}

	freopen(CONSOLE, "r", stdin);

	return ret;
}

//
static uint64 getGp(const uint64 n, int pn)
{
	assert(n % 2 == 0);

	if (n >= MINN) {
		double ts = getTime( );
		initGp(n);
		if (CHECK_FLAG(PRINT_LOG)) {
			printf("initGp time %.2lf ms\n", getTime() - ts);
		}
	}

	uint64 sgn = getSmallGp(n), gpn = 0;

	if (n < MINN) {
		Gp.Wheel = 0;
		return sgn;
	}

	if (pn <= 0 || pn > Gp.Patterns) {
		pn = Gp.Patterns;
	}
	int pbegi = 0;
	int pendi = pn;

	//load Last Task
	if (CHECK_FLAG(SAVE_TASK) && loadTask(LastTask) >= 0) {
		pbegi = LastTask.Pbegi;
		pendi = LastTask.Pendi;
		gpn = LastTask.Result;
	}

#if (OMP)
	omp_set_num_threads(Config.Threads);
#pragma omp parallel for reduction(+:gpn) if (n >= MINN * 3)
	for (int oi = 0; oi < Config.Threads; oi++) {
		int bi = pn / Config.Threads * oi;
		int ei = bi + pn / Config.Threads;
		if (oi == Config.Threads - 1) {
			ei = pn;
		}
		gpn += sievePattern(bi + 1, ei);
	}
#else
	if (pendi - pbegi > 10 && Config.Threads > 1) {
		gpn += startWorkThread(Config.Threads, pbegi, pendi);
	} else if (pendi > pbegi) {
		gpn += sievePattern(pbegi, pendi);
	}
#endif

	//save Current Task
	if (CHECK_FLAG(SAVE_TASK) && pbegi < pendi) {
		LastTask.Pbegi = pendi;
		LastTask.Result = gpn;
		saveTask(LastTask);
	}

	if (Pattern) {
		free(Pattern);
		Pattern = NULL;
	}

	Gp.Wheel = 0;

	return gpn + sgn;
}

static void printResult(const uint64 n, uint64 gpn, double ts)
{
	int pow10 = ilog10(n);
	if (n % ipow(10, pow10) == 0) {
		printf("G(%de%d) = %lld", (int)(n / ipow(10, pow10)), pow10, gpn);
	} else {
		printf(PrintFormat, n, gpn);
	}

	if (CHECK_FLAG(PRINT_TIME)) {
		printf(" (%.2lf sec)", (getTime() - ts) / 1000);
	}

	putchar('\n');
}

static uint64 GPartiton(const uint64 n, int pn)
{
	double ts = getTime( );

	uint64 gpn = getGp(n, pn);

	if (CHECK_FLAG(PRINT_RET)) {
		printResult(n, gpn, ts);
	}

	return gpn;
}

//test case number
static int startTest(int tesecase, bool rwflag)
{
	srand((uint)time(NULL));
	double ts = getTime( );

	if (rwflag) {
		if (!(freopen("prime.gp", "rb", stdin))) {
			printf("can not read test data file\n");
			freopen(CONSOLE, "r", stdin);
			SET_FLAG(PRINT_RET);
			return -1;
		}
		CLR_FLAG(PRINT_RET | PRINT_LOG);
		Config.PrintGap = 0;
	} else {
		if (!(freopen("prime.gp", "wb", stdout))) {
			puts("can not write test data file");
			freopen(CONSOLE, "r", stdin);
			SET_FLAG(PRINT_RET);
			return -2;
		}
		Config.PrintGap = 0;
		CLR_FLAG(PRINT_TIME);
	}

	int fails = 0;
	for (int i = 1; i <= tesecase; i++) {
		if (!rwflag) {
			uint64 n = rand( ) * rand( );
			n = (n + 4) * 2 + 0;
			if (n < 1000000) {
				n = 4 * n + 1000000;
			}
			GPartiton(n, 0);
		} else {
			uint64 res, n;
			char linebuf[256] = {0};
			gets(linebuf);
			if (sscanf(linebuf, PrintFormat, &n, &res) != 2) {
				printf("line %d is wrong data\n", i);
				if (fails++ > 30) {
					break;
				} else {
					continue;
				}
			}

			uint64 gpn = getGp(n, 0);
			if (gpn != res) {
				printf("case %d with wrong result %lld, ", i, gpn);
				printf(PrintFormat, n, res);
				putchar('\n');
			}
		}
		if ((i & 63) == 0) {
			printf("case pass %d%%\r", i * 100 / tesecase);
		}
	}

	printf("test case time use %.lf ms\n", getTime( ) - ts);

	SET_FLAG(PRINT_RET);
	freopen(CONSOLE, "w", stdout);
	freopen(CONSOLE, "r", stdin);

	return 0;
}

//list Gpk by the input Result start, end, step
static void listDiffGp(const char params[][80], int cmdi)
{
	double ts = getTime( );

	int ni = 1, step = 2;
	uint64 start = ipow(10, 9), end = start + 1000;
	uint64 buf[ ] = {0, start, end, step, 0};

	for (int i = cmdi; params[i][0] && ni < sizeof(buf) / sizeof(buf[0]); i++) {
		char c = params[i][0];
		if (isdigit(c) || toupper(c) == 'E') {
			buf[ni++] = atoint64(params[i], 10000);
		}
	}

	start = buf[1], end = buf[2], step = (int)buf[3];

	start += (start & 1);
	step += step & 1;

	if (step < 2) {
		step = 2;
	}
	if (start > end) {
		end = start + end * step - 1;
	}

	if (CHECK_FLAG(SAVE_RESUTL)) {
		Config.PrintGap = 0;
		CLR_FLAG(PRINT_TIME);
		freopen("batch.txt", "wb", stdout);
	}

	printf("%lld:%d:%d\n\n", start, (int)(2 + end - start) / 2, step);

	int pcnt = 0;
	uint64 allSum = 0;
	for (uint64 n = start; n <= end; n += step) {
		pcnt++;
		if (isdigit(params[4][0])) {
			printf("%d ", pcnt);
		}
		allSum += GPartiton(n, 0);
	}

	printf("average = %lld, ", allSum / pcnt);

	printf("all case time use %.lf ms\n", getTime( ) - ts);
	freopen(CONSOLE, "w", stdout);
}

static void listPowGp(const char params[][80], int cmdi)
{
	int m = atoint64(params[cmdi + 1], 10);
	int startindex = atoint64(params[cmdi + 2], 5);
	int endindex = atoint64(params[cmdi + 3], 10);

	printf("in %d^%d - %d^%d\n", m, startindex, m, endindex);

	if (m < 2 || m > 10000) {
		m = 10;
	}
	if (startindex > endindex) {
		startindex ^= (endindex ^= (startindex ^= endindex));
	}

	if (CHECK_FLAG(SAVE_RESUTL)) {
		Config.PrintGap = 0;
		CLR_FLAG(PRINT_TIME);
		freopen("batch.txt", "wb", stdout);
	}

	CLR_FLAG(PRINT_RET);
	for (uint64 i = startindex; i <= endindex; i++) {
		uint64 n = ipow(m, i);
		uint64 r = getGp(n, 0);
		printf(PrintFormat, i, r);
		putchar('\n');
	}
	SET_FLAG(PRINT_RET);
}

static void listPatterns(uint64 start, int count)
{
	Gp.Wheel = 0;
	getPrime((int)sqrt((double)start) + count * 4 + 2000);
	uint wheel = getWheel(start);
	printf("wheel = %u\n", wheel);
	for (int i = 0; i < count; i++) {
		int patterns = getGppattern(start, wheel, 0);
		printf("pattern(%lld) = %d\n", start, patterns);
		start += 2;
	}
}

//benchMark
static void benchMark(const char params[][80])
{
	SET_FLAG(PRINT_RET);
	uint64 start = atoint64("e10", 1000000);
	uint64 n = atoint64("e10", 0);

	for (int i = 0, ni = 1; params[i][0]; i++) {
		char c = params[i][0];
		if (isdigit(c) || toupper(c) == 'E') {
			uint64 tmp = atoint64(params[i], 0);
			if (ni++ == 1) {
				start = tmp;
			} else if (ni == 3) {
				n = tmp;
			}
		}
	}

	if (CHECK_FLAG(SAVE_RESUTL)) {
		freopen("benchmark.txt", "wb", stdout);
	}

	uint64 gap = start;
	for (; start * 10 <= n; ) {
		for (int j = 0; j < 9; j++) {
			GPartiton((j + 1) * gap, 0);
		}
		start *= 10;
		gap *= 10;
	}

	freopen(CONSOLE, "w", stdout);
}

//test pi, pi2, gp function
static void testGp( )
{
	const char* gpdata[][6] =
	{
		{"1e07", "7",  "000", "00038807"},
		{"1e08", "7",  "000", "00291400"},
		{"1e09", "11", "000", "02274205"},
		{"1e10", "13", "000", "18200488"},
		{"1e11", "13", "100", "15059883"},
		{"1e12", "17", "200", "16756829"},
		{"1e13", "17", "060", "14200452"},
		{"1e14", "19", "080", "14385672"},
		{"1e15", "19", "060", "14693072"},
		{"1e16", "23", "060", "16146185"}
	};

	Config.CpuL2Size = 1 << 10;
	Config.PrintGap = 0;
	SET_FLAG(PRINT_TIME);

	for (int i = 0; i < sizeof(gpdata) / sizeof(gpdata[0]); i++) {
		uint64 n = atoint64(gpdata[i][0], 0);
		int psize = atoi(gpdata[i][2]);
		Gp.Wheel = atoi(gpdata[i][0 + 1]);
		int gp = getGp(n, psize);
		if (atoi(gpdata[i][3]) != gp) {
			printf("%s : %s ~= %d fail !!!\n", gpdata[i][0], gpdata[i][3], gp);
		}
	}
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
		Config.CpuL1Size = 64 << 13;
	} else {
		Config.CpuL1Size = 32 << 13;
	}

//	Config.CpuL2Size = 256 << 13;

	return cpuinfo[2] >> 16;
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
#else
	Config.Threads = sysconf(_SC_NPROCESSORS_CONF);
#endif

	return Config.Threads;
}

//print the Gpk info
static void printInfo( )
{
	const char* const sepator =
		"-------------------------------------------------------------------------";
	puts(sepator);
	puts(sepator);

	printf("Goldbach partition (n < %s) version %s\n", MAXN, GVERSION);
	puts("Copyright (c) by Huang Yuanbing 2008 - 2014 bailuzhou@163.com");

#ifdef _MSC_VER
	printf("Compiled by MS/vc++ %d", _MSC_VER);
#else
	printf("Compiled by g++ %d.%d.%d", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#endif

#if _M_AMD64 || __x86_64__
	printf(" on x64");
#endif
	printf(" on %s %s\n", __TIME__, __DATE__);

	printf("[MARCO] Work threads = %d, POPCNT = %d, PRIME_DIFF = %d\n", Config.Threads, POPCNT, PRIME_DIFF);
	printf("[MARCO] L1 = %d k, OPT_L1/L2 = %d/%d, BSHIFT = %d\n", Config.CpuL1Size >> 13, OPT_L1CACHE, OPT_L2CACHE, BSHIFT);
	puts(sepator);
	puts(sepator);
}

static void doCompile()
{
	char exename[64];
	strcpy(exename, __FILE__);
	char* pdot = strchr(exename, '.');
	if (pdot) {
		strcpy(pdot, "_.exe");
		puts(exename);
	}

	const char* const cxxflag =
#ifdef _MSC_VER
		"cl /O2 /Oi /Ot /Oy /GT /GL %s %s";
#else
		"g++ -mpopcnt -mtune=native -O3 -s -pipe -fomit-frame-pointer -lpthread %s -o %s";
#endif

	char compileLine[256] = {0};
	sprintf(compileLine, cxxflag, __FILE__, exename);
	puts(compileLine);
	system(compileLine);
}

//
static int parseCmd(char params[][80])
{
	int cmdi = -1;

	for (int i = 0; params[i][0]; i++) {
		char c = params[i][0];
		int tmp = atoi(params[i] + 1);
		if (c >= 'a' && c <= 'z') {
			c += 'A' - 'a';
		}
		if (isdigit(c) || c == 'E') {
			if (cmdi < 0) {
				cmdi = i;
			}
			continue;
		}

		switch (c)
		{
			case 'A':
			case 'D':
			case 'P':
			case 'S':
			case 'R':
				Config.Flag ^= (1 << (c - 'A'));
				break;
			case 'C':
				if (tmp <= (MAX_L1SIZE >> 13) && tmp > 15) {
					Config.CpuL1Size = tmp << 13;
				} else if (tmp < 4000) {
					Config.CpuL2Size = tmp;
				}
				break;
			case 'F':
				if (tmp > 0 && tmp < 214748364) {
					Gp.Wheel = tmp;
				}
				break;
			case 'T':
				if (tmp < MAX_THREADS && tmp > 0) {
					Config.Threads = tmp;
				}
				break;
			case 'M':
				if (tmp >= 2 && tmp <= 30) {
					Config.PrintGap = (1 << tmp) - 1;
				} else if (tmp == 1) {
					Config.PrintGap = Config.PrintGap * 2 + 1;
				} else if (tmp == 0) {
					Config.PrintGap >>= 1;
				}
				break;
			case 'H':
				puts(HelpConfig);
				puts(HelpCmd);
				puts(HelpUse);
				break;
			default:
				cmdi = i;
				break;
		}
	}

	return cmdi;
}

//split cmd to params[] by white space and ';'
static int splitCmd(const char* ccmd, char params[][80])
{
	int ncmds = 0;

	for (int i = 0; ; i++) {
		while (isspace(*ccmd)) {
			ccmd++;
		}
		if (*ccmd == 0 || *ccmd == ';') {
			break;
		}
		char* pc = params[i];
		char c = *ccmd;
		bool isValid = false;
		//only space alnum and "^+-*;" is valid input
		while (isalnum(c) || c == '^' ||
				c == '+' || c == '-' || c == '*') {
			*pc++ = c;
			c = *++ccmd;
			isValid = true;
		}
		if (isValid)
			ncmds++;
		else
			ccmd++;
	}

	return ncmds;
}

//
static bool excuteCmd(const char* cmd)
{
	while (cmd) {

		// split each command by ';'
		char* pcmd = (char*) strchr(cmd, ';');
		char params[8][80] = {0};

		if (splitCmd(cmd, params) <= 0) {
			return false;
		}

		int cmdi = parseCmd(params);
		if (cmdi < 0) {
			return true;
		}

		char cmdc = toupper(params[cmdi][0]);
		uint64 n1 = atoint64(params[cmdi], 10000000000);
		uint64 n2 = atoint64(params[cmdi + 1]);
		if (!isdigit(cmdc) && cmdc != 'E') {
			n1 = n2;
			n2 = atoint64(params[cmdi + 2]);
		}

		if (cmdc == 'B') {
			benchMark(params);
		} else if (cmdc == 'U') {
			testGp( ); startTest(n1, n2 == 0);
		} else if (cmdc == 'L') {
			listDiffGp(params, cmdi);
		} else if (cmdc == 'I') {
			listPowGp(params, cmdi);
		} else if (cmdc == 'N') {
			listPatterns(n1, n2);
		} else if (cmdc == 'E' || isdigit(cmdc)) {
			GPartiton(n1, n2);
		}

		if (pcmd) {
			cmd = pcmd + 1;
		} else {
			break;
		}
	}

	return true;
}

//
static void initCache( )
{
	simpleEratoSieve(10000);
	initBitTable( );
	getSystemInfo();
	getCpuInfo();
}

int main(int argc, char* argv[])
{
	initCache( );

	if (argc < 2) {
		printInfo( );
	}

	for (int i = 1; i < argc; i++) {
		if (argv[i][0] == 'm')
			doCompile();
		else
			excuteCmd(argv[i]);
	}

	excuteCmd("1e10;1e11");
	excuteCmd("t1 d m7 1e14 200; e15 100");
//	excuteCmd("d m7 1e15");

	char ccmd[256] = {0};
	while (true) {
		printf("\n>> ");
		if (!gets(ccmd) || !excuteCmd(ccmd))
			break;
	}

	return 0;
}

/*************************************************************************
G(1e14) = 90350630388, 6481.20 s AMD X4 820

MINGW: gcc 4.6.3
CXXFLAGS: -Ofast -msse4 -s -pipe  -march=corei7 -funroll-loops
windows 7 64 bit, AMD Althon 2 X4 641 2.8G  / Intel i3 350M 2.26G
G(1e11) = 149091160        3.54  | 09.0 sec
G(1e12) = 1243722370       42.8  | 89.4 sec
G(1e13) = 10533150855      491.  | 990. sec
G(1e14) = 90350630388      7323  | 3.20 hour
G(1e15) = 783538341852     36.0  | 45.4 hour
G(1e16) =                  400   |

  c1600 m5 d t4 e15 1000
	need 28h amd phoenm x4 830
feature:
  1. win32 gui
  2. detial algorithm and readme

Linux g++:
  g++ -Wall -march=native -s -pipe -O2 -ffunction-sections -fomit-frame-pointer -lpthread GPartiton.cpp
Mingw/g++:
  g++ -Wall -mpopcnt -mtune=native -O2 -s -pipe GPartiton.cpp
MS vc++:
  cl /O2 /Os GPartiton.cpp

 ****************************************************/

