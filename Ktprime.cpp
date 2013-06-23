/************************************************************
copyright (C) 2008-2013 by Huang Yuanbing
mail to: bailuzhou@163.com
free use for non-commercial purposes


http://primes.utm.edu/glossary/xpage/PrimeKTuplet.html
Prime k-tuplet definition at the Prime Glossary.
http://mathworld.wolfram.com/PrimeConstellation.html

http://anthony.d.forbes.googlepages.com/ktuplets.htm
K-tuplet definition and records, maintained by Tony Forbes.

http://www.trnicely.net
Computational prime research page of Thomas R. Nicely. Tables of values of pi(x), pi2(x), pi3(x) and pi4(x).

http://www.ieeta.pt/~tos/primes.html
Computational prime research page of Tomás Oliveira e Silva. Tables of values of pi(x) and of pi2(x).

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

# define KVERSION       "9.5"
# define TABLE_GAP      "1e11"
# define MAXN           "1e16"
# define MINN           10000000

# define MAX_L1SIZE     (64 << 13)
# define MAX_L2SIZE     (256 << 13)
# define SEGMENT_SIZE   (510510 * 19)
# define MAX_THREADS    32
# define POPCNT         0

//SSE4 popcnt instruction, make sure your cpu support it
//use of the SSE4.2/ SSE4a POPCNT instruction for fast bit counting.
#if _MSC_VER > 1400
	# include <intrin.h>

	# pragma warning(disable: 4996 4244 4127 4505 4018)
	#if _MSC_VER > 1200
	# pragma warning (disable:6328 6031)
	#endif

#elif (__GNUC__ * 10 + __GNUC_MINOR__ > 44)
//	# include <popcntintrin.h>
#endif

# define TREE2           0

# define OMP             0
# if OMP
	#include <omp.h>
# endif

# define FAST_CROSS      1
# define OPT_L1CACHE     0
# define OPT_L2CACHE     1

#if defined _M_AMD64
	# define ASM_X86     0
#if _MSC_VER
	# define LINE_COVER  1
#endif
#elif _MSC_VER >= 1200
	# define ASM_X86     1
#else
	# define ASM_X86     0
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

# define MASK_N(n)          (1 << ((n) & MASK))
# define SET_BIT(a, n)      a[(n) >> BSHIFT] |= MASK_N(n)
//# define FLP_BIT (a, n)     a[(n) >> BSHIFT] ^= MASK_N(n)
//# define CLR_BIT(a, n)      a[(n) >> BSHIFT] &= ~MASK_N(n)
# define TST_BIT(a, n)      (a[(n) >> BSHIFT] & MASK_N(n))
# define TST_BIT2(a, n)     TST_BIT(a, (n) / 2)

# define CHECK_FLAG(flag)   (Config.Flag & flag)
# define SET_FLAG(flag)     Config.Flag |= flag
# define CLR_FLAG(flag)     Config.Flag &= ~(flag)

static const char* const HelpConfig = "\
	[P: Print time use]\n\
	[D: Debug log]\n\
	[S: Save result to file]\n\
	[R: Runtime check pattern]\n\
	[A: Save/Continue last task]\n\
	[K: Calculate of Prime/Twin[Cousin]/Ktuplet Prime k(0 - 8)]\n\
	[M: Monitor progress m(0 - 30)]\n\
	[H: Help]\n\
	[F: Factorial of whell prime factor f(7 - 29)]\n\
	[T: Threads number t(2 - 64)]\n\
	[C: Cpu L1/L2 data cache size (L1:16-128, L2:128-1024)]\n";

static const char* const HelpCmd = "\
	[B: Benchmark (start) (end) (gptk)]\n\
	[U: Unit test (n 1 - 10000) (gptk 0 - 2)]\n\
	[N: Number of patterns (start)\n\
	[I: List base pow index (powbase) (start) (end)]\n\
	[L: List multi gptk (start) (end/count) (step)]\n";

static const char* const HelpUse = "\n\
	All command/config as follow:\n\
	B, B e9 e10 123\n\
	C31, C128000\n\
	K41, K52, KK26, KK268 T2-32\n\
	H, Hk, A, D, S, R\n\
	U, U 1000 012, U 1000+2 2\n\
	N 120000*1000\n\
	I 2 10 20\n\
	L 2e9-100 1000 10\n\
	L e9-100 2e9*2 1e8+2";

static const char* const Kpattern = "\n\
	k11 0\n\
	k21 2\n\
	k31 2 6\n\
	k32 4 6\n\
	k41 2 6  8\n\
	k51 2 6  8  12\n\
	k52 4 6  10 12\n\
	k61 4 6  10 12 16\n\
	k71 2 6  8  12 18 20\n\
	k72 2 8  12 14 18 20\n\
	k81 2 6  8  12 18 20 26\n\
	k82 2 6  12 14 20 24 26\n\
	k83 6 8  14 18 20 24 26\n\
	k91 2 6  8  12 18 20 26 30\n\
	k92 2 6  12 14 20 24 26 30\n\
	k93 4 6  10 16 18 24 28 30\n\
	k94 4 10 12 18 22 24 28 30";

static const char* const KtupletName[] =
{
	"Prime numbers",
	"Twin primes",
	"Prime triplets",
	"Prime quadruplets",
	"Prime quintuplets",
	"Prime sextuplets",
	"Prime septuplets",
	"Prime octuplets",
	"Prime nonuplets"
};

static const char* const TaskFormat =
"[Task]\n\
Ptk = %d\n\
Kgap = %d\n\
Wheel = %d\n\
Patterns = %d\n\
Tasks = %d\n\
Pbegi = %d\n\
Pendi = %d\n\
N = %lld\n\
Result = %lld";

static const char* const PrintFormat[] =
{
#if _MSC_VER == 1200
	"PI(%I64d) = %I64d",
	"PI2(%I64d) = %I64d",
	"PI%d(%I64d) = %I64d"
#else
	"PI(%lld) = %lld",
	"PI2(%lld) = %lld",
	"PI%d(%lld) = %lld"
#endif
};

enum EFLAG
{
	PRINT_RET = 1 << 30,
	PRINT_TIME = 1 << ('P' - 'A'),
	PRINT_LOG = 1 << ('D' - 'A'),
	SAVE_RESUTL = 1 << ('F' - 'A'),
	CHCECK_PATTERN = 1 << ('R' - 'A'),
	SAVE_TASK = 1 << ('A' - 'A'),
};

enum CAMODE
{
	PI1N = 0, //prime number
	PI2N = 1, //twin prime
	PIKN = 2, //ktuplet prime
};

/************************************/
# define PRIME_NUMS (5761455 + 160)
//prime difference in [0, 10^8]
#define PRIME_DIFF 1
#if PRIME_DIFF
	static uchar Prime[PRIME_NUMS];
	#define PRIME_NEXT(p, j) p += Prime[++j]
#else
	static int Prime[PRIME_NUMS];
	#define PRIME_NEXT(p, j) p += Prime[++j]
#endif

//the smallest Moudle[i] * wheel % Prime[i] = 1
static uint Moudle[PRIME_NUMS];

typedef uchar ptype;
static ptype* Pattern = NULL; //[1658880 * 22 + 10];

//table of ktuplet
static uint64 Ktable[10000];

//CrossedTpl cross out prime <= 17
static utype CrossedTpl[(SEGMENT_SIZE >> (BSHIFT + 1)) + 100];

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

	//flag for PI2(n), PI(n), Pik(n)
	int Ptk;
	//ktuplet pattern gap, 2 for twin
	int Kgap;
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
	PRINT_RET | PRINT_TIME, PI2N, 2, MAX_L1SIZE, 1300, 4, (1 << 9) - 1
};

static struct
{
	int Wheel;

	int Factorial;
	int SqrtN;
	int Patterns;
	int firstIndex;
	int PatternDiff;

	bool UseKtable;
	int Kpattern[16];
	uint64 N;
	uint64 S;
}
KData =
{
	0, 8, 0, 0,
	0, 0, 0, {4, 2, 6, 8},
};

static struct Task
{
	int Ptk;
	int Wheel;
//	int Patterns;

	int Tasks;

	int Pbegi;
	int Pendi;

	uint64 N;
	uint64 S;
	uint64 Result;
} LastTask;

static struct ThreadInfo
{
	int Pbegi;
	int Pendi;
	uint64 Result;
} TData[MAX_THREADS];

static bool excuteCmd(const char* cmd);
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
	if (pendi > KData.Patterns) {
		pendi = KData.Patterns;
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
	uint64 gpts = 0;
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
		thandle[i] = CreateThread(NULL, 0, threadProc,
			(LPVOID)(&TData[i]), 0, &tid[i]);
		if (thandle[i] == NULL) {
			printf("create win32 thread error %ld\n", GetLastError());
		}
	}
	for (i = 0; i < threads; i++) {
		WaitForSingleObject(thandle[i], INFINITE);
		CloseHandle(thandle[i]);
	}
#else
	pthread_t tid[MAX_THREADS];
	for (i = 0; i < threads; i++) {
		int error = pthread_create(&tid[i], NULL, threadProc, &TData[i]);
		if (error != 0) {
			printf("create posix thread error %d\n", error);
		}
	}
	for (i = 0; i < threads; i++) {
		pthread_join(tid[i], NULL);
	}
#endif

	for (i = 0; i < threads; i++) {
		gpts += TData[i].Result;
	}

	return gpts;
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
        x, y := extended_gcd(b, a mod b)
        return y, x-y*(a / b)
*/
//get min Result y: ay % b = 1
//and param a, b : gcd(a, b) = 1
static int extendedEuclid(int a, int b, int &y)
{
	int x;
	if (a == 0) {
		y = x = 0;
		if (b == 1 || b == -1) {
			y = b;
		}
	} else {
		y = extendedEuclid(b % a, a, x);
		x -= b / a * y;
	}

	return x;
}

static uint64 ipow(uint64 x, uint n)
{
	uint64 result = 1;
	while (n != 0) {
		if (n & 1) {
			result *= x;
			n -= 1;
		}
		x *= x;
		n >>= 1;
	}

	return result;
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

//the last Qword(64 bit integer) which the bitleng poisition is filled with bit 1
static void packQwordBit1(utype* bitarray, const int bitleng)
{
	uint64* memstart = (uint64*)bitarray + (bitleng >> 6);
	memstart[0] |= ~(((uint64)1 << (bitleng % (1 << 6))) - 1);
	memstart[1] = (uint64)(~0);
}

//use asm to accelerate hot spot
#if _MSC_VER && ASM_X86
#define EBP_OFFSET	32
__declspec(naked)
#endif
inline static void
set2BitArray(utype bitarray[], int s1, int s2, const int step)
{
#if ASM_X86 == 0 || _MSC_VER == 0
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
#else
	__asm
	{
		pushad //store all register into stack, the esp will decrease 0x20

		mov ebp, dword ptr [esp + EBP_OFFSET +  4] //bitarray
		mov eax, dword ptr [esp + EBP_OFFSET +  8] //s2
		mov esi, dword ptr [esp + EBP_OFFSET + 12] //s1
		mov edx, dword ptr [esp + EBP_OFFSET + 16] //step
		cmp esi, eax

		jge $zero
		xchg eax, esi

$zero:
		cmp eax, 0
		jle $s1

$loop:
		mov edi, esi
		mov ecx, esi
		shr edi, BSHIFT
		mov ebx, 1
		and ecx, MASK
		shl ebx, cl

#if   BSHIFT == 3
		or byte ptr[ebp + edi], bl
#elif BSHIFT == 4
		or word ptr[ebp + 2 * edi], bx
#else
		or dword ptr[ebp + 4 * edi], ebx
#endif

		add esi, edx

		mov edi, eax
		mov ecx, eax
		shr edi, BSHIFT
		mov ebx, 1
		and ecx, MASK
		shl ebx, cl

#if   BSHIFT == 3
		or byte ptr[ebp + edi], bl
#elif BSHIFT == 4
		or word ptr[ebp + 2 * edi], bx
#else
		or dword ptr[ebp + 4 * edi], ebx
#endif
		add eax, edx
		jg $loop

$s1:
		cmp esi, 0
		jl $end

		mov edi, esi
		mov ecx, esi
		shr edi, BSHIFT
		mov ebx, 1
		and ecx, MASK
		shl ebx, cl

#if   BSHIFT == 3
		or byte ptr[ebp + edi], bl
#elif BSHIFT == 4
		or word ptr[ebp + 2 * edi], bx
#else
		or dword ptr[ebp + 4 * edi], ebx
#endif

$end:
		popad
		ret
	}
#endif
}

static void inline
setBitArray(utype bitarray[], int start, const int step)
{
#if (ASM_X86 == 0)
	for (; start > 0; start += step) {
		SET_BIT(bitarray, start);
	}
#elif defined _MSC_VER
	__asm
	{
#if 1
		mov esi, bitarray
		mov eax, start
		mov edx, step
$loop1:
		mov edi, eax
		mov ecx, eax
		shr edi, BSHIFT
		mov ebx, 1
		and ecx, MASK
		shl ebx, cl

	#if BSHIFT == 3
		or byte ptr[esi + edi], bl
	#elif BSHIFT == 4
		or word ptr[esi + 2 * edi], bx
	#else
		or dword ptr[esi + 4 * edi], ebx
	#endif

		add eax, edx
		jg $loop1
#else
		mov ecx, [bitarray]
		mov eax, start
		mov edx, step
$set:
		bts [ecx], eax
		add eax, edx
		jg $set
#endif
	}
#else
	__asm
	(
		"leal (%1), %%ecx\n"
		"movl %2, %%eax\n"
		"movl %3, %%edx\n"
"set:\n"
		"btsl %%eax, (%%ecx)\n"
		"addl %%edx, %%eax\n"
		"cmpl %%eax, $0\n"
		"jl set\n"
		: "=m" (step)
		: "r" (bitarray), "r" (start), "g" (step)
		: "ax", "memory", "%eax", "%ecx", "%edx"
	);
#endif
}

static uint64
set2BitArray(utype bitarray[], const uint64 start, const int step, const int bitleng)
{
	int s2 = (int)(start >> 32);
	int s1 = (int)start;

	if (s1 > s2) {
		s1 ^= (s2 ^= (s1 ^= s2));
	}

	for (; s2 < bitleng; ) {
		SET_BIT(bitarray, s2); s2 += step;
		SET_BIT(bitarray, s1); s1 += step;
	}
	if (s1 < bitleng) {
		SET_BIT(bitarray, s1); s1 += step;
	}

	s1 -= bitleng;
	s2 -= bitleng;

	return ((uint64)s2 << 32) | s1;
}

static uint64
setBitArray0(utype bitarray[], const uint64 start, const int step, const int leng)
{
	int s2 = start >> 32;
	int s1 = (int)start;

	if (s2 > s1) {
		s1 ^= (s2 ^= (s1 ^= s2));
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
//	assert(s1 >= 0 && s2 >= 0);

	return ((uint64)s2 << 32) | s1;
}

static inline int
setBitArray(utype bitarray[], int start, int step, int bitleng)
{
#if (ASM_X86 == 0)
	for (; start < bitleng; start += step) {
		SET_BIT(bitarray, start);
	}
#elif defined _MSC_VER
	__asm
	{
#if 1
		mov edx, step
		mov esi, bitarray
		mov eax, start
$loop1:
		mov edi, eax
		mov ecx, eax
		shr edi, BSHIFT
		mov ebx, 1
		and ecx, MASK
		shl ebx, cl

#if   BSHIFT == 3
		or byte ptr[esi + edi], bl
#elif BSHIFT == 4
		or word ptr[esi + 2 * edi], bx
#else
		or dword ptr[esi + 4 * edi], ebx
#endif

		add eax, edx
		cmp eax, bitleng
		jl $loop1
		mov start, eax
#else
		mov ebx, [bitarray]
		mov esi, bitleng
		mov eax, start
		mov edi, step
$loop2:
		bts [ebx], eax
		add eax, edi
		cmp eax, esi
		jl $loop2
		mov start, eax
#endif
	}
#else
	__asm
	(
		"leal (%1), %%eax\n"
		"movl %3, %%esi\n"
		"movl %2, %%ebx\n"
		"movl %4, %%edi\n"
		"movl %4, %%ecx\n"
"$loop1:\n"
		"btl %%ebx, (%%eax)\n"
		"jz set\n"
		"btsl %%ebx, (%%eax)\n"
"set:\n"

		"btsl %%ebx, (%%eax)\n"
		"addl %%edi, %%ebx\n"
		"cmpl %%esi, %%ebx\n"
		"jl $loop1\n"
		"movl %%eax, %3\n"
		: "=m"(start)
		: "r" (bitarray), "g" (start), "g" (bitleng),"g" (step)
		: "ax", "memory", "%eax", "%ebx", "%ecx", "%edx", "%edi"
	);
#endif

	return start;
}

//the ith bit of bitarray is map to start + 2 * i + 1
//it's difference with crossOutEvenFactor, only
//6k + 1, 6k + 5 number which is multiple of factor will be crossed out
//and gain performance 1/3 improvement
static void crossOutFactor(utype bitarray[], const uint64 start, const int leng, int factor)
{
	int s1 = factor - (int)(start % factor);
	if (s1 % 2 == 0) {
		s1 += factor;
	} else if (start <= factor) {
		s1 += 2 * factor;
	}

	if (s1 > leng) {
		return ;
	}

	const int bits = leng >> 1;
	if (factor < 7) {
		for (s1 >>= 1; s1 <= bits; s1 += factor) {
			SET_BIT(bitarray, s1);
		}
		return ;
	}


#if 0
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
			SET_BIT(bitarray, (p - start) / 2);
		}
	}
}

//sieve prime in [start, start + leng]
static void segmentedSieve1(utype bitarray[], const uint64 start, const int leng)
{
	const int sqrtn = (int)(sqrt((double)start + leng) + 0.1) + 1;
	memset(bitarray, 0, (leng >> 4) + 1);

	for (int i = 2, p = 3; p < sqrtn; PRIME_NEXT(p, i)) {
		crossOutFactor(bitarray, start, leng, p);
	}

	if (start == 0) {
		*(ushort*)bitarray = 0x3491;
	}
}

//another algorithm for sieve prime in [start, start + leng]
static void segmentedSieve2(utype bitarray[], const uint64 start, const int leng)
{
	assert(start % SEGMENT_SIZE == 0 && SEGMENT_SIZE % 19);
	const int sqrtn = (int)(sqrt((double)start + leng) + 0.1) + 1;
	memcpy(bitarray, CrossedTpl, (leng >> 4) + 1);

	for (int i = 8, p = 19; p < sqrtn; PRIME_NEXT(p, i)) {
		crossOutFactor(bitarray, start, leng, p);
	}

	if (start == 0) {
		//the first 7th bit poisition set 1
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
//!!! buffer of bitarray after position bitleng packeked with bit 1
static int countBitZeros(utype bitarray[], const int bitleng)
{
	int loops = bitleng >> 6;
	int bit0s = (1 + loops) << 6;

	packQwordBit1(bitarray, bitleng);
	for (uint64* psbuf = (uint64*) bitarray; loops >= 0; loops--) {
		bit0s -= countBitOnes(*psbuf++);
	}

	return bit0s;
}

//reverse word array bitarray with length = bitleng
static void reverseByteArray(ushort bitarray[], const int bitleng)
{
	assert(bitleng % 8 == 0);
	ushort* ps = bitarray;
	ushort* pe = (ushort*)((uchar*)ps + bitleng / 8 - 2);

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
static int shiftBitToLow(uchar bitarray[], const int bitleng, const int leftbitshift)
{
	uint* plowdword = (uint*)bitarray;
	uint* phigdword = (uint*)(bitarray + 2);

	const int rightbitshift = 16 - leftbitshift;

	//bugs for -O3 and utype = uchar
	for (int i = bitleng / 32 + 1; i > 0; i--) {
		*plowdword++ = ((ushort)(*plowdword >> leftbitshift)) |
			(*phigdword++ << rightbitshift);
	}

	return bitleng / 8;
}

//reverse bit array bitarray with length = bitleng, bugs
static void reverseBitArray(utype bitarray[], const int bitleng)
{
	const int bitmod = bitleng % 8;
	if (bitmod == 0) {
		reverseByteArray((ushort*)bitarray, bitleng);
	} else {
		reverseByteArray((ushort*)bitarray, bitleng + 8 - bitmod);
		shiftBitToLow((uchar*)bitarray, bitleng + bitmod, 8 - bitmod);
	}
	packQwordBit1(bitarray, bitleng);
}

//make sure no divide overflow
//improvement of 100%
static inline int
asmMulDiv(const uint startp, const uint pattern, uint p)
{
#ifdef LINE_COVER
	p = ((uint64)startp) * pattern % p;
#elif !defined _MSC_VER
	__asm
	(
#if 0
		"imul %%edx\n"
		"divl %%ecx\n"
		: "=d" (p)
		: "d"(startp), "a"(pattern), "c"(p)
#else
		"movl %1, %%eax\n"
		"imull %2\n"
		"divl %0\n"
		"movl %%edx, %0\n"
		: "+m" (p)
		: "g"(pattern), "g"(startp)
		: "%eax", "%edx"
#endif
	);
#else
	__asm
	{
		mov eax, pattern
		imul startp
		div p
		mov p, edx
	}
#endif

	return p;
}

static inline int
asmMulDivSub(const uint startp, const uint pattern, uint p, const uint bitleng)
{
#ifdef LINE_COVER
	p = (((uint64)startp) * pattern - bitleng) % p + bitleng;
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
		: "g"(pattern), "g"(startp), "g"(bitleng)
		: "%eax", "%ecx", "%edx"
#else
		"imull %3\n"

		"subl %%ecx, %%eax\n"
		"sbbl $0, %%edx\n"

		"idivl %4\n" //no overflow!!!
		"addl %%ecx, %%edx\n"
		: "=d" (p)
		: "a"(pattern), "c"(bitleng), "g"(startp), "g"(p)
#endif
	);
#else
	__asm
	{
		mov eax, pattern
		mov ecx, bitleng

		imul startp //startp * pattern
		sub eax, ecx // - bitleng
		sbb edx, 0

		idiv p //p
		add edx, ecx // + bitleng
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
			Prime[primes++] = p - lastprime;
			lastprime = p;

			if (p > sqrtn / p) {
				continue;
			}
			for (uint j = p * p / 2; j <= sqrtn / 2; j += p) {
				SET_BIT(bitarray, j);
			}
		}
	}

	//pack the last two number
	Prime[primes + 0] = Prime[primes + 1] = 255;

	if (0 && CHECK_FLAG(PRINT_LOG)) {
		printf("Prime[%d] = %d\n", primes, lastprime);
	}

	return primes;
}

//init bit tables
static void initBitTable()
{
	//1. init WordNumBit1 table in 0-2^16, can use popcnt replace it
	int i, nsize;

#if 0 == POPCNT && 0 == TREE2
	nsize = sizeof(WordNumBit1) / sizeof(WordNumBit1[0]);
	WordNumBit1[0] = 0;
	for (i = 1; i < nsize; i++) {
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
	for (i = 1; i < (1 << 16); i++) {
		WordReverse[i] = bytereverse[i >> 8] | (bytereverse[i & 255] << 8);
	}
#endif

#if 0 == FAST_CROSS
	//3. init CrossedTpl table, pre sieve the factor in array sievefactor
	sieveWheelFactor(CrossedTpl, 0, sizeof(CrossedTpl) * 16, SEGMENT_SIZE);
#endif
}

//
static int savePi1Pattern(const utype bitarray[], int leng, ptype pattern[])
{
	int pn = 0, lastpattern = 0;
	int diff = pattern[0];

	for (int p = 1; p < leng; p += 2) {
		if (!TST_BIT2(bitarray, p)) {
			pattern[pn++] = p - lastpattern;
			lastpattern = p;
		}
	}

	pattern[0] += diff;
	pattern[pn] = leng - lastpattern;

	return pn;
}

static int savePi2Pattern(const utype bitarray[], int leng, ptype pattern[])
{
	int pn = 0, lastpattern = 0;
	int diff = pattern ? pattern[0] : 0;

	for (int p = 1; p < leng; p += 2) {
		if (!TST_BIT2(bitarray, p) && !TST_BIT2(bitarray, p + Config.Kgap)) {
			if (pattern) {
				pattern[pn] = p - lastpattern;
				lastpattern = p;
			}
			pn++;
		}
	}

	if (pattern) {
		pattern[0] += diff;
		pattern[pn] = leng - lastpattern;
	}

	return pn;
}

static int savePikPattern(const utype bitarray[], int leng, ptype pattern[])
{
	int pn = 0, lastpattern = 0;
	int diff = pattern ? pattern[0] : 0;
	const int pis = KData.Kpattern[0];

	for (int p = 1; p < leng; p += 2) {
		if (TST_BIT2(bitarray, p))
			continue;

		int pasp = 0;
		for (int i = 1; i < pis; i++) {
			if (TST_BIT2(bitarray, p + KData.Kpattern[i])) {
				pasp = 1;
				break;
			}
		}

		if (pasp == 0) {
			if (pattern) {
				pattern[pn] = p - lastpattern;
				lastpattern = p;
				assert(p - lastpattern > 1 << (8 * sizeof(ptype)));
			}
			pn++;
		}
	}

	if (pattern) {
		pattern[0] += diff;
		pattern[pn] = leng - lastpattern;
	}

	return pn;
}

//
static int getPi2Pattern(const int factorial, ptype pi2pattern[])
{
	int sleng = SEGMENT_SIZE, patterns = 0;

	for (int start = 0; start < factorial; start += sleng) {
		utype bitarray[(SEGMENT_SIZE + 1000) >> (BSHIFT + 1)];
		if (start + sleng >= factorial)
			sleng = (int)(factorial - start);

		memset(bitarray, 0, (sleng >> 4) + 1 + (Config.Kgap >> 4));
		//memset(bitarray, 0, (sleng >> 4) + 1);
		sieveWheelFactor(bitarray, start, sleng + Config.Kgap, factorial);
		if (pi2pattern)
			patterns += savePi2Pattern(bitarray, sleng, pi2pattern + patterns);
		else
			patterns += savePi2Pattern(bitarray, sleng, 0);
	}

	if (pi2pattern) {
		if (CHECK_FLAG(PRINT_LOG)) {
			printf("factorial pattern pi2n = %d\n", patterns);
		}
	}

	return patterns;
}

//
static int getPi1Pattern(const int factorial, ptype pi1pattern[])
{
	int sleng = SEGMENT_SIZE, patterns = 0;
	for (int start = 0; start < factorial; start += sleng) {
		utype bitarray[(SEGMENT_SIZE + 1000) >> (BSHIFT + 1)];
		if (start + sleng >= factorial)
			sleng = (int)(factorial - start);

		memset(bitarray, 0, (sleng >> 4) + 1);
		sieveWheelFactor(bitarray, start, sleng, factorial);
		if (pi1pattern)
			patterns += savePi1Pattern(bitarray, sleng, pi1pattern + patterns);
		else
			patterns += countBitZeros(bitarray, sleng / 2);
	}

	if (pi1pattern) {
		if (CHECK_FLAG(PRINT_LOG)) {
			printf("factorial pattern pi1n = %d\n", patterns);
		}
		assert(1 == pi1pattern[0] && pi1pattern[patterns] == 1);
	}

	return patterns;
}

//
static int getPikPattern(const int factorial, ptype pikpattern[])
{
	int sleng = SEGMENT_SIZE, patterns = 0;

	for (int start = 0; start < factorial; start += sleng) {
		utype bitarray[(SEGMENT_SIZE + 1000) >> (BSHIFT + 1)];
		if (start + sleng >= factorial)
			sleng = (int)(factorial - start);

		memset(bitarray, 0, (sleng >> 4) + 1 + (Config.Kgap >> 4));
		sieveWheelFactor(bitarray, start, sleng + Config.Kgap, factorial);
		if (pikpattern) {
			patterns += savePikPattern(bitarray, sleng, pikpattern + patterns);
		} else {
			patterns += savePikPattern(bitarray, sleng, 0);
		}
	}

	if (pikpattern) {
		if (CHECK_FLAG(PRINT_LOG)) {
			printf("factorial pattern pi%d = %d\n", KData.Kpattern[0], patterns);
		}
	}

	return patterns;
}

static int initFirstPos(int s[], int startp, int p, int bitleng)
{
	int ns = KData.Kpattern[0];
	int si = s[1];
	int mod2p = startp * 2;
	if (mod2p > p) {
		mod2p -= p;
	}
	int mod4p = mod2p * 2;
	if (mod4p > p) {
		mod4p -= p;
	}

	for (int i = 1; i < ns; i++) {
		int pdiff = KData.Kpattern[i] - KData.Kpattern[i - 1];
		if (pdiff == 2) {
			si += mod2p;
		} else if (pdiff == 4) {
			si += mod4p;
		} else if (pdiff == 6) {
			si += mod4p + mod2p;
			if (si > bitleng)
				si -= p;
		} else {
			si = s[1] + startp * KData.Kpattern[i] % p;
		}
		if (si > bitleng)
			si -= p;
		s[i + 1] = si;
	}

	return ns;
}

static int sievePi1L1(utype bitarray[], const uint pattern, int bitleng, int& p)
{
	const int minp = KData.SqrtN < Config.CpuL1Size ? KData.SqrtN : Config.CpuL1Size;
	int k = KData.firstIndex + 1;
	int spos[MAX_L1SIZE / 11];

	for (; p <= minp; PRIME_NEXT(p, k)) {
		spos[k] = asmMulDiv(Moudle[k], pattern, p);
	}
	//		assert(k < sizeof(spos));
	for (int start = 0, sleng = Config.CpuL1Size; start < bitleng; start += Config.CpuL1Size) {

		k = KData.firstIndex;
		p = SmallPrime[k++];

		if (start + Config.CpuL1Size > bitleng) {
			sleng = bitleng - start;
		}

		for (; p <= minp; PRIME_NEXT(p, k)) {
			int npos = spos[k];
			if (npos < sleng) {
				npos = setBitArray(bitarray + (start >> BSHIFT), npos, p, sleng);
			}
			spos[k] = npos - sleng;
		}
	}

	return k;
}

#if 1
static int sievePi2L1(utype bitarray[], const uint pattern1, const uint pattern2, int bitleng, int& p)
{
	const int minp = KData.SqrtN < Config.CpuL1Size ? KData.SqrtN : Config.CpuL1Size;
	int k = KData.firstIndex + 1;
	uint64 spos[MAX_L1SIZE / 11];

	for (; p <= minp; PRIME_NEXT(p, k)) {
		const uint s1 = asmMulDiv(Moudle[k], pattern1, p);
		const uint s2 = asmMulDiv(Moudle[k], pattern2, p);
		spos[k] = (uint64)s2 << 32 | s1;
	}

	for (int start = 0, sleng = Config.CpuL1Size; start < bitleng; start += Config.CpuL1Size) {

		k = KData.firstIndex;
		p = SmallPrime[k++];

		if (start + Config.CpuL1Size > bitleng) {
			sleng = bitleng - start;
		}

		for (; p <= minp; PRIME_NEXT(p, k)) {
			spos[k] = set2BitArray(bitarray + (start >> BSHIFT), spos[k], p, sleng);
		}
	}

	return k;
}
#else
static int sievePi2L1(utype bitarray[], const uint pattern1, const uint pattern2, int bitleng, int& p)
{
	const int minp = KData.SqrtN < Config.CpuL1Size ? KData.SqrtN : Config.CpuL1Size;
	int k = KData.firstIndex + 1;
	int sleng = Config.CpuL1Size;
	uint64 spos[MAX_L1SIZE / 11];

	for (; p <= minp; PRIME_NEXT(p, k)) {
		int s1 = asmMulDivSub(Moudle[k], pattern1, p, bitleng);
		if (s1 > bitleng) {
			s1 -= p;
		}

		int s2 = asmMulDivSub(Moudle[k], pattern2, p, bitleng);
		if (s2 > bitleng) {
			s2 -= p;
		}

		s1 -= (bitleng - Config.CpuL1Size);
		s2 -= (bitleng - Config.CpuL1Size);

		spos[k] = (uint64)s2 << 32 | s1;
	}

	for (int start = bitleng; start > 0; start -= sleng) {

		k = KData.firstIndex;
		p = SmallPrime[k++];

		if (start < sleng) {
			start = sleng;
		}
		for (; p <= minp; PRIME_NEXT(p, k)) {
			spos[k] = setBitArray(bitarray + ((start - sleng) >> BSHIFT), spos[k], -p, sleng);
		}
	}

	return k;
}
#endif

static int sievePikL1(utype bitarray[], const uint pattern, int bitleng, int& p)
{
	int k = KData.firstIndex + 1;
#if 1
	int ns = KData.Kpattern[0];
	if (ns & 1) {
		k = sievePi1L1(bitarray, pattern, bitleng, p);
		ns -= 1;
	}
	for (int j = ns; j > 0; j -= 2) {
		k = KData.firstIndex;
		p = SmallPrime[k++];
		k = sievePi2L1(bitarray, pattern + KData.Kpattern[j], pattern + KData.Kpattern[j - 1], bitleng, p);
	}
#else

	const int minp = KData.SqrtN < Config.CpuL1Size ? KData.SqrtN : Config.CpuL1Size;
	const int ns = KData.Kpattern[0];
	int spos[MAX_L1SIZE / 11][8];

	for (; p <= minp; PRIME_NEXT(p, k)) {
		spos[k][1] = asmMulDiv(Moudle[k], pattern, p);
		for (int j = 1; j < ns; j ++) {
			spos[k][j + 1] = asmMulDiv(Moudle[k], KData.Kpattern[j] + pattern, p);
		}
	}

	for (int start = 0, sleng = Config.CpuL1Size; start < bitleng; start += sleng) {

		if (start + sleng > bitleng) {
			sleng = bitleng - start;
		}

		for (int j = 1; j <= ns; j++) {
			k = KData.firstIndex;
			p = SmallPrime[k++];
			for (; p <= minp; PRIME_NEXT(p, k)) {
				int npos = spos[k][j];
				if (npos < sleng) {
					npos = setBitArray(bitarray + (start >> BSHIFT), npos, p, sleng);
				}

				spos[k][j] = npos - sleng;
			}
		}
	}
#endif

	return k;
}

//30%
static int sievePi1(utype bitarray[], const int pattern)
{
	const int bitleng = 1 + (int)((KData.N - pattern) / KData.Wheel);
	const int sqrtn = KData.SqrtN;

	int k = KData.firstIndex;
	int p = SmallPrime[k++];

#if OPT_L1CACHE
	// performance improvement from 100 -> 61
	if (bitleng > Config.CpuL1Size)
		k = sievePi1L1(bitarray, pattern, bitleng, p);
#endif

	for (; p <= sqrtn; PRIME_NEXT(p, k)) {
#if 0
		int s1 = asmMulDiv(Moudle[k], pattern, p);
		if (s1 < bitleng)
			setBitArray(bitarray, s1, p, bitleng);
#else
		int s1 = asmMulDivSub(Moudle[k], pattern, p, bitleng);
		if (s1 > bitleng) {
			s1 -= p; //			if (s1 >= bitleng) s1 -= p;
		}

		if (s1 > 0)
			setBitArray(bitarray, s1, -p);
#endif
	}

	return bitleng;
}

//core part of the algorithm and code will be ported to CUDA
//set bitarray flag with k * wheel + pattern
static int sievePi2(utype bitarray[], const uint pattern)
{
	const int kgap = Config.Kgap;
	const int bitleng = 1 + ((KData.N - kgap - pattern) / KData.Wheel);
	const int sqrtn = KData.SqrtN;

	int k = KData.firstIndex;
	int p = SmallPrime[k++];

#if OPT_L1CACHE
	//performance improvement from 10 -> 72
	if (bitleng > Config.CpuL1Size) {
		k = sievePi2L1(bitarray, pattern, pattern + kgap, bitleng, p);
	}
#endif

	for (; p <= sqrtn; PRIME_NEXT(p, k)) {
		int moudle = Moudle[k];
		int s1 = asmMulDivSub(moudle, pattern, p, bitleng);
		if (s1 > bitleng) {
			s1 -= p;
		}
		int s2 = s1 + moudle * 2;
		if (kgap > 2) {
			s2 = s1 + moudle * kgap % p;
		} else if (s2 >= bitleng) {
			s2 -= p;
		}
		if (s2 >= bitleng)
			s2 -= p;
		set2BitArray(bitarray, s1, s2, -p);
	}

	return bitleng;
}

static int sievePik(utype bitarray[], const uint pattern)
{
	const int bitleng = 1 + (KData.N - Config.Kgap - pattern) / KData.Wheel;

	const int sqrtn = KData.SqrtN;
	int k = KData.firstIndex;
	int p = SmallPrime[k++];
	int s[16];

#if OPT_L1CACHE
	// performance improvement from 100 -> 61
	if (bitleng > Config.CpuL1Size) {
		k = sievePikL1(bitarray, pattern, bitleng, p);
	}
#endif

	for (; p <= sqrtn; PRIME_NEXT(p, k)) {
		int s1 = asmMulDivSub(Moudle[k], pattern, p, bitleng);
		if (s1 > bitleng) {
			s1 -= p;
		}

		s[1] = s1;
		int ns = initFirstPos(s, Moudle[k], p, bitleng);
		if (ns & 1)
			setBitArray(bitarray, s[ns--], -p);
		for (int j = ns; j > 0; j -= 2) {
			set2BitArray(bitarray, s[j], s[j - 1], -p);
		}
	}

	return bitleng;
}

//bad performance !!!!
static int countKtable(const ushort bitarray[], const int bitleng, const int pattern, uint64 tdata[])
{
	const uint64 base = atoint64(TABLE_GAP, 0);
	const int wordleng = bitleng / 16;
	const int diff = pattern + KData.Kpattern[KData.Kpattern[0]];
	int ktupels = 0;

	for (int b = 0; b <= wordleng; b ++) {
		ushort masks = (ushort)(~bitarray[b]);
		while (masks != 0) {
			const int bitindex = LeftMostBit1[masks];
			const uint64 P = (uint64)KData.Wheel * (b * 16 + bitindex) + diff;
			ktupels++;
			tdata[(int)(P / base)]++;
			masks &= masks - 1;
		}
	}

	return ktupels;
}

static ptype* getNextPattern(int& pattern, ptype* pnext)
{
	ptype pdiff = *pnext++;
	pattern += pdiff;
	if (pdiff == 0) {
		pattern += KData.PatternDiff;
		pnext = Pattern + 1;
	}

	return pnext;
}

static void printProgress(const int tid, const double tstart, int aves, int pcnt)
{
	double currper = 100.0 * pcnt / KData.Patterns;
	double totaltime = (getTime() - tstart) / (10 * currper);
	static uint64 lastValue = 0;
	if (pcnt <= Config.PrintGap) {
		lastValue = 0;
	}

	printf("thread(%d) %.2lf%%, ~= ", tid, currper);
	if (totaltime < 10000) {
		printf("%.2lf sec", totaltime);
	} else {
		printf("%.2lf hour", totaltime / 3600);
	}

	uint64 curValue = ((uint64)aves) * KData.Patterns;
	printf(", %s ~= %lld", KtupletName[Config.Ptk], curValue);
	if (lastValue > 0) {
		printf(", err ~= %.4lf%%%%",
				(curValue - lastValue) * 10000.0 / lastValue);
	}
	putchar('\n');

	lastValue = curValue;
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

	uint64 gpts = 0;
	double tstart = getTime();

	const uint sieve_byte = ((KData.N / KData.Wheel) / 8) + 1024 / 4;
	utype* bitarray = (utype*)malloc(sieve_byte);
	assert(bitarray && Pattern);

	uint64 *tdata = 0;
	if (KData.UseKtable) {
		tdata = (uint64*)malloc(100000 * sizeof(uint64));
		memset(tdata, 0, sizeof(tdata[0]) * 100000);
	}

	ptype* pnext = Pattern;
	int pattern = 0;
	for (int i = 0; i < pbegi; i++) {
		pnext = getNextPattern(pattern, pnext);
	}

	if (CHECK_FLAG(PRINT_LOG)) {
		printf("thread %d : pattern %d - %d\n", tid, pbegi, pendi);
	}

	for (int pcuri = pbegi; pcuri < pendi; pcuri++) {

		pnext = getNextPattern(pattern, pnext);
#if 1
		if (CHECK_FLAG(CHCECK_PATTERN)) {
			if (gcd(KData.Wheel, pattern) != 1)
				printf("error pattern = %d\n", pattern);
			continue;
		}
#endif
		memset(bitarray, 0, sieve_byte);
		int bitleng = 0;
		if (PI2N == Config.Ptk) {
			bitleng = sievePi2(bitarray, pattern);
		} else if (PI1N == Config.Ptk) {
			bitleng = sievePi1(bitarray, pattern);
		} else {
			bitleng = sievePik(bitarray, pattern);
		}

		bitarray[0] |= 1;
		gpts += countBitZeros(bitarray, bitleng);

		if (tdata) {
			countKtable((ushort*)bitarray, bitleng, pattern, tdata);
		}

#ifdef _M_AMD64
		InterlockedAdd((LONG*)(&scnt), 1);
#else
		scnt++;
#endif
		if ((scnt & Config.PrintGap) == Config.PrintGap - 1) {
			printProgress(tid, tstart, (int)(gpts / (1 + pcuri - pbegi)), scnt);
		}
	}

	free(bitarray);

	if (tdata) {
		for (int i = 0; tdata[i]; i++) {
			Ktable[i] += tdata[i];
		}
		free(tdata);
	}

	if (CHECK_FLAG(PRINT_LOG)) {
		printf("Thread %d: pattern[%3d - %3d] = %lld\n", tid, pbegi, pendi, gpts);
	}

	return gpts;
}

//factorial of prime factor of wheel
static int getFactorial(const int wheel)
{
	int factorial = 1;
	for (int i = 0, p = SmallPrime[i]; wheel % p == 0; p = SmallPrime[++i]) {
		factorial *= p;
	}

	return factorial;
}

//
static int countKpattern(int wheel, const int kgap)
{
	int patterns = 3;
	if (kgap == 1) {
		patterns = 8;
	}
	for (int i = 3, p = SmallPrime[i]; wheel % p == 0; p = SmallPrime[++i]) {
		wheel /= p;
		patterns *= p - kgap;
	}

	return patterns * (wheel / 30);
}

//get the frist prime index in Prime which is not a factor of wheel
static int getFirstPrime(const int wheel)
{
	int i = 0;
	for (int p = SmallPrime[i]; wheel % p == 0; p = SmallPrime[++i]) {

	}
	return i;
}

//16 bits number p1, p2 and gcd((p1 + p2), factorial) = 1
// ---------- gp ------------------*/
static int initPattern(const int factorial, const int wheel)
{
	double ts = getTime();
	int pns = countKpattern(factorial, 1);
	Pattern = (ptype*)malloc(10000 + (sizeof(ptype) * (pns + 1) >> Config.Ptk));
	Pattern[0] = 0;

	if (PI1N == Config.Ptk) {
		pns = getPi1Pattern(factorial, Pattern);
	} else if (PI2N == Config.Ptk) {
		pns = getPi2Pattern(factorial, Pattern);
		assert(pns == countKpattern(factorial, 2) && Config.Kgap < 6);
	} else {
		pns = getPikPattern(factorial, Pattern);
	}

	KData.PatternDiff = Pattern[pns] + Pattern[0];

	if (CHECK_FLAG(PRINT_LOG)) {
		const int ns = KData.Kpattern[0];
		printf("init pattern time use %.2lf ms, count %.3lf, %.2lf fast than pi\n",
				getTime() - ts, 100.0 * pns * ns / factorial,
				1.0 * countKpattern(factorial, 1) / (pns * ns));
	}
	if (CHECK_FLAG(CHCECK_PATTERN)) {
		printf("wheel pattern = %d\n", pns);
	}

	Pattern[pns] = 0;
	pns *= wheel / factorial;
	return pns;
}

static void initMoudle(const int wheel, int maxp, uint Moudle[])
{
	double ts = getTime();

	for (int j = 1, p = 2; p < maxp; PRIME_NEXT(p, j)) {
		int y = 0;
		y = extendedEuclid(-wheel % p, p, y);
		if (y < 0) {
			y += p;
		}
		Moudle[j] = y;
	}

	if (CHECK_FLAG(PRINT_LOG)) {
		printf("init startp time use %.2lf ms\n", getTime() - ts);
	}
}

static uint getDefaultWheel(const uint64 n)
{
	uint wheel = KData.Wheel;
	if (wheel > 30 && wheel % 30 == 0 && n > wheel) {
		return wheel;
	}

	wheel = 1;
	for (int i = 0; SmallPrime[i] <= KData.Wheel && SmallPrime[i] < 24;) {
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
		wheel = 9699690 * 23;
	} else {
		wheel = 223092870ul * 29;
	}

	return wheel;
}

//set sieve buffer size and adjust wheel based
//on cpu L2 cache size and n
/*
static int getSieveCacheSize(const uint64 n, int wheel)
{
	const int cachesize = n / wheel;
	return (cachesize >> BSHIFT) + 1;
}*/

//set sieve buffer size and adjust wheel based
//on cpu L2 cache size and n
static int getWheel(const uint64 n)
{
	int wheel = getDefaultWheel(n);

	int cachesize = n / wheel;

	int blocks = cachesize / (Config.CpuL2Size << 13);

	wheel *= (blocks + 1);

	return wheel;
}

//get prime number with diff Result in array Prime
//segmented sieve of EratoSieve to enum prime number
static int getPrime(const int s, const int n, uchar prime[])
{
	int pi1n = 1;
	if (prime) {
		prime[2] = 1 - 3;
	}

	for (int start = s, sleng = SEGMENT_SIZE; start < n; start += sleng) {
		utype bitarray[(SEGMENT_SIZE >> (BSHIFT + 1)) + 100];
		if (sleng >= n - start) {
			sleng = n - start;
		}
#if FAST_CROSS
		segmentedSieve1(bitarray, start, sleng + 16);
#else
		segmentedSieve2(bitarray, start, sleng + 16);
#endif
		if (prime) {
			pi1n += savePi1Pattern(bitarray, sleng, prime + pi1n + 1);
		} else {
			pi1n += countBitZeros(bitarray, sleng / 2);
		}
	}

	assert(pi1n < sizeof(Prime) / sizeof(Prime[0]));
	if (CHECK_FLAG(PRINT_LOG)) {
		printf("pi(%d) = %d\n", n, pi1n);
	}

	return pi1n;
}

static int getSegKtuplet(const int s, const int n)
{
	int ret = 0;

#if (OMP)
	omp_set_num_threads(Config.Threads);
	#pragma omp parallel for reduction(+:ret) if (n >= MINN * 3)
#endif
	for (int start = s; start < n; start += SEGMENT_SIZE) {
		utype bitarray[(SEGMENT_SIZE >> (BSHIFT + 1)) + 100];
		int sleng = SEGMENT_SIZE;
		if (sleng >= n - start) {
			sleng = n - start;
		}
#if FAST_CROSS
		segmentedSieve1(bitarray, start, sleng + 16);
#else
		segmentedSieve2(bitarray, start, sleng + 16);
#endif
		if (PI2N == Config.Ptk)
			ret += savePi2Pattern(bitarray, sleng, 0);
		else
			ret += savePikPattern(bitarray, sleng, 0);
	}

	return ret;
}

//get small parathion or ktuplet prime in range[0 - min(Wheel, sqrt(n))]
//if n is less than a small fix value MIN
static int getSmallGpt(const uint64 s, const uint64 n)
{
	double ts = getTime();

	//adjust leng for last few number
	int leng = 0;
	if (n < MINN) {
		if (PI1N == Config.Ptk) {
			leng = n;
		} else if (PI2N >= Config.Ptk) {
			leng = n - (Config.Kgap - 1);
		}
	} else {
		leng = KData.SqrtN;
		if (leng < KData.Wheel) {
			leng = KData.Wheel;
		} else if (PI2N <= Config.Ptk) {
			leng += 1;
		}
	}

	int ret = 0;
	if (PI1N == Config.Ptk) {
		ret = getPrime(s ,leng + 1, 0);
	} else if (PI2N <= Config.Ptk) {
		ret = getSegKtuplet(s, leng);
	}

	if (CHECK_FLAG(PRINT_LOG)) {
		printf("sieve small n = %lld, leng = %d", n, leng);
		printf("\nand small ret = %d, and time use %.lf ms\n", ret, getTime() - ts);
	}

	return ret;
}

//init Prime, Pattern, Moudle
static int initParam(const uint64 s, const uint64 n)
{
	const int wheel = getWheel(n - s);
	KData.N = n;
	KData.S = s;
	KData.Wheel = wheel;
	KData.SqrtN = (int)sqrt((double)n + 0.1);

	const int factorial = getFactorial(wheel);
	KData.Factorial = factorial;
	KData.firstIndex = getFirstPrime(wheel);
	KData.Patterns = initPattern(factorial, wheel);

	initMoudle(wheel, KData.SqrtN + 10, Moudle);

	LastTask.Ptk = Config.Ptk;
	LastTask.Wheel = KData.Wheel;
	LastTask.N = KData.N;
	LastTask.S = KData.S;

	if (CHECK_FLAG(PRINT_LOG)) {
		printf("wheel = %d * %d, ", factorial, wheel /factorial);
		printf("cachesize = %d k\n", (n / wheel) >> 13);
	}

	return wheel;
}

static void saveTask(struct Task& curtask)
{
	if (freopen("prime.ta", "rb", stdin)) {
		freopen(CONSOLE, "r", stdin);
		remove("prime.ta.bak");
		if (rename("prime.ta", "prime.ta.bak") != 0) {
			perror("bake data fail");
		}
	}

	freopen("prime.ta", "wb", stdout);

	if (LastTask.Tasks == 0)
		curtask.Pendi += KData.Patterns / 3;
	else
		curtask.Pendi += KData.Patterns / LastTask.Tasks;

	if (curtask.Pendi > KData.Patterns) {
		curtask.Pendi = KData.Patterns;
	}

	printf(TaskFormat, curtask.Ptk,
			Config.Kgap, curtask.Wheel, KData.Patterns,
			curtask.Tasks, curtask.Pbegi, curtask.Pendi,
			curtask.N, curtask.Result);

	freopen(CONSOLE, "r", stdin);
	freopen(CONSOLE, "w", stdout);
}

static int readTaskData(struct Task &curtask)
{
	int ret = 0;
	char linebuf[400] = {0};
	char taskdata[400] = {0};
	freopen("prime.ta", "rb", stdin);

	while (gets(linebuf)) {
		strcat(taskdata, strcat(linebuf, "\n"));
	}

	//read last task data
	struct Task tmp = {0};
	if (sscanf(taskdata, TaskFormat, &tmp.Ptk,
		&Config.Kgap, &tmp.Wheel, &KData.Patterns,
		&curtask.Tasks,	&curtask.Pbegi, &curtask.Pendi,
		&tmp.N, &curtask.Result) != 9) {
		printf("invalid task data format: %s : %s\n", taskdata, TaskFormat);
		ret = -1;
	}

	//check task data
	if (curtask.N != tmp.N ||
		curtask.Wheel != tmp.Wheel ||
		curtask.Ptk != tmp.Ptk) {
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

	if (readTaskData(curtask) != 0) {
		curtask.Result = curtask.Pbegi = curtask.Pendi = 0;
		saveTask(curtask);
	}

	if (curtask.Tasks > 0) {
		curtask.Pendi = curtask.Pbegi + KData.Patterns / curtask.Tasks + 1;
	}

	if (curtask.Pendi > KData.Patterns) {
		curtask.Pendi = KData.Patterns;
	}

	if (CHECK_FLAG(PRINT_LOG)) {
		printf("load last Task Data with pattern[%d - %d] ok\n", curtask.Pbegi, curtask.Pendi);
	}

	freopen(CONSOLE, "r", stdin);

	return ret;
}

//
static uint64 getGpKtPn(const uint64 s, const uint64 n, int pn, bool addsmall)
{
	assert(n < atoint64(MAXN, 100) + 1000 && s < n);
	assert(Config.Ptk >= 0 && Config.Ptk < 10);

	if (n >= MINN) {
		initParam(s, n);
	}

	uint64 sgn = 0, gptn = 0;
	if (addsmall) {
		sgn = getSmallGpt(s, n);
		if (KData.UseKtable) {
			Ktable[0] += sgn;
		}
	}

	if (n >= MINN) {
		if (pn <= 0 || pn > KData.Patterns) {
			pn = KData.Patterns;
		}
		int pbegi = 0;
		int pendi = pn;

		//load Last Task
		if (CHECK_FLAG(SAVE_TASK) && loadTask(LastTask) >= 0) {
			pbegi = LastTask.Pbegi;
			pendi = LastTask.Pendi;
			gptn = LastTask.Result;
		}

#if (OMP)
		omp_set_num_threads(Config.Threads);
		#pragma omp parallel for reduction(+:gptn) if (n >= MINN * 3)
		for (int oi = 0; oi < Config.Threads; oi++) {
			int bi = pn / Config.Threads * oi;
			int ei = bi + pn / Config.Threads;
			if (oi == Config.Threads - 1) {
				ei = pn;
			}
			gptn += sievePattern(bi + 1, ei);
		}
#else
		if (pendi - pbegi > 6 && Config.Threads > 1) {
			gptn += startWorkThread(Config.Threads, pbegi, pendi);
		} else if (pendi > pbegi) {
			gptn += sievePattern(pbegi, pendi);
		}
#endif

		//save Current Task
		if (CHECK_FLAG(SAVE_TASK) && pbegi < pendi) {
			LastTask.Pbegi = pendi;
			LastTask.Result = gptn;
			saveTask(LastTask);
		}
	}

	if (Pattern) {
		free(Pattern);
		Pattern = NULL;
	}

	KData.Wheel = 0;

	return gptn + sgn;
}

static void printResult(const uint64 n, uint64 gptn, double ts)
{
	if (KData.UseKtable) {
		putchar('\n');
		uint64 sum = 0;
		for (int i = 0; Ktable[i]; i++) {
			sum += Ktable[i];
			printf("%4d%s %lld\n", i + 1, TABLE_GAP, sum);
		}
	}

	if (Config.Ptk < PIKN) {
		printf(PrintFormat[Config.Ptk], n, gptn);
	} else {
		printf(PrintFormat[PIKN], KData.Kpattern[0], n, gptn);
	}

	if (CHECK_FLAG(PRINT_TIME)) {
		printf(" (%.2lf sec)", (getTime() - ts) / 1000);
	}

	//		if (CHECK_FLAG(PRINT_LOG)) {
	//			printf(", ave = %d", (int)(gptn / KData.Patterns));
	//		}
	putchar('\n');
}

//get twin primes, or ktuplet
static uint64 ktprime(const uint64 s, const uint64 n, int pn)
{
	double ts = getTime();

	uint64 gptn = 0;
	if (KData.Patterns == 0) {
		LastTask.Ptk = -1;
	}

	const int sqrtn = (int)sqrt((double)n + 0.1);
	if (n > KData.N + 1000) {
		getPrime(0, sqrtn + 1001, Prime);
	}

	int wheel = getWheel(n);
	const int maxn = sqrtn > wheel ? sqrtn : wheel;
	bool addsmall = s < maxn;
	//optimize for Ktuplet prime with small n > e7
	if (addsmall && n > atoint64("e13", 0)) {
		wheel = KData.Wheel;
		KData.Wheel = 0;
		gptn = getGpKtPn(s, maxn + Config.Kgap - 1, 0, addsmall);
		KData.Wheel = wheel;
		if (CHECK_FLAG(PRINT_LOG)) {
			//			printf("PI%d[%d] = %lld, and time use %.lf ms\n",
			//					Config.Ptk, maxn, gptn, getTime() - ts);
		}
		addsmall = false;
	}

	if (Config.Ptk == PIKN && pn == 0 && (n - s) > atoint64(TABLE_GAP, 0)) {
		memset(Ktable, 0, sizeof(Ktable));
		KData.UseKtable = true;
	} else {
		KData.UseKtable = false;
	}

	gptn += getGpKtPn(s, n, pn, addsmall);

	if (CHECK_FLAG(PRINT_RET)) {
		printResult(n, gptn, ts);
	}

	return gptn;
}

//test case number
static int startTest(int tesecase, bool rwflag)
{
	srand((uint)time(NULL));
	double ts = getTime();
	const char* dataFile[ ] = {"prime.gp", "prime.pi", "prime.pi2", "prime.pik"};
	printf("-----------start %s test -----------\n", dataFile[Config.Ptk]);

	if (rwflag) {
		if (!(freopen(dataFile[Config.Ptk], "rb", stdin))) {
			printf("can not read test data file %s\n", dataFile[Config.Ptk]);
			freopen(CONSOLE, "r", stdin);
			SET_FLAG(PRINT_RET);
			return -1;
		}
		CLR_FLAG(PRINT_RET | PRINT_LOG);
		Config.PrintGap = 0;
	} else {
		if (!(freopen(dataFile[Config.Ptk], "wb", stdout))) {
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
			uint64 n = rand() * rand();
			n = (n + 4) * 2 + 0;
			if (n < 1000000) {
				n = 4 * n + 1000000;
			}
			ktprime(0, n, 0);
		} else {
			uint64 res, n;
			char linebuf[256] = {0};
			gets(linebuf);
			if (sscanf(linebuf, PrintFormat[Config.Ptk], &n, &res) != 2) {
				printf("line %d is wrong data\n", i);
				if (fails++ > 30) {
					break;
				} else {
					continue;
				}
			}

			uint64 gptn = ktprime(0, n, 0);
			if (gptn != res) {
				printf("case %d with wrong result %lld, ", i, gptn);
				printf(PrintFormat[Config.Ptk], n, res);
				putchar('\n');
			}
		}
		if ((i & 63) == 0) {
			printf("case pass %d%%\r", i * 100 / tesecase);
		}
	}

	printf("test case time use %.lf ms\n", getTime() - ts);

	SET_FLAG(PRINT_RET);
	freopen(CONSOLE, "w", stdout);
	freopen(CONSOLE, "r", stdin);

	return 0;
}

//list Ptk by the input Result start, end, step
static void listDiffGpt(const char cmdparams[][80], int cmdi)
{
	double ts = getTime();

	int ni = 1, step = 2;
	uint64 start = ipow(10, 9), end = start + 1000;
	uint64 buf[ ] = {0, start, end, step, 0};

	for (int i = cmdi; cmdparams[i][0] && ni < sizeof(buf) / sizeof(buf[0]); i++) {
		char c = cmdparams[i][0];
		if (isdigit(c) || toupper(c) == 'E') {
			buf[ni++] = atoint64(cmdparams[i], 10000);
		}
	}

	start = buf[1], end = buf[2], step = buf[3];

	start += (start & 1);
	step += step & 1;

	if (step < 2) {
		step = 2;
	}
	if (start > end) {
		end = start + end * step - 1;
	}

	printf("calculating %s\n", KtupletName[Config.Ptk]);

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
		if (isdigit(cmdparams[4][0])) {
			printf("%d ", pcnt);
		}
		allSum += ktprime(0, n, 0);
	}

	printf("all case time use %.lf ms\n", getTime() - ts);
	freopen(CONSOLE, "w", stdout);
}

static void listPowGpt(const char cmdparams[][80], int cmdi)
{
	printf("calculating %s ", KtupletName[Config.Ptk]);

	int m = atoint64(cmdparams[cmdi + 1], 10);
	int startindex = atoint64(cmdparams[cmdi + 2], 5);
	int endindex = atoint64(cmdparams[cmdi + 3], 10);

	printf("in %d^%d - %d^%d\n", m, startindex, m, endindex);

	if (m < 2 && m > 10000) {
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
		uint64 n = (uint64)(pow((double)m, (int)i) + 0.01);
		uint64 r = ktprime(0, n, 0);
		printf(PrintFormat[Config.Ptk], m, i, r);
		putchar('\n');
	}
	SET_FLAG(PRINT_RET);
}

static void listPatterns(uint64 start, int count)
{
	getPrime(0, (int)sqrt((double)start) + count * 4 + 2000, Prime);
	int wheel = getWheel(start);
	printf("wheel = %d\n", wheel);
	for (int i = 0; i < count; i++) {
		int patterns = 0;
		printf("pattern(%lld) = %d\n", start, patterns);
		start += 2;
	}
}

//benchMark
static void benchMark(const char cmdparams[][80])
{
	SET_FLAG(PRINT_RET);
	uint64 start = atoint64("e10", 1000000);
	uint64 n = atoint64("e10", 0);
	char gptk[20] = {'0', '1', '2'};

	for (int i = 0, ni = 1; cmdparams[i][0]; i++) {
		char c = cmdparams[i][0];
		if (isdigit(c) || toupper(c) == 'E') {
			uint64 tmp = atoint64(cmdparams[i], 0);
			if (ni++ == 1) {
				start = tmp;
			} else if (ni == 2) {
				n = tmp;
			} else if (tmp < 222) {
				strcpy(gptk, cmdparams[i]);
			}
		}
	}

	if (CHECK_FLAG(SAVE_RESUTL)) {
		freopen("benchmark.txt", "wb", stdout);
	}

	uint64 gap = start;
	for (; start * 10 <= n; ) {
		for (int j = 0; j < 9; j++) {
			for (int k = 0; gptk[k]; k++) {
				Config.Ptk = gptk[k] - '0';
				ktprime(0, (j + 1) * gap, 0);
			}
			puts("");
		}
		start *= 10;
		gap *= 10;
	}

	freopen(CONSOLE, "w", stdout);
}

//test pi, pi2, gp function
static void testPik()
{
	const char* const pidata[][5] =
	{
		{"k31", "1e9 379508", "e10 2713347", "e11 20093124", "e12 152850135"},
		{"k32", "1e9 379748", "e10 2712226", "e11 20081601", "e12 152839134"},
		{"k41", "1e9 028388", "e10 0180529", "e11 01209318", "e12 008398278"},
		{"k51", "1e9 003633", "e10 0020203", "e11 00122457", "e12 000776237"},
		{"k52", "1e9 003588", "e10 0020211", "e11 00122855", "e12 000775986"},
		{"k61", "e10 001613", "e11 0008626", "e12 00050408", "e13 000303828"},
		{"k71", "e10 000234", "e11 0001183", "e12 00006056", "e13 000033395"},
		{"k72", "e10 000239", "e11 0001152", "e12 00005913", "e13 000033066"}
	};

	Config.CpuL2Size = 1 << 10;
	Config.PrintGap = 0;
	SET_FLAG(PRINT_TIME);

	for (int i = 0; i < sizeof(pidata) / sizeof(pidata[0]); i++) {
		excuteCmd(pidata[i][0]);

		for (int j = 1; j < sizeof(pidata[i])/sizeof(pidata[i][0]); j ++) {
			uint64 n = atoint64(pidata[i][j], 100);
			uint64 r = atoint64(pidata[i][j] + 3, 0);
			if (r != ktprime(0, n, 0)) {
				printf("%s != %lld fail\n", pidata[i][j], r);
			}
		}
		putchar('\n');
	}
}

//test pi2 function
static void testPi2()
{
	//twin Cousin Sexty prime data from
	//Formeln zur Berechnung der Anzahl
	const char* const pi2data[] =
	{
		"10^3  35 41 74 38 51 69",
		"10^4  205 203 411 208 270 404",
		"10^5  1224 1216 2447 1260 1624 2420",
		"10^6  8169 8144 16386 8242 10934 16378",
		"10^7  58980 58622 117207 58595 78211 117486",
		"10^8  440312 440258 879908 439908 586811 880196",
		"01e9  3424506 3424680 6849047 3426124 4567691 6847940",
		"02e9  6388041 6386968 12773727 6387840 8520754 12775879",
		"04e9  11944438 11946709 23886965 11946187 15926546 23892168",
		"08e9  22384176 22383210 44759392 22384852 29844274 44769861",
		"1e10  27412679 27409999 54818296 27411508 36548839 54822710",
		"2e10  51509099 51501940 102999661 51507299 68669678 103008732",
		"4e10  96956707 96957512 193905500 96960771 129275019 193928169",
		"8e10  182855913 182860246 365695515 182857130 243792036 365720264",

		"1e11  224376048 224373161 448725003 224365334 299140330 448749360",
		"2e11  424084653 424077104 848122150 424069491 565421597 848162636",
		"4e11  802817718 802794952 1605598668 802794667 1070396987 1605600054",
		"8e11  1521998439 1521984605 3044016884 1521998098 2029285571 3043915677",
		"1e12  1870585220 1870585459 3741217498 1870580394 2494056601 3741051790"
	};

	Config.PrintGap = 0;
	CLR_FLAG(PRINT_TIME);
	Config.Ptk = PI2N;

	for (int i = 0; i < sizeof(pi2data) / sizeof(pi2data[0]); i++) {
		const char* pdata = pi2data[i];
		uint64 n = atoint64(pdata, 0);
		pdata += 6;
		if (n < atoint64("e3", 0)) {
			continue;
		}
		for (int j = 1; j < 7; j ++) {
			Config.Kgap = 2 * j;
			while (isspace(*pdata))
				pdata++;

			uint64 ret = atoint64(pdata, 0);
			uint64 cal = ktprime(0, n, 0);
			if (cal != ret) {
				printf("error: PI2_%d(%lld), cal = %lld, table ret = %lld\n",
						Config.Kgap, n, cal, ret);
			}
			while (isdigit(*pdata))
				pdata++;
		}
		putchar('\n');
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

static int getSystemInfo()
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

//print the Ptk info
static void printInfo()
{
	puts("---------------------------------------------------------------");
	puts("---------------------------------------------------------------");

	printf("\
	1.%s PI(n)\n \
	2.Twin/Cousin/Sexy/p, p+2n/ prime pairs PI2(n)\n \
	3.Ktuplet prime PIk(n)\n (n < %s) version %s\n",
	KtupletName[0], KtupletName[1], MAXN, KVERSION);

	puts("Copyright (c) by Huang Yuanbing 2008 - 2013 bailuzhou@163.com");

#ifdef _MSC_VER
	printf("Compiled by MS/vc++ %d", _MSC_VER);
#else
	printf("Compiled by g++ %d.%d.%d", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#endif

#if _M_AMD64 || __x86_64__
	printf(" on x64");
#endif
	printf(" on %s %s\n", __TIME__, __DATE__);

	getSystemInfo();
	getCpuInfo();

	printf("Work threads = %d, POPCNT = %d, ASM_X86 = %d\n",
			Config.Threads, POPCNT, ASM_X86);
	printf("L1 Size = %d k, OPT_L1CACHE = %d, BIT BSHIFT = %d\n",
			 Config.CpuL1Size >> 13, OPT_L1CACHE, BSHIFT);
	puts("---------------------------------------------------------------");
	puts("---------------------------------------------------------------\n");
}

static bool setKpattern(const char* pstr)
{
	//user defined pattern
	if (pstr[1] != 'k') {
		int ns = KData.Kpattern[0] = pstr[1] - '0';
		pstr = strstr(Kpattern, pstr);
		if (!pstr)
			return false;

		pstr += 3;
		for (int j = 1; j < ns; j++) {
			int dig = 0;
			while (!isdigit(*pstr)) {
				pstr++;
			}
			while (isdigit(*pstr)) {
				dig = 10 * dig + *pstr++ - '0';
			}
			KData.Kpattern[j] = dig;
		}
		KData.Kpattern[0] = ns - 1;
		KData.Kpattern[ns] = 0;
		return true;
	}

	for (int j = 2, k = 1; pstr[j]; j++) {
		char dig = pstr[j] - '0';
		if (dig <= 0 || dig > 9) {
			continue;
		}
		int kdata = 0;
		if (dig % 2 == 0 && j < 8) {
			kdata = dig;
		} else {
			kdata = dig * 10 + pstr[++j] - '0';
		}
		assert(kdata % 2 == 0);
		if (CHECK_FLAG(PRINT_LOG)) {
			printf("%d\n", kdata);
		}
		KData.Kpattern[k] = kdata;
		KData.Kpattern[0] = k++;
	}

	return true;
}

static void printKpattern()
{
	Config.Kgap = KData.Kpattern[1];
	Config.Ptk = KData.Kpattern[0];

	int ktuplets = Config.Ptk;
	if (ktuplets < PIKN) {
		puts(KtupletName[ktuplets]);
	} else {
		if (ktuplets <= sizeof(KtupletName) / sizeof(KtupletName[0]))
			printf(KtupletName[ktuplets]);
		else
			printf("%d-tuples", ktuplets);
		printf(": p");
		for (int i = 1; i < ktuplets; i++) {
			int kp = KData.Kpattern[i];
			printf(", p + %d", kp);
			if (kp % 2 != 0 || (i > 1 && kp <= KData.Kpattern[i - 1])) {
				printf(" : invalid pattern[i] %d\n", i, kp);
			}
		}
		putchar('\n');
	}
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
static int parseCmd(char cmdparams[][80])
{
	int cmdi = -1;

	for (int i = 0; cmdparams[i][0]; i++) {
		char c = cmdparams[i][0];
		int tmp = atoi(cmdparams[i] + 1);
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
					KData.Wheel = tmp;
				}
				break;
			case 'K':
				if (cmdparams[i][1] && setKpattern(cmdparams[i])) {
					printKpattern();
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

//split cmd to cmdparams[] by white space and ';'
static int splitCmd(const char* ccmd, char cmdparams[][80])
{
	int ncmds = 0;

	for (int i = 0; ; i++) {
		while (isspace(*ccmd)) {
			ccmd++;
		}
		if (*ccmd == 0 || *ccmd == ';') {
			break;
		}
		char* pc = cmdparams[i];
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
		char cmdparams[8][80] = {0};

		if (splitCmd(cmd, cmdparams) <= 0) {
			return false;
		}

		int cmdi = parseCmd(cmdparams);
		if (cmdi < 0) {
			return true;
		}

		char cmdc = toupper(cmdparams[cmdi][0]);
		if (cmdc == 'H') {
			if (cmdparams[0][1] == 'k') {
				puts(Kpattern);
			} else {
				puts(HelpConfig);
				puts(HelpCmd);
				puts(HelpUse);
			}
		} else if (cmdc == 'B') {
			puts("----------- start benchmark ------------");
//			if (isdigit(cmdparams[cmdi + 1][0]))
//				excuteCmd("e15 0 1500 m6 d");
			benchMark(cmdparams);
		} else if (cmdc == 'U') {
			if (cmdparams[cmdi + 1][0] == 0) {
				testPi2();
				testPik();
			}
			for (int i = 0; cmdparams[cmdi + 2][i]; i++) {
				Config.Ptk = cmdparams[cmdi + 2][i] - '0';
				bool rwflag = true;
				if (Config.Ptk < 3 && Config.Ptk >= 0) {
					startTest(atoint64(cmdparams[cmdi + 1]), rwflag);
				}
			}
		} else if (cmdc == 'L') {
			puts("-----start list PI2(n)/PI(n) ------");
			listDiffGpt(cmdparams, cmdi);
		} else if (cmdc == 'I') {
			puts("-------list pow PI2(n)/PI(n) ------");
			listPowGpt(cmdparams, cmdi);
		} else if (cmdc == 'N') {
			uint64 n = atoint64(cmdparams[cmdi + 1], 1000000000);
			initParam(0, n);
			printf("patterns %lld = %d\n", n, KData.Patterns);
			KData.Wheel = 0;
		} else if (cmdc == 'E' || isdigit(cmdc)) {
			uint64 n = atoint64(cmdparams[cmdi], 1000000000);
			int pattern = atoint64(cmdparams[cmdi + 1], 0);
			ktprime(0, n, pattern);
		} else if (cmdc == 'Q') {
			return false;
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
static void initCache()
{
	initBitTable();
	simpleEratoSieve(10000);
	getSystemInfo();
	getCpuInfo();
}

int main(int argc, char* argv[])
{
	initCache();

	if (argc < 2) {
		printInfo();
	}

	for (int i = 1; i < argc; i++) {
		if (argv[i][0] == 'm')
			doCompile();
		else
			excuteCmd(argv[i]);
	}

	excuteCmd("k1 e10");
	excuteCmd("k21 e10");

	char ccmd[256] = {0};
	while (true) {
		printf("\n>> ");
		if (!gets(ccmd) || !excuteCmd(ccmd))
			break;
	}

	return 0;
}

/*************************************************************************
--------------------------------------------------------------------------
BLOCKSIZE  PI        PI2       PI3      PI4      PI5     PI6     PI7     PI8
------------------------------------------------------------------------------
210        48        15        8        3        2       1       1       1
(SP = 7)   22.86     14.29     11.43    5.71     2.86    2.85    3.33    3.84
------------------------------------------------------------------------------
2310       480       135       64       21       12      5       4       3
(SP = 11)  20.78     11.69     8.31     3.64     2.60    1.30    1.21    1.04
------------------------------------------------------------------------------
30030      5760      1485      640      189      96      35      24      18
(SP = 13)  19.18     9.98      6.39     2.52     1.60    0.70    0.56    0.48
------------------------------------------------------------------------------
510510     92160     22275     8960     2457     1152    385     240     162
(SP = 17)  18.05     8.73      5.27     1.93     1.13    0.45    0.33    0.25
------------------------------------------------------------------------------
9699690    1658880   378675    143360   36855    16128   5005    2880    1782
(SP = 19)  17.10     7.81      4.43     1.52     0.83    0.31    0.21    0.15
------------------------------------------------------------------------------
223092870  36495360  7952175   2867200  700245   290304  85085   46080   26730
(SP = 23)  16.36     7.13      3.86     1.26     0.65    0.23    0.14    0.10
--------------------------------------------------------------------------------
6469693230 1021870080,214708725 74547200 17506125 6967296 1956955 1013760 561330
(SP = 29)  15.80     6.64      3.46     1.08     0.54    0.18    0.11    0.07
---------------------------------------------------------------------------------
200560490130         6226553025,2087321600,472665375,181149696,48923875, 24330240
(SP = 31)            6.21      3.12     0.94     0.45    0.15    0.08    0.04
---------------------------------------------------------------------------------
********************************************************************************

PI(1e13) = 346065536839, time use 1681.128 s

feature:
  14. win32 gui
  15. remove unused marco

Linux g++:
  g++ -Wall -msse4 -O3 -march=native -s -pipe -ffunction-sections -fomit-frame-pointer -lpthread ktprime.cpp
Mingw/g++:
  g++ -Wall -mpopcnt -mtune=native -O2 -s -pipe -fomit-frame-pointer ktprime.cpp
MS vc++:
  cl /O2 /Os ktprime.cpp

command:
	D t2 C2000 M5 e15

c1600 m5 d t4 e15 0 1000
need 28h amd phoenm x4 830
 ****************************************************/

