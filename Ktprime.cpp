/************************************************************
copyright (C) 2007-2013 by Huang Yuanbing
mail to: bailuzhou@163.com
free use for non-commercial purposes

benchmark on ms windows
                    G(1e10)  PI2(1e10)
2.00G AMD   3600+   1.08     1.50  seconds
2.80G AMD   X4 820  0.41     0.56  seconds
2.80G AMD   X4 641  0.41     0.56  seconds
2.90G Intel E7500   0.75     1.14  seconds
1.66G Intel T5500   1.40     2.15  seconds
2.26G Intel I3 350M 0.85     1.25  seconds
2.96G Intel I5 560M 0.66     0.96  seconds

Patterns
p_i = p_1 + b_i, i = 1, 2, ..., k
s   = p_k - p_1
k   s : b
2   2 : 0  2
3   6 : 0  4  6
3   6 : 0  2  6
4   8 : 0  2  6  8
5  12 : 0  4  6  10  12
5  12 : 0  2  6  8   12
6  16 : 0  4  6  10  12  16

    * pi(x) : prime-counting function
    * pi2(x) : count of twin primes.
      twin prime = (p, p + 2)
    * pi3(x) : count of prime triplets.
      prime triplet = (p, p + 2, p + 6) or (p, p + 4, p + 6)
    * pi4(x) : count of prime quadruplets.
      prime quadruplet = (p, p + 2, p + 6, p + 8)
    * pi5(x) : count of prime quintuplets.
      prime quintuplet = (p, p + 2, p + 6, p + 8, p + 12) or (p, p + 4, p + 6, p + 10, p + 12)
    * pi6(x) : count of prime sextuplets.
      prime sextuplet = (p, p + 4, p + 6, p + 10, p + 12, p + 16)
    * pi7(x) : count of prime septuplets.
      prime septuplet = (p, p + 2, p + 6, p + 8, p + 12, p + 18, p + 20)

Goldbach partitions
Goldbach     partitions   r(n)  is  the  number     of
representations of an even number n as the sum of two primes.
A pair of primes that sum to an even integer are known as
a Goldbach partitions (Oliveira e Silva).
Letting denote the number of Goldbach partitions of
without regard to order, then the number of ways of writing
as a sum of two prime numbers taking the order of the two primes into account is

Let P be the set of primes. For an even number n ≥ 6, let n = p + q
with p, q ∈ P be a Goldbach partitions of n. Denote g(n) the number of the
unordered Goldbach partitions of n. Denote, furthermore, as in [11], N2(n)
the number of such partitions with the taking into account the order of
parts. Then, evidently,

Goldbach's comet is the name given to a plot of the function g(E),
the so-called Goldbach function
performance and throughput improvement by minimizing cache conflicts and misses
 in the last level caches of multi-cores processors.

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

# include <stdio.h>
# include <string.h>
# include <stdlib.h>
# include <ctype.h>
# include <time.h>
# include <memory.h>
# include <math.h>

# include <assert.h>

# define KVERSION       "9.4"
# define TABLE_GAP      "1e11"
# define MAXN           "1e16"
# define MINN           10000000

# define MAX_L1SIZE     (64 << 13)
# define SEGMENT_SIZE   (510510 * 4)
# define MAX_THREADS    32

//SSE4 popcnt instruction, make sure your cpu support it
//use of the SSE4.2/ SSE4a POPCNT instruction for fast bit counting.
#if _MSC_VER > 1400
  # define POPCNT      1
	# include <intrin.h>
#elif (__GNUC__ * 10 + __GNUC_MINOR__ > 44)
	# define POPCNT      0
	# include <popcntintrin.h>
#else
	# define POPCNT      0
#endif
	# define TREE2       1

#ifdef _MSC_VER
	# pragma warning(disable: 4996 4244 4127 4505 4018)
	#if _MSC_VER > 1200
	# pragma warning (disable:6328 6031)
	#endif
#endif

//#pragma pack (16)

# define OMP             0
# if OMP
	#include <omp.h>
# endif

# define FAST_CROSS      1
# define OPT_L1CACHE     1

#if defined _M_AMD64
	# define ASM_X86     0
#elif _MSC_VER >= 1200
	# define ASM_X86     1
#else
	# define ASM_X86     0
#endif

typedef unsigned char  uchar;
typedef unsigned short ushort;
typedef unsigned int   uint;

#ifdef _WIN32
	typedef unsigned __int64 uint64;
	typedef __int64 int64;
	#define CONSOLE "CON"
	#include <windows.h>
#else
	typedef unsigned long long uint64;
	typedef long long int64;
	#define CONSOLE "/dev/tty"
	#include <unistd.h>
	#include <sys/time.h>
	#include <pthread.h>
#endif

//amd 5, intel 5
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
# elif BSHIFT ==  6
	typedef uint64 utype;
	# define MASK 63
# endif

typedef uint64 stype;
# define SMOVE 6

# define MASK_N(n)         (1 << ((n) & MASK))
# define SET_BIT(a, n)     a[(n) >> BSHIFT] |= MASK_N(n)
# define FLP_BIT (a, n)    a[(n) >> BSHIFT] ^= MASK_N(n)
# define CLR_BIT(a, n)     a[(n) >> BSHIFT] &= ~MASK_N(n)
# define TST_BIT(a, n)     (a[(n) >> BSHIFT] & MASK_N(n))
# define TST_BIT2(a, n)    TST_BIT(a, (n) / 2)

#if !defined __CPLUSPLUS
	#define bool   int
	#define true   1
	#define false  0
#endif

static const char* const HelpConfig = "\
	[P: Print time use]\n\
	[D: Debug log]\n\
	[S: Save result to file]\n\
	[R: Runtime check pattern]\n\
	[A: Save/Continue last task]\n\
	[K: Calculate of Goldbach/Prime/Twin[Cousin]/Ktuplet Prime k(0 - 8)]\n\
	[M: Monitor progress m(0 - 30)]\n\
	[F: factorial of whell prime factor f(7 - 29)]\n\
	[T: Threads number t(2 - 64)]\n\
	[C: Cpu L1/L2 data cache size (L1:16-128, L2:128-1024)]\n";

static const char* const HelpCmd = "\
	[H: Help cmd [k]]\n\
	[B: Benchmark (start) (end) (gptk)]\n\
	[Q: Exit programming]\n\
	[U: Unit test (n 1 - 10000) (gptk 0 - 2)]\n\
	[N: Number of patterns (start) (count)]\n\
	[O: Optimaze best threads (count)]\n\
	[I: List base pow index (powbase) (start) (end)]\n\
	[L: List multi gptk (start) (end/count) (step)]\n";

static const char* const HelpUse = "\n\
	All command/config as follow:\n\
	B, B e9 e10 0123\n\
	C31, C128000\n\
	K41, K52, KK26, KK268 T2-32\n\
	H, Hk, A, D, S, R\n\
	U, U 1000 012, U 1000+2 2\n\
	N 2e8+20 100\n\
	M0-30, O E11 1-31\n\
	N 120000*1000 100\n\
	K0 2^31 012, K2 2e10*3\n\
	K1 10^11*2 0-3, K32 400000000+100\n\
	L 2e9-100 1000 10\n\
	I 2 10 20\n\
	L e9-100 2e9*2 1e8+2";

static const char* const Kpattern = "\n\
	k01 0\n\
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
	"Goldbach partition",
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

static const char* const TaskDataFormat =
"[Task]\n\
Gptk = %d\n\
Kgap = %d\n\
Wheel = %d\n\
Patterns = %d\n\
Tasks = %d\n\
Pbegi = %d\n\
Pendi = %d\n\
N = %lld\n\
Result = %lld";

static const char* const ConfigDataFormat =
	"[Cpu]\nWork thread = %d\nCpu L2Size = %d\nCmd = %s\n";

static const char* const PrintFormat[] =
{
#if _MSC_VER == 1200
	"G(%I64d) = %I64d",
	"PI(%I64d) = %I64d",
	"PI2(%I64d) = %I64d",
	"PI%d(%I64d) = %I64d"
#else
	"G(%lld) = %lld",
	"PI(%lld) = %lld",
	"PI2(%lld) = %lld",
	"PI%d(%lld) = %lld"
#endif
};

enum CAMODE
{
	GP0N = 0, //goldbach partition
	PI1N = 1, //prime number
	PI2N = 2, //twin
	PIKN = 3, //ktuplet prime
};

/************************************/
# define PRIME_NUMS 5761455 + 160
//prime difference in [0, 10^8]
static uchar Prime[PRIME_NUMS];

//the smallest Startp[i] * wheel % Prime[i] = 1
static uint Startp[PRIME_NUMS];

typedef ushort ptype;
static ptype Pattern[47713050 + 10];

//table of ktuplet
typedef int64 ltype;
static ltype Ktable[10000];

//CrossedTpl cross out prime <= 17
static utype CrossedTpl[(SEGMENT_SIZE >> (BSHIFT + 1)) + 100];

//bit 1 left most table
static uchar LeftMostBit1[1 << 16];

//number of bits 1 binary representation table in Range[0-2^16)
static uchar WordNumBit1[1 << 16];

//WordReverse[i] is equal to the bit reverse of i (i < 2^16)
static ushort WordReverse[1 << 16];

//the first 10 even prime numbers
static const uchar SmallPrime[ ] =
{
	3, 5, 7, 11, 13,
	17, 19, 23, 29, 31
};

//config
static struct
{
	//show result
	bool ShowResult;
	//show calculating time
	bool ShowTime;
	//show debug log
	bool ShowLog;
	//save result to file
	bool SaveToFile;
	//check the pattern
	bool CheckPattern;
	//save last task
	bool SaveTask;

	//flag for G(n), PI2(n), PI(n), Pik(n)
	int Gptk;
	//ktuplet pattern gap, 2 for twin
	int Kgap;
	//cpu L1/L2 size
	int CpuL1Size;
	int CpuL2Size;
	//work threads
	int Threads;
	//print the progress gap
	int PrintGap;
}
Config =
{
	true, true, false,
	false, false, false,
	0, 2, MAX_L1SIZE, 1100, 4, (1 << 9) - 1
};

static struct
{
	int Wheel;
	int SieveCache;

	int firstSievedIndex;

	int SqrtN;
	int Patterns;
	int Factorial;
	int PatternDiff;

	bool UseKtable;
	int Kpattern[16];
	ltype N;
}
KData =
{
	0, 0, 8, 0, 0,
	0, 0, 0, {4, 2, 6, 8},
	100000
};

static struct TaskConfig
{
	int Gptk;
	int Kgap;

	int Wheel;
	int Patterns;

	int Tasks;

	int Pbegi;
	int Pendi;

	ltype N;
	ltype Result;
}
LastTask =
{
	0, -1, -1, 0, 4
};

static struct ThreadInfo
{
	int Pbegi;
	int Pendi;
	ltype Result;
} TData[MAX_THREADS];

static ltype sievePattern(int, int);

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

static ltype startWorkThread(int threads, int pbegi, int pendi)
{
	ltype gpts = 0;
	int i;
	//assert(pendi >= pbegi && pbegi >= 0);
	if (threads > MAX_THREADS) {
		threads = 4;
	}
	if (pendi - pbegi < threads) {
		threads = 1;
	}

	devideTaskData(threads, pbegi, pendi);

#ifdef _WIN32
	HANDLE thandle[MAX_THREADS];
	DWORD tid[MAX_THREADS];
	for (i = 0; i < threads; i++) {
		thandle[i] = CreateThread(NULL, 0, threadProc,
			(LPVOID)(&TData[i]), 0, &tid[i]);
		Sleep(5);

		if (thandle[i] == NULL) {
			printf("create win32 thread error %ld\n", GetLastError( ));
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
static double getTime( )
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
        return {y, x-y*(a div b)}
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

//       n = 0     1
//x^n =  n%2 = 1   x * x*(n-1)
//       n%2 = 0   (x^n/2)^2
static ltype ipow(ltype x, uint n)
{
	ltype result = 1;
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

//convert str to ltype ((x)*E(y)+-*(z))
//invalid input format: 123456 1234-12 e9 2e7+2^30 2e10-2 10^11-25 2e6*2
static ltype atoint64(const char* str, ltype defaultn = 0)
{
	ltype n = 0;

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
	stype* memstart = (stype*)bitarray + (bitleng >> SMOVE);
	memstart[0] |= ~(((stype)1 << (bitleng % (1 << SMOVE))) - 1);
	memstart[1] = (stype)(~0);
}

//use asm to accelerate hot spot
#if _MSC_VER && ASM_X86
	#define EBP_OFFSET	32
__declspec(naked)
#endif
static void
set2BitArray(utype bitarray[], int s1, int s2, const int step)
{
#if ASM_X86 == 0 || _MSC_VER == 0
	if (s2 > s1) {
		s2 ^= (s1 ^= (s2 ^= s1));
	}

#if 0
	if (s1 & MASK == s2 & MASK) {
		for (; s2 > 0; ) {
			const utype mask = MASK_N(s1);
			bitarray[s1 >> BSHIFT] |= mask; s1 += step;
			bitarray[s2 >> BSHIFT] |= mask; s2 += step;
		}
	}
#endif

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

#if 1
static int64
set2BitArray(utype bitarray[], const int64 start, const int step, const int bitleng)
{
	int s2 = start >> 32;
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

	return ((int64)s2 << 32) | s1;
}

static int64
setBitArray0(utype bitarray[], const int64 start, const int step, const int leng)
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

	return ((int64)s2 << 32) | s1;
}

static int inline
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
#endif

//all old number in range[start, statr + leng] which
//is multiple of factor will be crossed
//out and the bitarray will the bit position 1
//the ith bit of bitarray is map to number start + 2 * i + 1
/**
static void crossOutEvenFactor(utype bitarray[], const ltype start,
						int leng, const int factor)
{
	int offset = factor - start % factor;
	if (offset % 2 == 0) {
		offset += factor;
	} else if (start <= factor) {
		offset += 2 * factor;
	}

	const int bits = leng >> 1;
	for (offset >>= 1; offset <= bits; offset += factor) {
		SET_BIT(bitarray, offset);
	}
}*/

//the ith bit of bitarray is map to start + 2 * i + 1
//it's difference with crossOutEvenFactor, only
//6k + 1, 6k + 5 number which is multiple of factor will be crossed out
//and gain performance 1/3 improvement
static void crossOutFactor(utype bitarray[], const ltype start,
						const int leng, int factor)
{
	int s1 = factor - start % factor;
	if (s1 % 2 == 0) {
		s1 += factor;
	} else if (start <= factor) {
		s1 += 2 * factor;
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
static void sieveWheelFactor(utype bitarray[], const ltype start, const int leng, const uint wheel)
{
	assert(wheel % 6 == 0);

	for (int i = 2, p = 3; wheel % p == 0; p += Prime[++i]) {
		crossOutFactor(bitarray, start, leng, p);
		if (start <= p) {
			SET_BIT(bitarray, (p - start) / 2);
		}
	}
}

//sieve prime in [start, start + leng]
static void segmentedEratoSieve1(utype bitarray[], const ltype start, const int leng, bool initzero)
{
	const int sqrtn = (int)(sqrt((double)start + leng) + 0.1) + 1;

	if (initzero)
		memset(bitarray, 0, (leng >> 4) + 1);

	for (int i = 2, p = 3; p < sqrtn; p += Prime[++i]) {
		crossOutFactor(bitarray, start, leng, p);
	}
	if (start == 0) {
		*(ushort*)bitarray = 0x3491;
	}
}

//another algorithm for sieve prime in [start, start + leng]
static void segmentedEratoSieve2(utype bitarray[], const ltype start, const int leng, bool initzero)
{
	const int sqrtn = (int)(sqrt((double)start + leng) + 0.1) + 1;

	assert(start % SEGMENT_SIZE == 0 && SEGMENT_SIZE % 19);

	if (initzero) {
		memcpy(bitarray, CrossedTpl, (leng >> 4) + 1);
	} else {
	//		bitarray[i] |= CrossedTpl[i];
	}

	for (int i = 8, p = 19; p < sqrtn; p += Prime[++i]) {
		crossOutFactor(bitarray, start, leng, p);
	}

	if (start == 0) {
		//the first 7th bit poisition set 1
		*(ushort*)bitarray = 0x3491;
	}
}

//reverse bit order of a byte with binary representation
static uchar reverseByte(const uchar c)
{
	uchar n =
		(c & 0x55) << 1 | (c & 0xAA) >> 1;
	n = (n & 0x33) << 2 | (n & 0xCC) >> 2;
	n = (n & 0x0F) << 4 | (n & 0xF0) >> 4;
	return n;
}

static inline int countBit1s(stype n)
{
#if POPCNT
	//popcnt instruction : INTEL i7/SSE4.2, AMD Phonem/SSE4A
	#if _M_AMD64 || __x86_64__
		return _mm_popcnt_u64(n);
	#elif (SMOVE == 5)
		return _mm_popcnt_u32(n);
	#else
		return _mm_popcnt_u32(n) + _mm_popcnt_u32(n >> 32);
	#endif
#elif TREE2 == 0
	#if SMOVE == 5
		return //WordNumBit1[(ushort)n] + WordNumBit1[n >> 16];
				WordNumBit1[n & 0xffff] + WordNumBit1[n >> 16];
	#else
		uint hig = n >> 32, low = (uint)n;
		return WordNumBit1[(ushort)low] + WordNumBit1[low >> 16] +
				WordNumBit1[(ushort)hig] + WordNumBit1[hig >> 16];
	#endif
#else
	#if SMOVE == 6
		n -= (n >> 1) & 0x5555555555555555ull;
		n = (n & 0x3333333333333333ull) + ((n >> 2) & 0x3333333333333333ull);
		n = (n + (n >> 4)) & 0x0F0F0F0F0F0F0F0Full;
		n += n >> 8;
		n += n >> 16;
		n += n >> 32;
		return (n & 0x00000000FF);
	#else
		n -= (n >> 1) & 0x55555555;
		n = (n & 0x33333333) + ((n >> 2) & 0x33333333);
		n = (n + (n >> 4)) & 0x0F0F0F0F;
		n += n >> 8;
		n += n >> 16;
		return (n & 0x0000003F);
	#endif
#endif
}

//count number of bit 0 in binary representation
//!!! buffer of bitarray after position bitleng packeked with bit 1
static int countZeroBitsArray(utype bitarray[], const int bitleng)
{
	int bit1s = 0;
	int loops = bitleng >> SMOVE;

	packQwordBit1(bitarray, bitleng);
	for (stype* psbuf = (stype*) bitarray; loops >= 0; loops--) {
		bit1s += countBit1s(*psbuf++);
	}

	return ((1 + (bitleng >> SMOVE)) << SMOVE) - bit1s;
}

//reverse word array bitarray with length = bitleng
static void reverseByteArray(ushort bitarray[], const int bitleng)
{
	assert(bitleng % 8 == 0);
	ushort* ps = bitarray;
	ushort* pe = (ushort*)((uchar*)ps + bitleng / 8 - 2);

	while (ps < pe) {
		ushort tmp = WordReverse[*ps];
		*ps++ = WordReverse[*pe];
		*pe-- = tmp;
	}

	if (ps == pe) {
		*ps = WordReverse[*ps];
	} else if ((uchar*)pe + 1 == (uchar*)ps) {
		*((uchar*)ps) = WordReverse[*ps] >> 8;
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

/**
static inline int test_and_set_bit(int nr, volatile void *addr)
{
	int oldbit;
	asm volatile("bts %2,%1\n\t"
			"sbb %0,%0"
			: "=r" (oldbit), ADDR
			: "Ir" (nr) : "memory");
	return oldbit;
}

static inline int test_and_clear_bit(int nr, volatile void *addr)
{
	int oldbit;
	asm volatile("btr %2,%1\n\t"
			"sbb %0,%0"
			: "=r" (oldbit), ADDR
			: "Ir" (nr) : "memory");
	return oldbit;
}*/

//make sure no divide overflow
//improvement of 100%
static inline int
asmMulDiv(const uint startp, const uint pattern, uint p)
{
#ifdef LINE_COVER
	p = ((ltype)startp) * pattern % p;
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
	p = (((ltype)startp) * pattern - bitleng) % p + bitleng;
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
Prime[i] is the difference of the adjacent prime:
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

	if (Config.ShowLog) {
		printf("Prime[%d] = %d\n", primes, lastprime);
	}

	return primes;
}

//init bit tables
static void initBitTable( )
{
	//1. init WordNumBit1 table in 0-2^16, can use popcnt replace it
	int nsize = sizeof(WordNumBit1) / sizeof(WordNumBit1[0]);
	int i;

#if 0 == POPCNT
	WordNumBit1[0] = 0;
	for (i = 1; i < nsize; i++) {
		WordNumBit1[i] = WordNumBit1[i >> 1] + (i & 1);
	}
#endif

	//2. init bit WordReverse table
	uchar bytereverse[256] = {0};
	nsize = sizeof(WordReverse) / sizeof(WordReverse[0]);
	//reverse bit order of byte(with 8 bit) in [0, 2^8)
	for (i = 1; i < (1 << 8); i++) {
		bytereverse[i] = reverseByte((uchar)i);
	}
	//reverse bit order of short(with 16 bit) in [0, 2^16)
	for (i = 1; i < nsize; i++) {
		WordReverse[i] = bytereverse[i >> 8] | (bytereverse[i & 255] << 8);
	}

	//3. init LeftMostBit1 table
	for (int m = 2; m < (1 << 16); m += 2) {
		LeftMostBit1[m + 0] = LeftMostBit1[m >> 1] + 1;
		LeftMostBit1[m + 1] = 0;
	}

#if 0 == FAST_CROSS
	//4. init CrossedTpl table, pre sieve the factor in array sievefactor
	sieveWheelFactor(CrossedTpl, 0, sizeof(CrossedTpl) * 16, SEGMENT_SIZE);
#endif
}

//
static int savePrimeDiff(const utype bitarray[], const int start,
						const int bitleng, uchar primediff[])
{
	int primes = 0;
	static int lastprime = 2;
	if (start == 0) {
		lastprime = 2;
	}

	for (int p = 1; p < bitleng; p += 2) {
		if (!TST_BIT2(bitarray, p)) {
			primediff[primes++] = start + p - lastprime;
			lastprime = start + p;
		}
	}

	return primes;
}

static int savePi1Pattern(const utype bitarray[], int start, int leng,
		int& lastpattern, ptype pi1pattern[])
{
	int pi1n = 0;

	for (int p = 1; p < leng; p += 2) {
		if (!TST_BIT2(bitarray, p)) {
			if (pi1pattern) {
				pi1pattern[pi1n] = p + start - lastpattern;
				lastpattern = p + start;
/*				if (Config.CheckPattern) {
					assert(gcd(lastpattern, KData.Wheel) == 1);
				} */
//				printf("p[%d] = %d\n", pi1n, lastpattern);
			}
			pi1n++;
		}
	}

	return pi1n;
}

static int savePi2Pattern(const utype bitarray[], int start, int leng,
		int& lastpattern, ptype pi2pattern[])
{
	int pi2n = 0;

	for (int p = 1; p < leng; p += 2) {
		if (!TST_BIT2(bitarray, p) && !TST_BIT2(bitarray, p + Config.Kgap)) {
			if (pi2pattern) {
				pi2pattern[pi2n] = p + start - lastpattern;
				lastpattern = p + start;
				//printf("p2[%d] = %d\n", pi2n, lastpattern);
			}
			pi2n++;
		}
	}

	return pi2n;
}

static int savePikPattern(const utype bitarray[], int start, int leng,
		int& lastpattern, ptype pikpattern[])
{
	int pikn = 0;
	const int pis = KData.Kpattern[0];
	const int maxpd = 1 << (8 * sizeof(ptype));

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
			if (pikpattern) {
				assert(p - lastpattern < maxpd);
				pikpattern[pikn] = p + start - lastpattern;
				lastpattern = p + start;
				//printf("pk[%d] = %d\n", pikn, lastpattern);
			}
			pikn++;
		}
	}

	return pikn++;
}

//
static int getPi2Pattern(const int factorial, ptype pi2pattern[])
{
	int lastpattern = 0;
	int sleng = SEGMENT_SIZE, patterns = 0;

	for (int start = 0; start < factorial; start += sleng) {
		utype bitarray[(SEGMENT_SIZE + 1000) >> (BSHIFT + 1)];
		if (start + sleng >= factorial)
			sleng = (int)(factorial - start + 1);

		memset(bitarray, 0, (sleng >> 4) + 1 + (LastTask.Kgap >> 4));
		//memset(bitarray, 0, (sleng >> 4) + 1);

		sieveWheelFactor(bitarray, start, sleng + Config.Kgap, KData.Wheel);
		if (pi2pattern)
			patterns += savePi2Pattern(bitarray, start, sleng,
						lastpattern, pi2pattern + patterns);
		else
			patterns += savePi2Pattern(bitarray, start, sleng, lastpattern, 0);
	}

	if (pi2pattern) {
		if (Config.ShowLog) {
			printf("factorial pattern pi2n = %d\n", patterns);
		}
		assert(patterns < sizeof(Pattern));
		KData.PatternDiff = factorial - lastpattern;
		pi2pattern[patterns] = 0;
		patterns *= KData.Wheel / factorial;
	}

	return patterns;
}

//
static int getPi1Pattern(const int factorial, ptype pi1pattern[])
{
	int lastpattern = 0;
	int sleng = SEGMENT_SIZE, patterns = 0;

	for (int start = 0; start < factorial; start += sleng) {
		utype bitarray[(SEGMENT_SIZE + 1000) >> (BSHIFT + 1)];
		if (start + sleng >= factorial)
			sleng = (int)(factorial - start + 1);

		memset(bitarray, 0, (sleng >> 4) + 1);

		sieveWheelFactor(bitarray, start, sleng, KData.Wheel);
		if (pi1pattern)
			patterns += savePi1Pattern(bitarray, start, sleng,
						lastpattern, pi1pattern + patterns);
		else
			patterns += savePi1Pattern(bitarray, start, sleng, lastpattern, 0);
	}

	if (pi1pattern) {
		if (Config.ShowLog) {
			printf("factorial pattern pi1n = %d\n", patterns);
		}
		KData.PatternDiff = factorial - lastpattern;
		pi1pattern[patterns] = 0;
		patterns *= KData.Wheel / factorial;
	}

	return patterns;
}

//
static int getPikPattern(const int factorial, ptype pikpattern[])
{
	int lastpattern = 0;
	int sleng = SEGMENT_SIZE, patterns = 0;

	for (int start = 0; start < factorial; start += sleng) {
		utype bitarray[(SEGMENT_SIZE + 1000) >> (BSHIFT + 1)];
		if (start + sleng >= factorial)
			sleng = (int)(factorial - start + 1);

		memset(bitarray, 0, (sleng >> 4) + 1 + (LastTask.Kgap >> 4));

		sieveWheelFactor(bitarray, start, sleng + LastTask.Kgap, KData.Wheel);
		if (pikpattern)
			patterns += savePikPattern(bitarray, start, sleng,
						lastpattern, pikpattern + patterns);
		else
			patterns += savePikPattern(bitarray, start, sleng, lastpattern, 0);
	}

	if (pikpattern) {
		if (Config.ShowLog) {
			printf("factorial pattern pi%d = %d\n", KData.Kpattern[0], patterns);
		}
		KData.PatternDiff = factorial - lastpattern;
		pikpattern[patterns] = 0;
		patterns *= KData.Wheel / factorial;
	}

	return patterns;
}

//min prime in range [start, sum / 2]
static int getSegpattern(const int sum, int start, ptype pattern[])
{
	int gpn = 0;
	int lastpattern = start;
	const int bitleng = sum / 2 + ((sum / 2) & 1);

	for (int sleng = SEGMENT_SIZE; start < bitleng; start += sleng) {
		utype bitarray[(SEGMENT_SIZE + 1000) >> (BSHIFT + 1)];
		if (sleng >= bitleng - start) {
			sleng = bitleng - start;
		}

		memset(bitarray, 0, (sleng >> 4) + 1);

		sieveWheelFactor(bitarray, sum - start - sleng, sleng + 16, KData.Wheel);
		reverseBitArray(bitarray, sleng >> 1);
		sieveWheelFactor(bitarray, start, sleng + 16, KData.Wheel);

		gpn += savePi1Pattern(bitarray, start, sleng, lastpattern, pattern + gpn);
	}

	return gpn;
}

//optimize for memory
static int getGppattern(const ltype n, const int factorial, ptype gppattern[])
{
	const int sumgp1 = n % factorial;
	const int sumgp2 = sumgp1 + factorial;

	int gppn = getSegpattern(sumgp1, 0, gppattern);
	gppattern[gppn++] = 0;
	KData.PatternDiff = sumgp1;
	gppn += getSegpattern(sumgp2, sumgp1, gppattern + gppn);
	gppn --;

	if (Config.ShowLog) {
		printf("factorial pattern gppn = %d\n", gppn);
	}

	return gppn;
}

//
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
	int k = KData.firstSievedIndex;
	int spos[MAX_L1SIZE / 11];

	for (; p <= minp; p += Prime[++k]) {
		spos[k] = asmMulDiv(Startp[k], pattern, p);
	}
	//		assert(k < sizeof(spos));
	for (int start = 0, sleng = Config.CpuL1Size; start < bitleng; start += Config.CpuL1Size) {

		k = KData.firstSievedIndex;
		p = SmallPrime[k - 2];

		if (start + Config.CpuL1Size > bitleng) {
			sleng = bitleng - start;
		}

		for (; p <= minp; p += Prime[++k]) {
			int npos = spos[k];
			if (npos < sleng) {
				npos = setBitArray(bitarray + (start >> BSHIFT), npos, p, sleng);
			}
			spos[k] = npos - sleng;
		}
	}

	return k;
}

static int sieveGPnL1(utype bitarray[], const uint pattern1, int bitleng, int& p)
{
	int pattern2 = LastTask.Kgap - pattern1;
	if (pattern2 < 0) {
		pattern2 += KData.Wheel;
	}
	int k = KData.firstSievedIndex;
	const int minp = KData.SqrtN < Config.CpuL1Size ? KData.SqrtN : Config.CpuL1Size;

	//186k
	int64 spos[MAX_L1SIZE / 11];

	for (; p <= minp; p += Prime[++k]) {
		spos[k] = asmMulDiv(Startp[k], pattern1, p);
		//assert(spos[k][0] * wheel + pattern1 % p == 0)
		int s2 = bitleng - asmMulDivSub(Startp[k], pattern2, p, bitleng);
		if (s2 < 0)
			s2 += p;
		spos[k] |= (int64)s2 << 32;
	}

	for (int start = 0, sleng = Config.CpuL1Size; start < bitleng; start += Config.CpuL1Size) {

		k = KData.firstSievedIndex;
		p = SmallPrime[k - 2];

		if (start + Config.CpuL1Size > bitleng) {
			sleng = bitleng - start;
		}

		//40%
		for (; p <= minp; p += Prime[++k]) {
			spos[k] = set2BitArray(bitarray + (start >> BSHIFT), spos[k], p, sleng);
		}
	}

	return k;
}

#if 1
static int sievePi2L1(utype bitarray[], const uint pattern1, const uint pattern2, int bitleng, int& p)
{
	const int minp = KData.SqrtN < Config.CpuL1Size ? KData.SqrtN : Config.CpuL1Size;
	int k = KData.firstSievedIndex;
	int64 spos[MAX_L1SIZE / 11];

	for (; p <= minp; p += Prime[++k]) {
		const uint s1 = asmMulDiv(Startp[k], pattern1, p);
		const uint s2 = asmMulDiv(Startp[k], pattern2, p);
		spos[k] = (int64)s2 << 32 | s1;
	}

	for (int start = 0, sleng = Config.CpuL1Size; start < bitleng; start += Config.CpuL1Size) {

		k = KData.firstSievedIndex;
		p = SmallPrime[k - 2];

		if (start + Config.CpuL1Size > bitleng) {
			sleng = bitleng - start;
		}

		for (; p <= minp; p += Prime[++k]) {
			spos[k] = set2BitArray(bitarray + (start >> BSHIFT), spos[k], p, sleng);
		}
	}

	return k;
}
#else
static int sievePi2L1(utype bitarray[], const uint pattern1, const uint pattern2, int bitleng, int& p)
{
	const int minp = KData.SqrtN < Config.CpuL1Size ? KData.SqrtN : Config.CpuL1Size;
	int k = KData.firstSievedIndex;
	int sleng = Config.CpuL1Size;
	int64 spos[MAX_L1SIZE / 11];

	for (; p <= minp; p += Prime[++k]) {
		int s1 = asmMulDivSub(Startp[k], pattern1, p, bitleng);
		if (s1 > bitleng) {
			s1 -= p;
		}

		int s2 = asmMulDivSub(Startp[k], pattern2, p, bitleng);
		if (s2 > bitleng) {
			s2 -= p;
		}

		s1 -= (bitleng - Config.CpuL1Size);
		s2 -= (bitleng - Config.CpuL1Size);

		spos[k] = (int64)s2 << 32 | s1;
	}

	for (int start = bitleng; start > 0; start -= sleng) {

		k = KData.firstSievedIndex;
		p = SmallPrime[k - 2];

		if (start < sleng) {
			start = sleng;
		}
		for (; p <= minp; p += Prime[++k]) {
			spos[k] = setBitArray(bitarray + ((start - sleng) >> BSHIFT), spos[k], -p, sleng);
		}
	}

	return k;
}
#endif

static int sievePikL1(utype bitarray[], const uint pattern, int bitleng, int& p)
{
	int k = KData.firstSievedIndex;
#if 1
	int ns = KData.Kpattern[0];
	if (ns & 1) {
		k = sievePi1L1(bitarray, pattern, bitleng, p);
		ns -= 1;
	}
	for (int j = ns; j > 0; j -= 2) {
		k = KData.firstSievedIndex;
		p = SmallPrime[k - 2];
		k = sievePi2L1(bitarray, pattern + KData.Kpattern[j],
				pattern + KData.Kpattern[j - 1], bitleng, p);
	}
#else

	const int minp = KData.SqrtN < Config.CpuL1Size ? KData.SqrtN : Config.CpuL1Size;
	const int ns = KData.Kpattern[0];
	int spos[MAX_L1SIZE / 11][8];

	for (; p <= minp; p += Prime[++k]) {
		spos[k][1] = asmMulDiv(Startp[k], pattern, p);
		for (int j = 1; j < ns; j ++) {
			spos[k][j + 1] = asmMulDiv(Startp[k], KData.Kpattern[j] + pattern, p);
		}
	}

	for (int start = 0, sleng = Config.CpuL1Size; start < bitleng; start += sleng) {

		if (start + sleng > bitleng) {
			sleng = bitleng - start;
		}

		for (int j = 1; j <= ns; j++) {
			k = KData.firstSievedIndex;
			p = SmallPrime[k - 2];
			for (; p <= minp; p += Prime[++k]) {
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

	int k = KData.firstSievedIndex;
	int p = SmallPrime[k - 2];

#if OPT_L1CACHE
	// performance improvement from 100 -> 61
	if (bitleng > Config.CpuL1Size)
		k = sievePi1L1(bitarray, pattern, bitleng, p);
#endif

	for (; p <= sqrtn; p += Prime[++k]) {
#if 0
		int s1 = asmMulDiv(Startp[k], pattern, p);
		if (s1 < bitleng)
			setBitArray(bitarray, s1, p, bitleng);
#else
		int s1 = asmMulDivSub(Startp[k], pattern, p, bitleng);
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
	const int bitleng = 1 + ((KData.N - Config.Kgap - pattern) / KData.Wheel);
	const int sqrtn = KData.SqrtN;

	int k = KData.firstSievedIndex;
	int p = SmallPrime[k - 2];

#if MAX_L1SIZE
	//performance improvement from 10 -> 72
	if (bitleng > Config.CpuL1Size) {
		k = sievePi2L1(bitarray, pattern, pattern + Config.Kgap, bitleng, p);
	}
#endif

	for (; p <= sqrtn; p += Prime[++k]) {
		int s1 = asmMulDivSub(Startp[k], pattern, p, bitleng);
		if (s1 > bitleng) {
			s1 -= p;
		}
		int s2 = s1 + Startp[k] * 2;
		if (Config.Kgap > 2) {
			s2 = s1 + Startp[k] * Config.Kgap % p;
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
	int k = KData.firstSievedIndex;
	int p = SmallPrime[k - 2];
	int s[16];

#if OPT_L1CACHE
	// performance improvement from 100 -> 61
	if (bitleng > Config.CpuL1Size)
		k = sievePikL1(bitarray, pattern, bitleng, p);
#endif

	for (; p <= sqrtn; p += Prime[++k]) {
		int s1 = asmMulDivSub(Startp[k], pattern, p, bitleng);
		if (s1 > bitleng) {
			s1 -= p;
		}

		s[1] = s1;
		int ns = initFirstPos(s, Startp[k], p, bitleng);
		if (ns & 1)
			setBitArray(bitarray, s[ns--], -p);
		for (int j = ns; j > 0; j -= 2) {
			set2BitArray(bitarray, s[j], s[j - 1], -p);
		}
	}

	return bitleng;
}

//set goldbach partition (k*wheel + pattern1) with pattern = pattern1
static int sieveGpn(utype bitarray[], const int pattern1)
{
	int pattern2 = LastTask.Kgap - pattern1;
	if (pattern2 < 0) {
		pattern2 += KData.Wheel;
	}
	int bitleng = (int)((KData.N - pattern1 - pattern2) / KData.Wheel);

	const int sqrtn = KData.SqrtN;

	int k = KData.firstSievedIndex;
	int p = SmallPrime[k - 2];

#if OPT_L1CACHE
	//performance improvement from 100 -> 72
	if (bitleng > Config.CpuL1Size)
		k = sieveGPnL1(bitarray, pattern1, bitleng, p);
#endif

	for (; p <= sqrtn; p += Prime[++k]) {
		int s2 = bitleng - asmMulDiv(Startp[k], pattern2, p);
		//int s1 = bitleng - asmMulDivSub(Startp[k], pattern2, p, 0);
		int s1 = asmMulDivSub(Startp[k], pattern1, p, bitleng);
		if (s1 > bitleng) {
			s1 -= p;
		}
		set2BitArray(bitarray, s1, s2, -p);
	}

	if (pattern2 == pattern1) {
		bitleng = 1 + bitleng / 2;
	}

	return bitleng;
}

//bad performance !!!!
static int countKtable(const ushort bitarray[], const int bitleng,
		const int pattern, int64 tdata[])
{
	const int64 base = atoint64(TABLE_GAP, 0);
	const int wordleng = bitleng / 16;
	const int diff = pattern + KData.Kpattern[KData.Kpattern[0]];
	int ktupels = 0;

	for (int b = 0; b <= wordleng; b ++) {
		ushort masks = (ushort)(~bitarray[b]);
		while (masks != 0) {
			const int bitindex = LeftMostBit1[masks];
			const int64 P = (int64)KData.Wheel * (b * 16 + bitindex) + diff;
			ktupels++;
			tdata[(int)(P / base)]++;
			masks &= masks - 1;
		}
	}

	return ktupels;
}

//
static int getNextPattern(int pattern, ptype** pnext)
{
	if (**pnext != 0) {
		pattern += **pnext;
	} else if (GP0N == Config.Gptk) {
		pattern = KData.PatternDiff + *(++*pnext);
	} else {
		pattern += KData.PatternDiff + Pattern[0];
		*pnext = Pattern;
	}

	(*pnext)++;

	return pattern;
}

static void printProgress(const int tid, const double tstart, int aves, int pcnt)
{
	double currper = 100.0 * pcnt / KData.Patterns;
	double totaltime = (getTime( ) - tstart) / (10 * currper);
	static ltype lastValue = 0;
	if (pcnt <= Config.PrintGap) {
		lastValue = 0;
	}

	printf("thread(%d) %.2lf%%, ~= ", tid, currper);
	if (totaltime < 10000) {
		printf("%.2lf s", totaltime);
	} else {
		printf("%.2lf h", totaltime / 3600);
	}

	ltype curValue = ((ltype)aves) * KData.Patterns;
	printf(", %s ~= %lld", KtupletName[Config.Gptk], curValue);
	if (lastValue > 0) {
		printf(", err ~= %.4lf%%%%",
				(curValue - lastValue) * 10000.0 / lastValue);
	}
	putchar('\n');

	lastValue = curValue;
}

//thread call: get result form pattern pbegi to pendi
//calcultate wheel * k + Pattern[pbegi, pendi] <= n
static ltype sievePattern(const int pbegi, const int pendi)
{
	static int stid = 0;
	static int scnt = 0;

	if (pbegi == 0) {
		stid = scnt = 0;
	}

	int tid = ++stid;

	ltype gpts = 0;

	utype sbuffer[(256 << 13) >> BSHIFT];
	utype* bitarray = sbuffer;
	ltype *tdata = 0;
	if (KData.UseKtable) {
		tdata = (ltype*)malloc(100000 * sizeof(ltype));
		memset(tdata, 0, sizeof(tdata[0]) * 100000);
	}

	if (KData.SieveCache > sizeof(sbuffer) / sizeof(sbuffer[0])) {
		bitarray = (utype*)malloc((KData.SieveCache << (BSHIFT - 3)) + 100);
	}

	double tstart = getTime( );
	ptype* pnext = Pattern;
	int pattern = 0;
	for (int i = 0; i < pbegi; i++) {
		pattern = getNextPattern(pattern, &pnext);
	}

	if (Config.ShowLog) {
		printf("thread %d : pattern %d - %d\n",
				tid, pbegi, pendi);
	}

	for (int pcuri = pbegi; pcuri < pendi; pcuri++) {

		pattern = getNextPattern(pattern, &pnext);
#if 0
		if (Config.CheckPattern) {
			if (gcd(KData.Wheel, pattern) != 1)
				printf("error pattern = %d\n", pattern);
			continue;
		}
#endif
		memset(bitarray, 0, 8 + KData.SieveCache * sizeof(bitarray[0]));
//		memset(bitarray, 0xA0A0A0A0, 8 + KData.SieveCache * sizeof(bitarray[0]));

		int bitleng = 0;
		if (GP0N == Config.Gptk) {
			bitleng = sieveGpn(bitarray, pattern);
		} else if (PI2N == Config.Gptk) {
			bitleng = sievePi2(bitarray, pattern);
		} else if (PI1N == Config.Gptk) {
			bitleng = sievePi1(bitarray, pattern);
		} else {
			bitleng = sievePik(bitarray, pattern);
		}

		{
			bitarray[0] |= 1;
			gpts += countZeroBitsArray(bitarray, bitleng);
		}
		if (KData.UseKtable) {
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

	if (bitarray != sbuffer) {
		free(bitarray);
	}

	if (KData.UseKtable) {
		for (int i = 0; tdata[i]; i++) {
			Ktable[i] += tdata[i];
		}
		free(tdata);
	}

	if (Config.ShowLog) {
		printf("Thread %d: pattern[%3d - %3d] = %lld\n",
				tid, pbegi, pendi, gpts);
	}

	return gpts;
}

//factorial of prime factor of wheel
static int getFactorial(const int wheel)
{
	int factorial = 2;
	for (int i = 0, p = SmallPrime[i]; wheel % p == 0; p = SmallPrime[++i]) {
		factorial *= p;
	}

	return factorial;
}

//get the frist prime index in Prime which is not a factor of wheel
static int getFirstPrime(const int wheel)
{
	int sievePrimeIndex = 1;
	for (int p = 2; wheel % p == 0; p += Prime[++sievePrimeIndex]) {

	}

	return sievePrimeIndex;
}

//
static int countKpattern(int wheel, const int kgap)
{
	int patterns = 3;
	if (kgap == 1) {
		patterns = 8;
	}
	for (int i = 2, p = SmallPrime[i]; wheel % p == 0; p = SmallPrime[++i]) {
		wheel /= p;
		patterns *= p - kgap;
	}

	return patterns * (wheel / 30);
}

//16 bits number p1, p2 and gcd((p1 + p2), factorial) = 1
// ---------- gp ------------------*/
static int initPattern(const ltype n, const int wheel)
{
	int pns = KData.Patterns;
	double ts = getTime( );

	if (GP0N == Config.Gptk) {
		if (LastTask.Wheel != wheel ||
			LastTask.Kgap != n % wheel ||
			LastTask.Gptk != Config.Gptk) {
			LastTask.Kgap = n % wheel;
			pns = getGppattern(n, wheel, Pattern);
		}
	} else if (PI1N == Config.Gptk) {
		if (LastTask.Wheel != wheel ||
			LastTask.Gptk != Config.Gptk) {
			pns = getPi1Pattern(KData.Factorial, Pattern);
			if (pns > 10) {
				assert(pns == countKpattern(wheel, 1));
			}
		}
	} else if (PI2N == Config.Gptk) {
		if (LastTask.Wheel != wheel ||
			LastTask.Kgap != Config.Kgap ||
			LastTask.Gptk != Config.Gptk) {
			LastTask.Kgap = Config.Kgap;
			pns = getPi2Pattern(KData.Factorial, Pattern);
			if (Config.Kgap < 6 && pns > 10) {
				assert(pns == countKpattern(wheel, 2));
			}
		}
	} else {
		LastTask.Kgap = KData.Kpattern[KData.Kpattern[0] - 1];
		pns = getPikPattern(KData.Factorial, Pattern);
	}

	if (Config.ShowLog) {
		const int ns = Config.Gptk == GP0N ? 2 : KData.Kpattern[0];
		printf("init pattern time use %.2lf ms, count %.3lf, %.2lf fast than pi\n",
				getTime( ) - ts, 100.0 * pns * ns / wheel,
				1.0 * countKpattern(wheel, 1) / (pns * ns));
	}

	if (Config.CheckPattern) {
		printf("wheel multiple = %d\n", pns);
		if (GP0N == Config.Gptk)
			assert(pns < sizeof(Pattern));
	}

	return pns;
}

static void initStartp(const int wheel)
{
	const int sqrtn = KData.SqrtN + 10;
	if (wheel == LastTask.Wheel && Startp[1] >= sqrtn) {
		return;
	}

	double ts = getTime( );

	for (int j = 1, p = 2; Prime[j] && p < sqrtn; p += Prime[++j]) {
		int y = 0;
		y = extendedEuclid(-wheel % p, p, y);
		if (y < 0) {
			y += p;
		}
		Startp[j] = y;
		//		assert(Startp[j] < p && Startp[j] >= 0);
	}

	Startp[1] = sqrtn;

	if (Config.ShowLog) {
		printf("init startp time use %.2lf ms\n", getTime( ) - ts);
	}
}

static uint getDefaultWheel(const ltype n)
{
	uint wheel = KData.Wheel;
	if (wheel > 30 && wheel % 30 == 0 && n > wheel) {
		return wheel;
	}

	wheel = 2;
	for (int i = 0; SmallPrime[i] <= KData.Wheel && SmallPrime[i] < 24;) {
		wheel *= SmallPrime[i++];
	}

	if (wheel > 30 && n > wheel) {
		return wheel;
	}

	const int powten = (int)(log((double)n) / log(10.0) + 0.1);
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
		wheel = 223092870u * 29;
	}

	return wheel;
}

//set sieve buffer size and adjust wheel based
//on cpu L2 cache size and n
/*
static int getSieveCacheSize(const ltype n, int wheel)
{
	const int cachesize = n / wheel;
	return (cachesize >> BSHIFT) + 1;
}*/

//set sieve buffer size and adjust wheel based
//on cpu L2 cache size and n
static int getWheel(const ltype n)
{
	int wheel = getDefaultWheel(n);

	int cachesize = n / wheel;

	int blocks = cachesize / (Config.CpuL2Size << 13);

	wheel *= (blocks + 1);

	return wheel;
}

//get prime number with diff Result in array Prime
//segmented sieve of EratoSieve to enum prime number
static int getPrime(const int n, uchar prime[])
{
	int pi1n = 1;
	if (prime) {
		prime[1] = 2;
	}

	for (int start = 0, sleng = SEGMENT_SIZE ; start < n; start += sleng) {
		utype bitarray[(SEGMENT_SIZE >> (BSHIFT + 1)) + 100];
		if (sleng >= n - start) {
			sleng = n - start;
		}
#if FAST_CROSS
		segmentedEratoSieve1(bitarray, start, sleng + 16, true);
#else
		segmentedEratoSieve2(bitarray, start, sleng + 16, true);
#endif
		if (prime) {
			pi1n += savePrimeDiff(bitarray, start, sleng, prime + pi1n + 1);
		} else {
			pi1n += countZeroBitsArray(bitarray, sleng >> 1);
		}
	}

	assert(pi1n < sizeof(Prime) / sizeof(Prime[0]));
	if (Config.ShowLog) {
		printf("pi(%d) = %d\n", n, pi1n);
	}

	return pi1n;
}

static int getSegKtuplet(const int n)
{
	int pi2n = 0;

#if (OMP)
	omp_set_num_threads(Config.Threads);
	#pragma omp parallel for reduction(+:pi2n) if (n >= MINN * 3)
#endif
	for (int start = 0; start < n; start += SEGMENT_SIZE) {
		utype bitarray[(SEGMENT_SIZE >> (BSHIFT + 1)) + 100];
		int sleng = SEGMENT_SIZE;
		if (sleng >= n - start) {
			sleng = n - start;
		}
#if FAST_CROSS
		segmentedEratoSieve1(bitarray, start, sleng + 16, true);
#else
		segmentedEratoSieve2(bitarray, start, sleng + 16, true);
#endif
		if (PI2N == Config.Gptk)
			pi2n += savePi2Pattern(bitarray, start, sleng, sleng, 0);
		else
			pi2n += savePikPattern(bitarray, start, sleng, sleng, 0);
	}

	return pi2n;
}

//bad performance!!!
static int getSegPartition(const ltype sum, int bitleng)
{
	int gp = 0;
	if (bitleng > sum / 2) {
		bitleng = sum / 2;
	}

	bitleng += bitleng & 1;

#if (OMP)
	omp_set_num_threads(Config.Threads);
	#pragma omp parallel for if (sum >= MINN * 3)
#endif
	for (int start = 0; start < bitleng; start += SEGMENT_SIZE) {
		utype bitarray[(SEGMENT_SIZE >> (BSHIFT + 1)) + 100];
		int sleng = SEGMENT_SIZE;
		if (sleng >= bitleng - start) {
			sleng = bitleng - start;
		}
		segmentedEratoSieve1(bitarray, start, sleng + 16, true);
		reverseBitArray(bitarray, sleng >> 1);
		segmentedEratoSieve1(bitarray, sum - start - sleng, sleng + 16, false);

		gp += countZeroBitsArray(bitarray, (sleng >> 1));
	}

	return gp;
}

//get small parathion or ktuplet prime in range[0 - min(Wheel, sqrt(n))]
//if n is less than a small fix value MIN
static int getSmallGpt(const ltype n)
{
	double ts = getTime( );

	//adjust leng for last few number
	int leng = 0;
	if (n < MINN) {
		if (GP0N == Config.Gptk) {
			leng = n / 2;
		} else if (PI1N == Config.Gptk) {
			leng = n;
		} else if (PI2N >= Config.Gptk) {
			leng = n - (Config.Kgap - 1);
		}
	} else {
		leng = KData.SqrtN;
		if (leng < KData.Wheel) {
			leng = KData.Wheel;
		} else if (PI2N <= Config.Gptk) {
			leng += 1;
		}
	}

	int ret = 0;
	if (GP0N == Config.Gptk) {
		ret = getSegPartition(n, leng);
	} else if (PI1N == Config.Gptk) {
		ret = getPrime(leng + 1, 0);
	} else if (PI2N <= Config.Gptk) {
		ret = getSegKtuplet(leng);
	}

	if (Config.ShowLog) {
		printf("sieve small n = %lld, leng = %d", n, leng);
		printf("\nand small ret = %d, and time use %.lf ms\n", ret, getTime( ) - ts);
	}

	return ret;
}

//init Prime, Pattern, Startp
static int initParam(const ltype n)
{
	KData.N = n;

	KData.SqrtN = (int)sqrt((double)n + 0.1);

#if TEST_BENCHMARK
	if (n > atoint64("e12", 1000000) && Config.Gptk == GP1N)
		n = 1000000000;
#endif

	const int wheel = getWheel(n);

	KData.Wheel = wheel;

	KData.SieveCache = ((n / wheel) >> BSHIFT) + 1;

	const int factorial = getFactorial(wheel);

	KData.Factorial = factorial;

	KData.firstSievedIndex = getFirstPrime(wheel);

	KData.Patterns = initPattern(n, wheel);

	initStartp(wheel);

	LastTask.Gptk = Config.Gptk;
	LastTask.Wheel = KData.Wheel;
	LastTask.N = KData.N;
	LastTask.Patterns = KData.Patterns;

	if (Config.ShowLog) {
		printf("wheel = %d * %d, ",
				factorial, wheel /factorial);
		printf("cachesize = %d k\n", (KData.SieveCache << BSHIFT) >> 13);
	}

	return wheel;
}

static bool excuteCmd(const char*);
//
static int loadConfig(const char* configName)
{
	char linebuf[400];
	char configdata[400] = {0};
	char cmd[400];

	if (!freopen(configName, "rb", stdin)) {
		printf("create a default config file %s\n", configName);
		freopen(configName, "wb", stdout);
		printf(ConfigDataFormat, Config.Threads, Config.CpuL2Size, "E10 012");
		freopen(CONSOLE, "w", stdout);
		freopen(CONSOLE, "r", stdin);
		remove(configName);
		return 0;
	}

	while (gets(linebuf)) {
		strcat(configdata, strcat(linebuf, "\n"));
	}

	if (sscanf(configdata, ConfigDataFormat,
				&Config.Threads, &Config.CpuL2Size, cmd) != 3) {
		printf("invalid config data format: \n%s\n%s\n",
				configdata, ConfigDataFormat);
		puts(configdata);
		freopen(CONSOLE, "r", stdin);
		return 1;
	}
	//read last task data
	if (Config.ShowLog) {
		puts(configdata);
	}

	if (cmd[0]) {
		excuteCmd(cmd);
	}

	freopen(CONSOLE, "r", stdin);

	return 0;
}

static void saveCurrentTask(struct TaskConfig& curtask)
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

	printf(TaskDataFormat, curtask.Gptk,
			curtask.Kgap, curtask.Wheel, curtask.Patterns,
			curtask.Tasks, curtask.Pbegi, curtask.Pendi,
			curtask.N, curtask.Result);

	freopen(CONSOLE, "r", stdin);
	freopen(CONSOLE, "w", stdout);
}

static int readTaskData(struct TaskConfig &curtask)
{
	int ret = 0;
	char linebuf[400] = {0};
	char taskdata[400] = {0};
	freopen("prime.ta", "rb", stdin);

	while (gets(linebuf)) {
		strcat(taskdata, strcat(linebuf, "\n"));
	}

	//read last task data
	struct TaskConfig tmp = {0};
	if (sscanf(taskdata, TaskDataFormat, &tmp.Gptk,
		&tmp.Kgap, &tmp.Wheel, &tmp.Patterns,
		&tmp.Tasks,	&tmp.Pbegi, &tmp.Pendi,
		&tmp.N, &tmp.Result) != 9) {
		printf("invalid task data format: %s : %s\n", taskdata, TaskDataFormat);
		ret = -1;
	}

	//check task data
	if (curtask.N != tmp.N || curtask.Wheel != tmp.Wheel ||
		curtask.Gptk != tmp.Gptk || curtask.Kgap != tmp.Kgap) {
		printf("last task \n%s is not match with current task\n", taskdata);
		ret = -2;
	}

	if (ret == 0) {
		curtask = tmp;
	}

	freopen(CONSOLE, "r", stdin);
	return ret;
}

//
static int loadLastTask(struct TaskConfig &curtask)
{
	int ret = 0;
	if (!freopen("prime.ta", "rb", stdin)) {
		puts("create a default task file prime.ta\n");
		saveCurrentTask(curtask);
	}

	struct TaskConfig tmp = curtask;
	if (readTaskData(tmp) != 0) {
		curtask.Result = curtask.Pbegi = curtask.Pendi = 0;
		saveCurrentTask(curtask);
		tmp = curtask;
	}

	curtask.Result = tmp.Result;
	curtask.Pbegi = tmp.Pbegi;
	if (tmp.Tasks > 0) {
		curtask.Pendi = tmp.Pbegi + tmp.Patterns / tmp.Tasks + 1;
	} else {
		curtask.Pendi = tmp.Pendi;
	}

	curtask.Tasks = tmp.Tasks;

	if (curtask.Pendi > KData.Patterns) {
		curtask.Pendi = KData.Patterns;
	}

	if (Config.ShowLog) {
		printf("load last Task Data with pattern[%d - %d] ok\n",
				curtask.Pbegi, curtask.Pendi);
	}

	freopen(CONSOLE, "r", stdin);
	return ret;
}

//
static ltype getGpKtPn(const ltype n, int pn, bool addsmall)
{
	assert(n < atoint64(MAXN, 100) + 1000);
	assert(Config.Gptk >= 0 && Config.Gptk < 10);

	if (GP0N == Config.Gptk) {
		assert(n % 2 == 0);
	}

	if (n >= MINN) {
		initParam(n);
	}

	ltype sgn = 0, gptn = 0;
	if (addsmall) {
		sgn = getSmallGpt(n);
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
		if (Config.SaveTask && loadLastTask(LastTask) >= 0) {
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
		if (Config.SaveTask && pbegi < pendi) {
			LastTask.Pbegi = pendi;
			LastTask.Result = gptn;
			saveCurrentTask(LastTask);
		}
	}

	KData.Wheel = 0;

	return gptn + sgn;
}

static void printResult(const ltype n, ltype gptn, double ts)
{
	if (KData.UseKtable) {
		putchar('\n');
		ltype sum = 0;
		for (int i = 0; Ktable[i]; i++) {
			sum += Ktable[i];
			printf("%4d%s %lld\n", i + 1, TABLE_GAP, sum);
		}
	}

	if (Config.Gptk < PIKN) {
		printf(PrintFormat[Config.Gptk], n, gptn);
	} else {
		printf(PrintFormat[PIKN], KData.Kpattern[0], n, gptn);
	}

	if (Config.ShowTime) {
		double timeuse = getTime( ) - ts;
		if (timeuse > 10000) {
			printf(", time use %.3lf sec", timeuse / 1000);
		} else {
			printf(", time use %.2lf ms", timeuse);
		}
	}

	//		if (Config.ShowLog) {
	//			printf(", ave = %d", (int)(gptn / KData.Patterns));
	//		}
	putchar('\n');
}

//get goldbach partition, or twin primes, or ktuplet
static ltype Ktprime(const ltype n, int pn)
{
	double ts = getTime( );

	ltype gptn = 0;
	bool addsmall = true;
	if (KData.Patterns == 0) {
		LastTask.Gptk = -1;
	}

	if (n > KData.N + 1000) {
		const int sqrtn = (int)sqrt((double)n) + 1001;
		getPrime(sqrtn, Prime);
	}

	//optimize for Ktuplet prime with small n > e7
	if (Config.Gptk != GP0N && n > atoint64("e13", 0)) {
		int sqrtn = (int)sqrt((double)n + 0.1);
		int wheel = getWheel(n);
		int maxn = sqrtn > wheel ? sqrtn : wheel;
		wheel = KData.Wheel;
		KData.Wheel = 0;
		gptn = getGpKtPn(maxn + Config.Kgap - 1, 0, addsmall);
		//restore the wheel
		KData.Wheel = wheel;
		if (Config.ShowLog) {
			//			printf("PI%d[%d] = %lld, and time use %.lf ms\n",
			//					Config.Gptk, maxn, gptn, getTime( ) - ts);
		}
		addsmall = !addsmall;
	}

	if (Config.Gptk == PIKN && pn == 0 && n > atoint64(TABLE_GAP, 0)) {
		memset(Ktable, 0, sizeof(Ktable));
		KData.UseKtable = true;
	} else {
		KData.UseKtable = false;
	}

	gptn += getGpKtPn(n, pn, addsmall);

	if (Config.ShowResult) {
		printResult(n, gptn, ts);
	}

	return gptn;
}

//test case number
static int startTest(int tesecase, bool rwflag)
{
	srand((uint)time(NULL));
	double ts = getTime( );
	const char* dataFile[ ] = {"prime.gp", "prime.pi", "prime.pi2", "prime.pik"};
	printf("-----------start %s test -----------\n", dataFile[Config.Gptk]);

	if (rwflag) {
		if (!(freopen(dataFile[Config.Gptk], "rb", stdin))) {
			printf("can not read test data file %s\n", dataFile[Config.Gptk]);
			freopen(CONSOLE, "r", stdin);
			Config.ShowResult = true;
			return -1;
		}
		Config.ShowResult = Config.PrintGap = Config.ShowLog = false;
	} else {
		if (!(freopen(dataFile[Config.Gptk], "wb", stdout))) {
			puts("can not write test data file");
			freopen(CONSOLE, "r", stdin);
			Config.ShowResult = true;
			return -2;
		}
		Config.PrintGap = Config.ShowTime = false;
	}

	int fails = 0;
	for (int i = 1; i <= tesecase; i++) {
		if (!rwflag) {
			ltype n = rand( ) * rand( );
			n = (n + 4) * 2 + 0;
			if (n < 1000000) {
				n = 4 * n + 1000000;
			}
			Ktprime(n, 0);
		} else {
			ltype res, n;
			char linebuf[256] = {0};
			gets(linebuf);
			if (sscanf(linebuf, PrintFormat[Config.Gptk], &n, &res) != 2) {
				printf("line %d is wrong data\n", i);
				if (fails++ > 30) {
					break;
				} else {
					continue;
				}
			}

			ltype gptn = Ktprime(n, 0);
			if (gptn != res) {
				printf("case %d with wrong result %lld, ", i, gptn);
				printf(PrintFormat[Config.Gptk], n, res);
				putchar('\n');
			}
		}
		if ((i & 63) == 0) {
			printf("case pass %d%%\r", i * 100 / tesecase);
		}
	}

	printf("test case time use %.lf ms\n", getTime( ) - ts);

	Config.ShowResult = true;
	freopen(CONSOLE, "w", stdout);
	freopen(CONSOLE, "r", stdin);

	return 0;
}

//list Gptk by the input Result start, end, step
static void listDiffGpt(const char cmdparams[][80], int cmdi)
{
	double ts = getTime( );

	int ni = 1, step = 2;
	ltype start = 1000000000, end = start + 1000;
	ltype buf[ ] = {0, start, end, step, 0};

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

	printf("calculating %s\n", KtupletName[Config.Gptk]);

	if (Config.SaveToFile) {
		Config.PrintGap = Config.ShowTime = false;
		freopen("batch.txt", "wb", stdout);
	}

	printf("%lld:%d:%d\n\n", start, (int)(2 + end - start) / 2, step);

	int pcnt = 0;
	ltype allSum = 0;
	for (ltype n = start; n <= end; n += step) {
		pcnt++;
		if (isdigit(cmdparams[4][0])) {
			printf("%d ", pcnt);
		}
		allSum += Ktprime(n, 0);
	}

	if (GP0N == Config.Gptk) {
		printf("average = %lld, ", allSum / pcnt);
	}

	printf("all case time use %.lf ms\n", getTime( ) - ts);
	freopen(CONSOLE, "w", stdout);
}

static void listPowGpt(const char cmdparams[][80], int cmdi)
{
	printf("calculating %s ", KtupletName[Config.Gptk]);

	int m = atoint64(cmdparams[cmdi + 1], 10);
	int startindex = atoint64(cmdparams[cmdi + 2], 5);
	int endindex = atoint64(cmdparams[cmdi + 3], 10);

	printf("in %d^%d - %d^%d\n", m, startindex,
			m, endindex);

	if (m < 2 && m > 10000){
		m = 10;
	}
	if (startindex > endindex) {
		startindex ^= (endindex ^= (startindex ^= endindex));
	}

	if (Config.SaveToFile) {
		Config.PrintGap = Config.ShowTime = false;
		freopen("batch.txt", "wb", stdout);
	}

	Config.ShowResult = false;
	for (ltype i = startindex; i <= endindex; i++) {
		ltype n = (ltype)(pow((double)m, (int)i) + 0.01);
		ltype r = Ktprime(n, 0);
		printf(PrintFormat[Config.Gptk], m, i, r);
		putchar('\n');
	}
	Config.ShowResult = true;
}

//benchMark
static void benchMark(const char cmdparams[][80])
{
	Config.ShowResult = true;
	ltype start = atoint64("e10", 1000000);
	ltype n = atoint64("e10", 0);
	char gptk[20] = {'0', '1', '2', '3'};

	for (int i = 0, ni = 1; cmdparams[i][0]; i++) {
		char c = cmdparams[i][0];
		if (isdigit(c) || toupper(c) == 'E') {
			ltype tmp = atoint64(cmdparams[i], 0);
			if (ni++ == 1) {
				start = tmp;
			} else if (ni == 3) {
				n = tmp;
			} else if (tmp < 222) {
				strcpy(gptk, cmdparams[i]);
			}
		}
	}

	if (Config.SaveToFile) {
		freopen("benchmark.txt", "wb", stdout);
	}

	ltype gap = start;
	for (; start * 10 <= n; ) {
		for (int j = 0; j < 9; j++) {
			for (int k = 0; gptk[k]; k++) {
				Config.Gptk = gptk[k] - '0';
				Ktprime((j + 1) * gap, 0);
			}
			puts("");
		}
		start *= 10;
		gap *= 10;
	}

	freopen(CONSOLE, "w", stdout);
}

//auto set best threads
static void benchMarkByThreads(const ltype n, int threads)
{
	int mints = 1234567890;
	int bestthreads = 4;
	Config.ShowResult = true;
	for (int i = 1; i <= threads; i++) {
		Config.Threads = i;
		double ts = getTime();
		Ktprime(n, 1000);
		ts = getTime() - ts;
		if (mints > ts) {
			bestthreads = Config.Threads;
			mints = ts;
		}
	}
	Config.Threads = bestthreads;
	printf("best thread = %d\n", bestthreads);
}

//test pi, pi2, gp function
static void testGpt( )
{
	const char* const gtpdata[][6] =
	{
		{"e07", "210",   "000", "00038807",  "00664579",  "00058980"},
		{"e14", "19",    "080", "14385672",  "78519928",  "14446361"},
		{"e08", "210",   "000", "00291400",  "05761455",  "00440312"},
		{"e09", "2310",  "000", "02274205",  "50847534",  "03424506"},
		{"e10", "30030", "000", "18200488",  "455052511", "27412679"},
		{"e11", "30030", "100", "15059883",  "71518765",  "15114393"},
		{"e12", "17",    "200", "16756829",  "81689179",  "16801492"},
		{"e13", "17",    "060", "14200452",  "75326388",  "14236925"},
		{"e15", "19",    "060", "14693072",  "90213054",  "14892236"}
	};

	Config.CpuL2Size = 1 << 10;
	Config.PrintGap = false;
	Config.ShowTime = true;

	for (int i = 0; i < sizeof(gtpdata) / sizeof(gtpdata[0]); i++) {
		ltype n = atoint64(gtpdata[i][0], 0);
		int psize = atoi(gtpdata[i][2]);
		for (int j = 0; j < 3; j ++) {
			KData.Wheel = atoi(gtpdata[i][0 + 1]);
			Config.Gptk = j;
			if (atoi(gtpdata[i][j + 3]) != Ktprime(n, psize)) {
				printf("%s : %s fail !!!\n", gtpdata[i][0], gtpdata[i][j + 3]);
			}
		}
		putchar('\n');
	}
}

//test pi, pi2, gp function
static void testPik( )
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
	Config.PrintGap = false;
	Config.ShowTime = true;

	for (int i = 0; i < sizeof(pidata) / sizeof(pidata[0]); i++) {
		excuteCmd(pidata[i][0]);

		for (int j = 1; j < sizeof(pidata[i])/sizeof(pidata[i][0]); j ++) {
			ltype n = atoint64(pidata[i][j], 100);
			ltype r = atoint64(pidata[i][j] + 3, 0);
			if (r != Ktprime(n, 0)) {
				printf("%s != %lld fail\n", pidata[i][j], r);
			}
		}
		putchar('\n');
	}
}

//test pi2 function
static void testPi2( )
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

	Config.PrintGap = false;
	Config.ShowTime = false;
	Config.Gptk = PI2N;

	for (int i = 0; i < sizeof(pi2data) / sizeof(pi2data[0]); i++) {
		const char* pdata = pi2data[i];
		ltype n = atoint64(pdata, 0);
		pdata += 6;
		if (n < atoint64("e3", 0)) {
			continue;
		}
		for (int j = 1; j < 7; j ++) {
			Config.Kgap = 2 * j;
			while (isspace(*pdata))
				pdata++;

			ltype ret = atoint64(pdata, 0);
			ltype cal = Ktprime(n, 0);
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
//	Config.CpuL2Size = cpuinfo[2] >> 16;

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
#else
	Config.Threads = sysconf(_SC_NPROCESSORS_CONF);
#endif

	return Config.Threads;
}

//print the Gptk info
static void printInfo( )
{
	puts("---------------------------------------------------------------");
	puts("---------------------------------------------------------------");

	printf("\
	1.%s G(n)\n \
	2.%s PI(n)\n \
	3.Twin/Cousin/Sexy/p, p+2n/ prime pairs PI2(n)\n \
	4.Ktuplet prime PIk(n)\n (n < %s) version %s\n",
	KtupletName[0], KtupletName[1], MAXN, KVERSION);

	puts("Copyright (c) by Huang Yuanbing 2008 - 2013 bailuzhou@163.com");

#ifdef _MSC_VER
	printf("Compiled by MS/vc++ %d", _MSC_VER);
#else
	printf("Compiled by Mingw g++ %d.%d.%d",
			__GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#endif

#if _M_AMD64 || __x86_64__
	printf(" on x64");
#endif
	printf(" on %s %s\n", __TIME__, __DATE__);

	getSystemInfo();
	getCpuInfo();

	printf("Work threads = %d, POPCNT = %d, ASM_X86 = %d\n", \
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
		KData.Kpattern[0] = ns;
		KData.Kpattern[ns] = 0;
		return true;
	}

	for (int j = 2, k = 2; pstr[j]; j++) {
		char dig = pstr[j] - '0';
		if (dig <= 0 || dig > 9) {
			continue;
		}
		if (dig % 2 == 0 && j < 8) {
			KData.Kpattern[k - 1] = dig;
		} else {
			KData.Kpattern[k - 1] = dig * 10 + pstr[++j] - '0';
		}
		assert(KData.Kpattern[k - 1] % 2 == 0);
		if (Config.ShowLog) {
			printf("%d\n", KData.Kpattern[k - 1] );
		}
		KData.Kpattern[0] = k++;
	}

	return true;
}

static void printKpattern( )
{
	Config.Kgap = KData.Kpattern[1];
	Config.Gptk = KData.Kpattern[0];

	int ktuplets = KData.Kpattern[0];
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
#if _MSC_VER
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
static int parseConfig(char cmdparams[][80])
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
				Config.SaveTask = !Config.SaveTask;
				break;
			case 'D':
				Config.ShowLog = !Config.ShowLog;
				break;
			case 'P':
				Config.ShowTime = !Config.ShowTime;
				break;
			case 'S':
				Config.SaveToFile = !Config.SaveToFile;
				break;
			case 'R':
				Config.CheckPattern = !Config.CheckPattern;
				break;
			case 'C':
				if (tmp <= MAX_L1SIZE && tmp > 15) {
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
			default:
				cmdi = i;
				break;
		}
	}

	return cmdi;
}

//split cmd to cmdparams[] by white space and ';'
static int splitCmdParams(const char* ccmd, char cmdparams[][80])
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

		if (splitCmdParams(cmd, cmdparams) <= 0) {
			return false;
		}

		int cmdi = parseConfig(cmdparams);

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
				testPi2( );
				testGpt( );
				testPik( );
			}
			for (int i = 0; cmdparams[cmdi + 2][i]; i++) {
				Config.Gptk = cmdparams[cmdi + 2][i] - '0';
				bool rwflag = true;
				if (Config.Gptk < 4 && Config.Gptk >= 0) {
					startTest(atoint64(cmdparams[cmdi + 1]), rwflag);
				}
			}
		} else if (cmdc == 'O') {
			puts("------start optimize threads -----------");
			ltype n = atoint64(cmdparams[cmdi + 1], 2000000000);
			int count = atoint64(cmdparams[cmdi + 2], MAX_THREADS / 4);
			benchMarkByThreads(n, count);
		} else if (cmdc == 'L') {
			puts("-----start list G(n)/PI2(n)/PI(n) ------");
			listDiffGpt(cmdparams, cmdi);
		} else if (cmdc == 'I') {
			puts("-------list pow G(n)/PI2(n)/PI(n) ------");
			listPowGpt(cmdparams, cmdi);
		} else if (cmdc == 'N') {
			puts("---------start list patterns -----------");
			ltype n = atoint64(cmdparams[cmdi + 1], 1000000000);
			int count = atoint64(cmdparams[cmdi + 2], 100);
			for (int i = 0; i < 2 * count; i += 2) {
				initParam(n + i);
				printf("patterns %lld = %d\n", n + i, KData.Patterns);
				if (KData.Patterns > sizeof(Pattern)) {
					puts("error : patterns size is small");
				}
			}
			KData.Wheel = 0;
		} else if (cmdc == 'E' || isdigit(cmdc)) {
			puts("---------start G(n)/PI2(n)/PI(n) -------");
			ltype n = atoint64(cmdparams[cmdi], 1000000000);
			bool retry = true;
			for (int i = 0; cmdparams[cmdi + 1][i]; i++) {
				Config.Gptk = cmdparams[cmdi + 1][i] - '0';
				if (Config.Gptk < 4 && Config.Gptk >= 0) {
					int pattern = atoint64(cmdparams[cmdi + 2], 0);
					Ktprime(n, pattern);
					retry = false;
				} else {
					Config.Gptk = GP0N;
				}
			}
			if (retry) {
				Ktprime(n, 0);
			}
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
static void initCache( )
{
	double ts = getTime( );

	simpleEratoSieve(10000);

	initBitTable( );

	KData.N = atoint64("e12", 0);
	const int sqrtn = (int)sqrt((double)KData.N) + 1001;
	int primes = getPrime(sqrtn, Prime);

	if (Config.ShowLog) {
		printf("init Pi(%d) = %d, table use %.2lf ms\n",
				sqrtn, primes, getTime( ) - ts);
	}
}

int main(int argc, char* argv[])
{
	if (argc < 2) {
		printInfo( );
	}

	initCache( );

	for (int i = 1; i < argc; i++) {
		if (argv[i][0] == 'c')
			loadConfig("Kconfig.ini");
		else if (argv[i][0] == 'm')
			doCompile();
		else
			excuteCmd(argv[i]);
	}

	excuteCmd("t4 1e10 01");
//	excuteCmd("d 1e15 0 t4");

	char ccmd[256] = {0};
	while (true) {
		printf("\n[input command] : ");
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
G(1e14) = 90350630388, gcc 4.6.1 0510 -march=native -O2 time use 6481.020 s AMD X4 820
G(1e15) = 783538341852, time use 44.17h -->32h -->24

  n     G(10^n)      PI2(10^n)
  7     38807        664579
  8     291400       5761455
  9     2274205      3424506
  10    18200488     27412679
  11    149091160    224376048
  12    1243722370   1870585220
  13    10533150855  15834664872
  14    90350630388  135780321665
  15    783538341852
  16

feature:
  1   multi thread support x
  2   detail comment x
  3   clear marco/function name x
  4   reduce memory use for initPattern x
  5   split task and save to file
  6   formula calculate
  7   Bx
  8   List G x
  9   config.ini x
  10  print ktuplet table x
  11  optimize threads and cache size x

  12. calculate range [e16, e16+ e9] !!!!
  13. calculate Goldbach partitions split from the code
  14. win32 gui
  15. remove unused marco

Linux g++:
  g++ -Wall -msse4 -O3 -march=native -s -pipe -ffunction-sections -fomit-frame-pointer -lpthread Ktprime.cpp
Mingw/g++:
  g++ -Wall -mpopcnt -mtune=native -O2 -s -pipe -fomit-frame-pointer Ktprime.cpp
MS vc++:
  cl /O2 /Os Ktprime.cpp

command:
	D t2 C2000 M5 e15

c1600 m5 d t4 e15 0 1000
need 28h amd phoenm x4 830
 ****************************************************/

