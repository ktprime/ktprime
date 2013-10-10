/************************************************************
Calculate number of tuplet T(n): (p1, p2, p3) with n = p1 + p2 + p3
(p1 <= p2 <= p3), p1, p2 and p3 are all odd prime number;
copyright (C) 2008 - 2013 by Huang Yuanbing
version 7.0
mail to: bailuzhou@163.com
free use for non-commercial purposes

T(10000001) = 6154337544

2.00G AMD  Althon	 3600+  1555 ms
2.80G AMD  phonem X4 820    120 ms  SSE4A
2.80G AMD  Althon X4 641    140 ms  SSE4A

2.00G Intel croe 2   e2180	937 ms
1.66G Intel core 2   t5500	1045 ms
2.27G Intel core i3  350M	200 ms  SSE4.2 x64

Algorithm:
	for an odd number n = p1 + p2 + p3, p1 <= p2 <= p3, p1, p2, p3 all are primes,
define T(n) as number of 3-tuples (p1, p2, p3),
and G(n) is defined as number of prime Pairs (p1, p2) with n = p1 + p2, p1 <= p2
for fast calculating T(n):
					 	    n
	T(n) = sigma(p1 + p2) sigma
						  p3 = n / 3

	if the max value of prime p3 >= n / 2
		T(n) = sigma G(n - p3)(p3 is in [n / 2, n])
	else if p3 > 3N/7
		T(n) = sigma G(n - p3)(p3 is in [3N / 7, n / 2]) with p1 < p3 or p2 < p3
	else if p3 >= n / 3
		T(n) = sigma G(n - p3)(p3 is in [n / 3, 3N / 7]) with p1 < p3 or p2 < p3

http://tieba.baidu.com/f?kz=830665538
**************************************************************/

# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <math.h>
# include <ctype.h>
# include <time.h>
# include <assert.h>

# define MAXN           10020000 // 587732175
# define PRIME_DIFF     1
# define TABLE_GP       1
# define LOAD_SEG       0
# define OMP            0
# define BEST_RATION    2 / 7
# define GSINDEX        10
# define MAX_THREADS    31

//use of the SSE4.2/ SSE4a  POPCNT instruction for fast bit counting.
#if _MSC_VER > 1400
	# define POPCNT      1
	# include <intrin.h>
#elif  (__GNUC__ * 10 + __GNUC_MINOR__ > 44)
	# define POPCNT      1
	//# include <x86intrin.h>
	# include <popcntintrin.h>
#else
	# define POPCNT      0
#endif

# if OMP
	#include <omp.h>
# endif

typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;

#ifdef _MSC_VER
	typedef __int64 int64;
	typedef unsigned __int64 uint64;
	# pragma warning(disable: 4996 4244 4127 4505 4018 6328 6031)
#else
	typedef long long int64;
	typedef unsigned long long uint64;
#endif

#if MAXN > 210000000
	# define FACTP      2310
	# define PATTERNS   480
#else
	# define FACTP      210
	# define PATTERNS   48
#endif

# define BIT 1
# if BIT == 1
	typedef uchar utype;
	# define MOVE 3
# elif BIT == 2
	typedef ushort utype;
	# define MOVE 4
# elif BIT == 3
	typedef uint utype;
	# define MOVE 5
# endif

static const uint MASK = (1 << MOVE) - 1;

# define TEST_BIT(a, n)    (a[(n) >> MOVE] & ((utype)1 << ((n) & MASK)))
# define TEST_BIT2(a, n)   TEST_BIT(a, n / 2)
# define SET_BIT1(a, n)    a[(n) >> MOVE] |= (utype)1 << ((n) & MASK)
# define SET_BIT0(a, n)    a[(n) >> MOVE] &= ~((utype)1 << ((n) & MASK))
# define MIN_DIF_BITS      (2 * 64)

static const char* HelpInfo = "\
	[B: Benchmark]\n\
	[H: Help ccmd]\n\
	[U: Unit test from prime.tn]\n\
	[T: Threads number (n 2 - 64)]\n\
	[P: Print time]\n\
	[A: ASC_DIREC flag]\n\
	[D: Debug log]\n\
	[C: Clear cache]\n\
	[Z: Exit programming]\n\
	[S: Save result to file]\n\
	[L: List Tn (beg) (end/count) (step)]\n\
	[G: List gp (beg) (end/count) (step)]";

static const char* HelpCommand = "\n\
	All command/param as follow:\n\
	B, B 10000\n\
	C\n\
	T2-32\n\
	H, A, S, P\n\
	U 1000, U 1000 10*10\n\
	N e8, N 2e8+20\n\
	R 120000*1000 100\n\
	G 2^31 1000, G 2e10*3 3e3\n\
	G 10^11*2 100, G 400000000+100 2000\n\
	L 2e9-100 1000 10]\n\
	L e8-100 e8 e3]";

//config
static int MapIndex[FACTP];
//pattern
static int Pattern[PATTERNS];
//pair
static int Pairs[FACTP / 2][PATTERNS / 2 + 4];

//pair mask
//static uint PairMask[FACTP / 2][(PATTERNS / 32) + 8];

//the first 7 prime numbers
static const uchar SmallPrime[ ] =
{
	3, 5, 7, 11, 13, 17, 19
};

//config
static struct
{
	//print result
	bool ShowRet;
	//create g(n) table
	bool CreateGpt;
	//print calculating time
	bool ShowTime;

	//save result to log file
	bool FileSave;
	//show Debug Log
	bool ShowLog;

	//set sieve direction
	bool EnumMax;
	//number of threads
	int Threads;

	int TestCase;
}
Config =
{
	true, false, true,
	false, false, true, 2
};

const int MOVEARRAYS = ((MAXN / FACTP * BEST_RATION) >> (MOVE - 1)) + 100;

#if TABLE_GP
#if MAXN < 10040000
	#if LOAD_SEG == 0
		static ushort GPT[MAXN * BEST_RATION + 1000];
	#else
		static ushort GPT[MAXN * 10 + 1000];
	#endif
#else
	#if LOAD_SEG == 0
		static int GPT[MAXN * BEST_RATION + 1000];
	#else
		static int GPT[MAXN / 16 + 1000];
	#endif
#endif
#endif

/** src bit prime array **/
static utype SrcBitArray[MAXN >> MOVE];
/** bit reverse of SrcBitArray */
static utype RevBitArray[MAXN >> MOVE];

/** left shift move array of SrcBitArray **/
#if PRIME_DIFF != 2
	static utype LeftShiftArray[PATTERNS][MASK + 1][MOVEARRAYS];
#else
	static utype LeftShiftArray[PATTERNS][MOVEARRAYS];
#endif

/** right shift move array of RevBitArray **/
#if PRIME_DIFF == 0
	static utype RightShiftArray[PATTERNS][MASK + 1][MOVEARRAYS];
#else
	static utype RightShiftArray[PATTERNS][MOVEARRAYS];
#endif

//WordReverse[i] is equal to the bit reverse of i (i < 2^16)
static ushort WordReverse[1 << 16];
//WordNumBit0[i] is bits number of binary representation of i
static uchar WordNumBit0[1 << 16];

static int getGpFromTable(int);
static int64 coreSieve(int, int, int);
static int createGpTable(const int n, int offset, int step);

static struct ThreadInfo
{
	int n;
	int	offset;
	int step;
	int64 tn;
} Tparam[MAX_THREADS];

#ifdef _WIN32
	# include <windows.h>
	# define CONSOLE "CON"
static DWORD WINAPI threadProcTn(void* ptinfo)
#else
	# include <pthread.h>
	# include <sys/time.h>
	# define CONSOLE "/dev/tty"
static void* threadProcTn(void* ptinfo)
#endif
{
	ThreadInfo* pThreadInfo = (ThreadInfo*)(ptinfo);
	pThreadInfo->tn = coreSieve(pThreadInfo->n, pThreadInfo->offset,
			pThreadInfo->step);

	return 0;
}

#ifdef _WIN32
static DWORD WINAPI threadProcGp(void* ptinfo)
#else
static void* threadProcGp(void* ptinfo)
#endif
{
	ThreadInfo* pThreadInfo = (ThreadInfo*)(ptinfo);
#if TABLE_GP
	createGpTable(pThreadInfo->n, pThreadInfo->offset / 2,
			pThreadInfo->step / 2);
#endif

	return 0;
}

static int64 startWorkThread(int threads, const int n, int* proc = 0)
{
	int64 tn = 0;
	int i;

	if (threads > MAX_THREADS) {
		threads = 8;
	}

	for (i = 0; i < threads; i++) {
		Tparam[i].n = n;
		Tparam[i].step = threads * 2;
		Tparam[i].offset = i * 2 + 1;
	}

#ifdef _WIN32
	HANDLE thandle[MAX_THREADS];
	DWORD tid[MAX_THREADS];
	for (i = 0; i < threads; i++) {
		thandle[i] = CreateThread(NULL, 0,
				(LPTHREAD_START_ROUTINE)proc,
				(LPVOID)(&Tparam[i]), 0, &tid[i]);
		Sleep(5);
		if (thandle[i] == NULL) {
			printf("create win32 thread error %ld\n", GetLastError());
		}
	}
	WaitForMultipleObjects(threads, thandle, true, INFINITE);
	for (i = 0; i < threads; i++) {
		if (thandle[i]) {
			CloseHandle(thandle[i]);
		}
	}
#else
	pthread_t tid[MAX_THREADS];
	for (i = 0; i < threads; i++) {
		int error = pthread_create(&tid[i], NULL, void* (void*)(proc), &Tparam[i]);
		if (error != 0) {
			printf("create posix thread error %d\n", error);
		}
	}
	for (i = 0; i < threads; i++) {
		pthread_join(tid[i], NULL);
	}
#endif

	for (i = 0; i < threads; i++)
		tn += Tparam[i].tn;

	return tn;
}


#ifdef _WIN32
static LONG WINAPI filterFunc(DWORD dwExceptionCode)
{
	return((dwExceptionCode == STATUS_ILLEGAL_INSTRUCTION)
		? EXCEPTION_EXECUTE_HANDLER : EXCEPTION_CONTINUE_SEARCH);
}
#endif

static int getCpuInfo(  )
{
#ifdef _WIN32
	SYSTEM_INFO si;
	GetSystemInfo(&si);

	if (si.dwNumberOfProcessors > Config.Threads) {
		Config.Threads = si.dwNumberOfProcessors;
	}

	int popcnt = 0;
#if POPCNT && _MSC_VER
	__try {
		popcnt = _mm_popcnt_u32(1);
	}
	__except (filterFunc(GetExceptionCode())) {
		popcnt = 0;
	}
#elif POPCNT
	popcnt = _mm_popcnt_u32(1);
#endif

	printf("Cpu cores = %ld, SSE4 with Popcnt = %d\n",
		si.dwNumberOfProcessors, popcnt);

	return si.dwNumberOfProcessors;
#endif
}

//get current time
static double getTime( )
{
#ifdef _WIN32
	static LARGE_INTEGER s_freq;
	LARGE_INTEGER performanceCount;
	if (s_freq.QuadPart == 0 && !QueryPerformanceFrequency(&s_freq))
		return -1;
	QueryPerformanceCounter(&performanceCount);
	return 1000.* performanceCount.QuadPart / (double)s_freq.QuadPart;
#else
	struct timeval tmVal;
	//struct timezone tmZone;
	gettimeofday(&tmVal, NULL);
	return tmVal.tv_sec * 1000. + tmVal.tv_usec / 1000.;
#endif
}

/**
  the sieve of Eratosthenes implementation by bit packing
all prime less than 2^16 will be saved in prime buffer List
Prime[0] is the first sieve prime, Prime[i] is the difference
of the adjacent prime, Prime[i] = Prime[i] - Prime[i - 1];
*/
static int eratosSieve(const int maxn)
{
	int primes = 1;
	SET_BIT1(SrcBitArray, 0);

	for (int p = 3; p < maxn; p += 2) {
		if (!TEST_BIT2(SrcBitArray, p)) {
			primes++;
			if (p > maxn / p)
				continue;
			for (int j = p * p / 2; 2 * j < maxn; j += p)
				SET_BIT1(SrcBitArray, j);
		}
	}

	if (Config.ShowLog)
		printf("sieve : PI[%d] = %d\n", primes, maxn);

	return primes;
}

//reverse bit of a byte with binary representation
static uchar reverseByte(const uchar c)
{
	uchar n = (c & 0x55) << 1 | (c & 0xAA) >> 1;
	n = (n & 0x33) << 2 | (n & 0xCC) >> 2;
	n = (n & 0x0F) << 4 | (n & 0xF0) >> 4;
	return n;
}

//reverse utype data by bit binary representation
static utype reverseUtype(const utype n)
{
	utype ret = 0;
	int sizeutype = sizeof(n);
	if (sizeutype == sizeof(uint))
		ret = WordReverse[n >> 16] | (WordReverse[(ushort)n] << 16);
	else if (sizeutype == sizeof(ushort))
		ret = WordReverse[n];
	else if (sizeutype == sizeof(uchar))
		ret = WordReverse[n] >> 8;
	else if (sizeutype == sizeof(uint64)) {
		uint low = n >> 32, hig = n;
		low = WordReverse[low >> 16] | (WordReverse[(ushort)low] << 16);
		hig = WordReverse[hig >> 16] | (WordReverse[(ushort)hig] << 16);
		ret = (((utype)hig) << 32) | low;
	}
	return ret;
}

//reverse and adjust LeftShiftArray and RightShiftArray
//the most complicated code here
static void adjustLRevBitArray(const int n)
{
	memset(LeftShiftArray, (uint)(-1), sizeof(LeftShiftArray));

	int st = 11;
#if (FACTP % 7 != 0)
		st = 7;
#elif (FACTP % 11 == 0)
		st = 13;
#endif

	for (int p = st; p < (n * BEST_RATION) * 2; p += 2) {
		if (0 == TEST_BIT2(SrcBitArray, p)) {
			int mapi = MapIndex[p % FACTP];
#if PRIME_DIFF != 2
			SET_BIT0(LeftShiftArray[mapi][0], p / FACTP);
			assert(p / FACTP < sizeof(LeftShiftArray[0][0]) * (MASK + 1));
#else
			SET_BIT0(LeftShiftArray[mapi], p / FACTP);
			assert(p / FACTP < sizeof(LeftShiftArray[0]) * (MASK + 1));
#endif
		}
	}

#if PRIME_DIFF != 0
	const int rlen = sizeof(RightShiftArray[0]) / sizeof(utype) - 1;
#else
	const int rlen = sizeof(RightShiftArray[0][0]) / sizeof(utype) - 1;
#endif

	//the core move part code
	for (int i = 0; i < PATTERNS; i++) {

#if PRIME_DIFF != 2
		for (int k = 0; k < MASK + 1; k++)
#endif
		for (int j = 0; j <= rlen; j++) {

#if PRIME_DIFF != 2
			utype movevalue = (LeftShiftArray[i][0][j] >> k) |
				(LeftShiftArray[i][0][j + 1] << (MASK + 1 - k));
			if (k == 0)
				movevalue = LeftShiftArray[i][0][j];
			LeftShiftArray[i][k][j] = movevalue;
#else
			utype movevalue = LeftShiftArray[i][j];
			RightShiftArray[i][rlen - j] = reverseUtype(movevalue);
#endif

#if PRIME_DIFF == 0
			if (k == 0)
				RightShiftArray[i][k][rlen - j] = reverseUtype(movevalue);
			else
				RightShiftArray[i][k][j] = (RightShiftArray[i][0][j] >> k) |
					(RightShiftArray[i][0][j + 1] << (MASK + 1 - k));
#elif PRIME_DIFF == 1
			if (k == 0)
				RightShiftArray[i][rlen - j] = reverseUtype(movevalue);
#endif
		}
	}
}

//get g(n) with the small number in SmallPrime list
static int getSmallGp(const int n)
{
	int gp = 0;
	for (int j = 0; FACTP % SmallPrime[j] == 0; j++) {
		int nextPrime = n - SmallPrime[j];
//		if (nextPrime < 0)
//			break;
		if (!TEST_BIT2(SrcBitArray, SmallPrime[j]) &&
			!TEST_BIT2(SrcBitArray, nextPrime))
			gp ++;
	}
	return gp;
}

//get g(n) with beg_bit + LIM < end1_bit
static int getBitsGp(int beg_bit, int end_bit)
{
	int gp = 0;
	while (beg_bit <= end_bit) {
		if (!TEST_BIT(SrcBitArray, end_bit) &&
			!TEST_BIT(SrcBitArray, beg_bit))
			gp++;
		beg_bit++, end_bit--;
	}
	return gp;
}

static inline int countBitOnes(const uint n)
{
#if POPCNT
	return 32 - _mm_popcnt_u32(n);
#else
	return 32 - WordNumBit0[(ushort)n] - WordNumBit0[n >> 16];
#endif
}

static inline int countZeroBits(uint64 pbeg[], uint64 pend[], const int bitleng)
{
	int loops = bitleng >> 6;
	int bits0 = 0;

	while (loops-- > 0) {
		uint64 n = *pbeg++ | *pend++;
		//60% cpu time
	#if POPCNT
		//popcnt instruction in intel i7/SSE4.2 and AMD phonem/SSE4A
		#if _M_AMD64 || __x86_64__
		bits0 += _mm_popcnt_u64(n);
		#else
		bits0 += _mm_popcnt_u32(n) + _mm_popcnt_u32(n >> 32);
		#endif
	#else
		uint hig = n >> 32, low = n;
		bits0 += WordNumBit0[(ushort)hig] + WordNumBit0[hig >> 16] +
			  WordNumBit0[(ushort)low] + WordNumBit0[low >> 16];
	#endif
	}

	bits0 = ((bitleng >> 6) << 6) - bits0;

	return bits0;
}

//the core part of this algorithm
//count gp from two arrays, the first array with pattern m1,
// bit position start form beg1_bit to end1_bit
//the second array with part m2 and bit position from
//end2_bit with reverse direction
static int countPairSum(int beg1_bit, int end1_bit, int end2_bit, int m1, int m2)
#if PRIME_DIFF == 2
{
	uint gp = 0, mask;

	end2_bit = sizeof(RightShiftArray[0]) * 8 - end2_bit - 1;

	while ((beg1_bit & 31) && (beg1_bit <= end1_bit)) {
		if (!TEST_BIT(LeftShiftArray[m1], beg1_bit) &&
			!TEST_BIT(RightShiftArray[m2], end2_bit))
			gp++;
		beg1_bit++; end2_bit++;
	}

	if (beg1_bit > end1_bit)
		return gp;

#if 1
	uint* pbeg = (uint*)LeftShiftArray[m1] + (beg1_bit >> 5);
	uint* pend = (uint*)RightShiftArray[m2] + (end2_bit >> 5);

	int shiftBits = end2_bit & 31;
	int loops = (end1_bit - beg1_bit) >> 5;

	//asm to hot spot
	while (loops-- > 0) {
		mask = *pbeg++ | (uint)(*((uint64*)(pend++)) >> shiftBits);
		gp += countBitOnes(mask);
	}

	uint us = *pbeg | (uint)(*((uint64*)pend) >> shiftBits);
	mask = ~((1u << ((1 + end1_bit) & 31)) - 1);
	if (mask != (uint) ~0)
		mask |= us;
	else
		mask = us;
	gp += countBitOnes(mask);

#else
	ushort* pbeg = (ushort*)LeftShiftArray[m1] + (beg1_bit >> 4);
	ushort* pend = (ushort*)RightShiftArray[m2] + (end2_bit >> 4);

	int shiftBits = end2_bit & 15;

	int loops = (end1_bit - beg1_bit) >> 4;
	//asm to hot spot
	while (loops-- > 0) {
		ushort mask16 = *pbeg++ | (*(uint*)pend++ >> shiftBits);
		gp += countBitOnes(mask16 | 0xffff0000);
	}
	ushort us = *pbeg | (*((uint*)pend) >> shiftBits);
	mask = ~((1 << ((1 + end1_bit - beg1_bit) & 15)) - 1);
	if ((ushort)mask != (ushort) ~0)
		mask |= us;
	else
		mask = us;
	gp += countBitOnes(mask | 0xffff0000);
/**
	int loops = (end1_bit - beg1_bit) >> 5;

	//asm to hot spot
	while (loops-- > 0) {
		ushort mask1 = *pbeg++ | (*(uint*)(pend++) >> shiftBits);
		ushort mask2 = *pbeg++ | (*(uint*)(pend++) >> shiftBits);
		gp += WordNumBit0[mask2] + WordNumBit0[mask1];
	}
	uint us = (*((uint*)pbeg)) | (*((uint64*)(pend)) >> (end2_bit & 31));
	mask = ~((1u << ((1 + end1_bit) & 31)) - 1);
	if (mask != (uint) ~0)
		mask |= us;
	else
		mask = us;
	gp += WordNumBit0[(ushort)mask] + WordNumBit0[mask >> 16];
*/
#endif

	return gp;
}
#else
{
	uint gp = 0;

#if PRIME_DIFF == 0
	if (beg1_bit > end1_bit)
		return 0;
	end2_bit = sizeof(RightShiftArray[0][0]) * 8 - end2_bit - 1;
	uint64* pbeg = (uint64*)(LeftShiftArray[m1][beg1_bit & MASK] + (beg1_bit >> MOVE));
	uint64* pend = (uint64*)(RightShiftArray[m2][end2_bit & MASK] + (end2_bit >> MOVE));
#else
	end2_bit = sizeof(RightShiftArray[0]) * 8 - end2_bit - 1;

	int shiftbits = end2_bit & MASK;
	utype* pleft = LeftShiftArray[m1][beg1_bit & MASK] + (beg1_bit >> MOVE);
	utype* pright = RightShiftArray[m2] + (end2_bit >> MOVE);
	//assert((end2_bit >> MOVE) <= sizeof(RightShiftArray[0]) / sizeof(utype));

	//pack the word/dword for align
	utype umask = ~(((utype)1 << (MASK + 1 - shiftbits)) - 1);
	umask |= *pleft | (*pright >> shiftbits);

# if MOVE >= 5
	if (0 == shiftbits)
		umask = *pleft | *pright;
#endif

	if (end1_bit < MASK + beg1_bit)
		umask |= ~(((utype)1 << (end1_bit - beg1_bit + 1)) - 1);

#if MOVE == 6
	gp = countBitOnes(umask >> 32) + countBitOnes((uint)umask);
#elif MOVE == 5
	gp = countBitOnes(umask);
#elif MOVE == 4
	gp = countBitOnes(umask | 0xffff0000);
#else
	gp = countBitOnes(umask | 0xffffff00);
#endif

	beg1_bit += MASK + 1 - shiftbits;
	if (beg1_bit > end1_bit)
		return gp;

	uint64* pbeg = (uint64*)(LeftShiftArray[m1][beg1_bit & MASK] + (beg1_bit >> MOVE));
	uint64* pend = (uint64*)(pright + 1);
#endif

	int bitleng = end1_bit - beg1_bit;

	gp += countZeroBits(pbeg, pend, bitleng);

	pbeg += bitleng >> 6, pend += bitleng >> 6;

	uint64 smask = ~(((uint64)1 << ((1 + bitleng) % (1 << 6) )) - 1);
	if (smask != (uint64) ~0)
		smask |= *pbeg | *pend;
	else
		smask = *pbeg | *pend;

	gp += countBitOnes(smask >> 32) + countBitOnes((uint)smask);

	return gp;
}
#endif

static int countPairMin(const int sum23, const int p1, const int m12, const int cacheGpt)
{
	int gp = 0;
	const int m1 = (ushort)m12;
	const int m2 = m12 >> 16;

	const int sum_bit = (sum23 - Pattern[m1] - Pattern[m2]) / FACTP;
	int beg1_bit = p1 / FACTP;
	if (beg1_bit * FACTP + Pattern[m1] < p1)
		beg1_bit++;

	int end1_bit = (sum23 - p1) / FACTP;
	if (end1_bit * FACTP + Pattern[m1] > sum23 - p1)
		end1_bit--;

	int end2_bit = sum_bit - beg1_bit;
	if (end2_bit * FACTP + Pattern[m2] > sum23 - p1) {
		end2_bit--;
		beg1_bit++;
	}

//	if (end1_bit < 0)
//		return 0;
	//p1 + p1 < (p1 + sum23) - 3*p1
#if TABLE_GP
	if (cacheGpt) {
		gp = -countPairSum(0, beg1_bit - 1, sum_bit, m1, m2);
		if (m1 != m2)
			gp -= countPairSum(end1_bit + 1, sum_bit, sum_bit - end1_bit - 1, m1, m2);
	} else {
		if (m1 == m2)
			end1_bit = (beg1_bit + end1_bit) / 2;
		gp = countPairSum(beg1_bit, end1_bit, end2_bit, m1, m2);
	}
#else
	if (m1 == m2)
		end1_bit = (beg1_bit + end1_bit) / 2;
	gp = countPairSum(beg1_bit, end1_bit, end2_bit, m1, m2);
#endif

	return gp;
}

//
static int countPairMax(const int sum12, const int p3, const int m12, const int cacheGpt)
{
	//beg1_bit * FACTP + Pattern[m1] + end2_bit * FACTP + Pattern[m2] = sum12;
	const int m1 = (ushort)m12;
	const int m2 = m12 >> 16;
	const int sum_bit = (sum12 - Pattern[m1] - Pattern[m2]) / FACTP;

	int beg1_bit = 0, gp = 0;
	int minprime = p3 < sum12 ? p3 : sum12;

	//for creating GPT table and get G(sum12)
	if (Config.CreateGpt)
		minprime = sum12;

	int end1_bit = minprime / FACTP;
	if (end1_bit * FACTP + Pattern[m1] > minprime)
		end1_bit--;

	int end2_bit = end1_bit;
	if (end2_bit * FACTP + Pattern[m2] > minprime)
		end2_bit--;
	beg1_bit = sum_bit - end2_bit;

#if TABLE_GP == 0
	if (beg1_bit > end1_bit)
		return 0;
#else
	//	if (end1_bit < 0)
	//		return 0;
#endif

#if TABLE_GP
	//(3 + 3 / 7) / (2 + 4 / 7) = 4 / 3
	//if (p3 * 4 > 3 * sum12 && GPT[sum12 / 2]) {
	if (cacheGpt) {
		gp = -countPairSum(0, beg1_bit - 1, sum_bit, m1, m2);
		if (m1 != m2)
			gp -= countPairSum(end1_bit + 1, sum_bit, sum_bit - end1_bit - 1, m1, m2);
	} else {
		if (m1 == m2)
			end1_bit = (beg1_bit + end1_bit) / 2;
		gp = countPairSum(beg1_bit, end1_bit, end2_bit, m1, m2);
	}
#else
	if (m1 == m2)
		end1_bit = (beg1_bit + end1_bit) / 2;
	gp = countPairSum(beg1_bit, end1_bit, end2_bit, m1, m2);
#endif

	return gp;
}

//T(n = sum2 + Num) = sigma gp(sum2) with Num is prime
static int countGp1(const int sum2, const int p13)
{
	int gp = 0;

	int cacheGpt = 0;
#if TABLE_GP
	//(p13 + sum2) / BEST_RATION = sum2
	if (Config.EnumMax) {
#if LOAD_SEG == 0
		if (2 * (p13 + sum2) * BEST_RATION > sum2 && GPT[sum2 / 2]) {
			cacheGpt = 1;
			gp = GPT[sum2 / 2];
			if (p13 >= sum2)
				return gp;
		}
#else
		if (2 * (p13 + sum2) * BEST_RATION > sum2) {
			gp = getGpFromTable(sum2 / 2);
			if (gp) {
				cacheGpt = 1;
				if (p13 >= sum2)
					return gp;
			}
		}
#endif
	} else {
		if (p13 * 4 < sum2 && sum2 / 2 < (sizeof(GPT) / sizeof(GPT[0]))
			&& GPT[sum2 / 2]) {
			cacheGpt = 1;
			gp = GPT[sum2 / 2];
		}
	}
#endif

	const int pindex = sum2 % FACTP / 2;
	for (int i = 0; Pairs[pindex][i] >= 0; i++) {
		int m12 = Pairs[pindex][i];
		if (Config.EnumMax)
			gp += countPairMax(sum2, p13, m12, cacheGpt);
		else
			gp += countPairMin(sum2, p13, m12, cacheGpt);
	}

	return gp;
}

//
static int countGp2(const int sum23, const int p1)
{
	int gp = 0;
	int beg_bit = p1 >> 1;
	int end_bit = (sum23 - p1) / 2;

	if (beg_bit + MIN_DIF_BITS > end_bit)
		return getBitsGp(beg_bit, end_bit);

	while (beg_bit & 31) {
		if (!TEST_BIT(SrcBitArray, end_bit) &&
			!TEST_BIT(SrcBitArray, beg_bit))
			gp++;
		beg_bit++; end_bit--;
	}

	int loops = (end_bit - beg_bit) >> 6;
	uint* pbeg = ((uint*)SrcBitArray) + (beg_bit >> 5);

	end_bit = sizeof(RevBitArray) * 8 - end_bit - 1;
	int shiftBits = end_bit & 31;
	uint* pend = ((uint*)RevBitArray) + (end_bit >> 5);

	while (loops-- > 0) {
		uint mask = *pbeg++ | (uint)(*((uint64*)pend) >> shiftBits);
		gp += countBitOnes(mask);
		pend++;
	}

	beg_bit = (int)(pbeg - (uint*)SrcBitArray) << 5;
	end_bit = sum23 / 2 - beg_bit - 1;
	//assert(beg_bit <= end_bit);
	gp += getBitsGp(beg_bit, end_bit);
	return gp;
}

//init number[0 - 2^16] of bit 0s by binary representation
//and word(16 bits) reverse of binary representation
static void initBitTable( )
{
	//init WordNumBit0 table in 0-2^16
	int nbitsize = sizeof(WordNumBit0) / sizeof(WordNumBit0[0]);
	int i = 0;
	WordNumBit0[0] = 0;
	for (i = 1; i < nbitsize; i++)
		WordNumBit0[i] = WordNumBit0[i >> 1] + (i & 1);
//	for (i = 0; i < nbitsize; i++)
//		WordNumBit0[i] = (1 << 4) - WordNumBit0[i];

	//reverse word by binary representation
	uchar bytereverse[256] = {0};
	nbitsize = sizeof(WordReverse) / sizeof(WordReverse[0]);
	for (i = 1; i < (1 << 8); i++)
		bytereverse[i] = reverseByte(i);
	for (i = 1; i < nbitsize; i++)
		WordReverse[i] = bytereverse[i >> 8] | (bytereverse[i & 255] << 8);
}

//init Pairs: Pairs[i] divide into two
//16 bits number p1, p2 and gcd((p1 + p2), FACTP) = 1
static void initPattern(  )
{
	int pas = 0;
	memset(MapIndex, -1, sizeof(MapIndex));

	//init Pattern
	for (int i = 1; i < FACTP; i += 2) {
		int j = 0, ok = true;
		for (; FACTP % SmallPrime[j] == 0; j++) {
			if (i % SmallPrime[j] == 0) {
				ok = false;
				break;
			}
		}
		if (ok) {
			Pattern[pas] = i;
			MapIndex[i] = pas++;
		}
	}

	int pindex[FACTP] = {0};
	memset(Pairs, -1, sizeof(Pairs));
	for (int k = 0; k < PATTERNS; k++) {
		for (int j = k; j < PATTERNS; j++) {
			int sum = (Pattern[k] + Pattern[j]) % FACTP / 2;
			Pairs[sum][pindex[sum]++] = (j << 16) + k;
//			PairMask[sum][j >> 5] |= (1u << (j & 31));
//			PairMask[sum][k >> 5] |= (1u << (k & 31));
		}
	}

#if FACTP == 30
	for (int n = 0; n < FACTP; n += 2) {
//		printf("sum = %2d, mask = %u\n", n, PairMask[n >> 1][0]);
//		printf("%d, ", PairMask[n >> 1][0]);
	}
//	putchar('\n');
#endif

}

//
static void initTnTable( )
{
	eratosSieve(MAXN - 10 * FACTP);

	initBitTable( );

	initPattern( );

	adjustLRevBitArray(MAXN - 10 * FACTP);

	const int rsize = sizeof(RevBitArray) / sizeof(RevBitArray[0]);
	for (int i = 0; i < rsize; i++)
		RevBitArray[rsize - i - 1] = reverseUtype(SrcBitArray[i]);
}

#if TABLE_GP
static bool checkDataType(FILE * file)
{
	int check[2];
	fread(check, sizeof(check[0]), 2, file);
	if (check[0] != sizeof(GPT[0])) {
		puts("data file with wrong size type");
		fclose(file);
		return false;
	}
	return true;
}

//bugs
#if LOAD_SEG
static int getGpFromTable(const int n)
{
	static FILE *file = 0;
	static int currpos = 0;
	static int lastpos = 0;
	static int gntsize = sizeof(GPT) / sizeof(GPT[0]);

	if (file == 0 || n < 0) {

		if (lastpos == 0 && file && n < 0)
			return 0;

		if (file)
			fclose(file);
#if MAXN > 20000000
		file = fopen("primei.ttn", "rb");
#else
		file = fopen("prime.ttn", "rb");
#endif
		if (file == 0 || checkDataType(file) == false)
			return -1;
		if (n < 0) {
			lastpos = 0;
			memset(GPT, 0, sizeof(GPT));
			int nw = fread(GPT + GSINDEX - 1, sizeof(GPT[0]), gntsize - GSINDEX, file);
			currpos = nw + GSINDEX - 1;
			return -1;
		}
	} else if (n >= currpos) {
		lastpos = currpos;
		int nw = fread(GPT, sizeof(GPT[0]), gntsize, file);
		if (nw < gntsize)
			memset(GPT + nw, 0, sizeof(GPT) - nw * sizeof(GPT[0]));
		currpos += nw;
	}
	int gp = GPT[n - lastpos] - getSmallGp(n * 2);
	if (gp < 0 && sizeof(GPT[0]) == 2)
		gp = 0;
	return gp;
}
#endif

//
static bool loadGpTable(const int n)
{
	double ts = getTime();

	int gntsize = sizeof(GPT) / sizeof(GPT[0]);
	if (gntsize > n * BEST_RATION)
		gntsize = n * BEST_RATION;

#if MAXN > 20000000
	FILE *file = fopen("primei.ttn", "rb");
#else
	FILE *file = fopen("prime.ttn", "rb");
#endif

	if (file && checkDataType(file)) {
		fread(GPT + GSINDEX - 1, sizeof(GPT[0]), gntsize, file);
		fclose(file);
#if OMP
		#pragma omp parallel for schedule(dynamic, 100) if (n > 1000000)
#endif
		for (int i = GSINDEX; i < gntsize; i++) {
			int gp = (int)(GPT[i] -= getSmallGp(i * 2));
# if MAXN < 20000000
			if (gp < 0 && sizeof(GPT[0]) == 2)
				GPT[i] = 0;
# endif
		}
		printf("load gp Table time use %.lf ms\n", getTime() - ts);
		return true;
	}

	return false;
}

static int createGpTable(const int n, int offset, int step)
{
	double ts = getTime();

	int gntsize = sizeof(GPT) / sizeof(GPT[0]);
	if (gntsize > n)
		gntsize = n;

	Config.CreateGpt = true;

	//create table GPT
#if OMP
	omp_set_num_threads(Config.Threads);
	#pragma omp parallel for schedule(dynamic, 100) if (n > 1000000)
#endif
	for (int i = GSINDEX + offset; i < gntsize; i += step) {
		int gp = countGp1(i * 2, 0) + getSmallGp(i * 2);
#if MAXN < 20000000
		if (sizeof(GPT[0]) == 2 && gp >= (1 << 16))
			gp = 0;
#endif
		GPT[i] = gp;
		if (i % (1 << 17) == 0)
			printf("finish %d%%\r", 100 * i / n );
	}

	return gntsize;
}

//
static bool createGpTable(const int n)
{
	double ts = getTime();

	int gntsize = sizeof(GPT) / sizeof(GPT[0]);
	if (gntsize > n * BEST_RATION)
		gntsize = n * BEST_RATION;

	printf("wait for creating goldbach partion(%d) table mode %d\n",
		2 * gntsize, FACTP);

	if (gntsize > 1000000 && Config.Threads > 1)
		startWorkThread(Config.Threads, gntsize, (int*)threadProcGp);
	else
		createGpTable(gntsize, 0, 1);

#if MAXN > 20000000
	FILE *file = fopen("primei.ttn", "wb");
#else
	FILE *file = fopen("prime.ttn", "wb");
#endif

	//save table GPT to file
	if (file) {
		//save check the size type
		int check[2] = {sizeof(GPT[0]), gntsize * sizeof(GPT[0]) + 8};
		fwrite(check, sizeof(check[0]), 2, file);
		fwrite(GPT + GSINDEX - 1, sizeof(GPT[0]), gntsize, file);
		fclose(file);
#if OMP
		#pragma omp parallel for if (n > 400000)
#endif
		for (int i = GSINDEX; i < gntsize; i++) {
			if (GPT[i])
				GPT[i] -= getSmallGp(i * 2);
		}
		printf("create G(%d) table time use %.lf s\n",
			gntsize * 2, (getTime() - ts) / 1000);
	} else {
		Config.CreateGpt = false;
		return false;
	}

	Config.CreateGpt = false;
	return true;
}

//
static int loadGpTable2(const int n)
{
	double ts = getTime();
	Config.CreateGpt = true;

	//use too much memory here!!!
	int gntsize = sizeof(GPT) / sizeof(GPT[0]);
	if (gntsize > n * BEST_RATION)
		gntsize = n * BEST_RATION;

	if (freopen("prime.ttn", "rb", stdin)) {
		int cnt = GSINDEX - 1, gp;
		while (scanf("%d\n", &gp) != EOF && cnt++ < gntsize) {
#if MAXN < 20000000
			if (gp < (1 << 16) || sizeof(GPT[0]) == 4)
#endif
				GPT[cnt] = gp - getSmallGp(cnt * 2);
		}

		//restore io to console
		freopen(CONSOLE, "r", stdin);
		Config.CreateGpt = false;
		printf("load gp Table time use %.lf ms\n", getTime() - ts);
		return 1;
	}

	printf("please wait for creating gp Table mod %d\n", FACTP);

#if OMP
	#pragma omp parallel for schedule(dynamic, 100) if (n > 1000000)
#endif

	for (int i = GSINDEX; i < gntsize; i++) {
		int gp = countGp1(i * 2, 0) ;
		int sn = getSmallGp(i * 2);
#if MAXN < 20000000
		if (gp + sn < (1 << 16) || sizeof(GPT[0]) == 4)
#endif
			GPT[i] = gp + sn;
		if (i % (1 << 16) == 0)
			printf("finish %d%%\r", 100 * i / (n * BEST_RATION));
	}

	freopen("primex.tn", "wb", stdout);
	for (int j = GSINDEX; j < n * BEST_RATION; j += 1)
		printf("%d\n", GPT[j]);

	Config.CreateGpt = false;
	//restore io to console
	freopen(CONSOLE, "w", stdout);
	printf("create G(%d) table time use %.lf ms\n",
		2 * n * BEST_RATION, getTime() - ts);

	return 0;
}
#endif

/************************************************************************/
/*                                                                      */
/************************************************************************/
static int64 coreSieve(int n, int offset, int step)
{
	int64 tn = 0;
	int MaxMinNum = n / 3, gp = 0;

	if (Config.ShowLog) {
		printf("%d, %d, %d\n", n, offset, step);
	}

	int sp = 0;
	for (int j = 0; FACTP % SmallPrime[j] == 0; j++)
		sp = SmallPrime[j];


	if (Config.EnumMax) {
/**	#if (OMP)
		omp_set_num_threads(Config.Threads + 1);
		#pragma omp parallel for reduction(+: tn) schedule(dynamic, 2048) if (n > 300000)
	#endif
**/
		for (int p3 = (n - 14) / 2 * 2 + offset; p3 >= MaxMinNum; p3 -= step) {
			if (!TEST_BIT2(SrcBitArray, p3))
				gp = countGp1(n - p3, p3), tn += gp;// printf("%d, %d\n", p3, gp);
#if MAXN > 100000000
			if (p3 % (1 << 19) == 1)
				printf("finish %.2lf%%, tn = %I64d\r", p3 * 100. / n, tn);
#endif
			//printf("n = %d, p3 = %d, tn = %I64d\n", n - p3, p3, tn);
		}
	} else {
/**	#if (OMP)
		omp_set_num_threads(Config.Threads + 1);
		#pragma omp parallel for reduction(+: tn) schedule(dynamic, 2048) if (n > 300000)
	#endif
*/
		for (int p1 = sp + 1 + offset; p1 <= MaxMinNum; p1 += step) {
			if (!TEST_BIT2(SrcBitArray, p1)) {
				gp = countGp1(n - p1, p1); tn += gp;
				//printf("offset = %d, n = %d, p1 = %d, tn = %I64d\n", offset, n - p1, p1, tn);
			}
		}
	}

	return tn;
}

/************************************************************************
  the basic of algorithm :
G(n)
/|\
 |
 |                   +\
 |                   |  \
 |                   |    \
 |                   |    |#\
 |                   |    |## \
 |                   |    |#### \
 |                   |    |###/ |-\
 |                   |    |##/  |---\
 |                   |    |#/   |-----\
 |                   |    |/    |-------\
 |                   |    /     |---------\
 |                   |   /|     |-----------\
 |                   |  /.|     |-------------\
 |                   | /..|     |---------------\
 |                   |/...|     |-----------------\
 |                   .....|     |-------------------\
 |-------------------------------------------------------------------> maxn
 0                  N/3  3N/7   N/2                  N

- mean use table to get G(n) with maxn > N / 2
# mean use table and calculation to get G(n) with maxn > 3N / 7
. mean use calculation to get G(n) with maxn > N / 3
************************************************************************/

static int64 getTn(const int n)
{
	assert(n % 2 == 1 && n < (MAXN - 10 * FACTP));
	int64 tn = 0;

	for (int j = 0, sp = 0; FACTP % SmallPrime[j] == 0; j++) {
		sp = SmallPrime[j];
		tn += countGp2(n - sp, sp);
	}

	double ts = getTime();

#if LOAD_SEG
	getGpFromTable(-1);
#endif

	if (n > 1000000 && Config.Threads > 1)
		tn += startWorkThread(Config.Threads, n, (int*)threadProcTn);
	else
		tn += coreSieve(n, 1, 2);

	//print the output
	if (Config.ShowRet) {
		printf("T(%d) = %I64d", n, tn);
		if (Config.ShowTime) {
			if (n < 20000000)
				printf(", time use %.lf ms", getTime() - ts);
			else
				printf(", time use %.2lf s", (getTime() - ts) / 1000.);
		}
		putchar('\n');
	}

	return tn;
}

//test case code, test data from third party

//convert str to int64 ((x)*E(y)+/-(z))
//e9 2e10+-2   10^11-25
static uint atouint(const char* str)
{
	uint ret = 0;

	while (isspace(*str)) {
		str++;
	}

	while (isdigit(*str)) {
		ret = ret * 10 + *str++ - '0';
	}

	if (*str && isdigit(str[1])) {
		if (*str == '^') {
			ret = (uint)(pow((int)ret, (double)atoi(str + 1)) + 0.01);
		} else if (*str == 'e' || *str == 'E') {
			if (ret == 0) {
				ret = 1;
			}
			ret *= (uint)(pow((int)10, (double)atoi(str + 1)) + 0.01);
		}
	}

	const char* ps = str;
	if (ps = strchr(str, '+')) {
		ret += atoi(ps + 1);
	} else if (ps = strchr(str, '-')) {
		ret -= atoi(ps + 1);
	} else if (ps = strchr(str, '*')) {
		ret *= atoi(ps + 1);
	}

	return ret;
}

#define TEST_FILEDATA 1
static int startTest(int t = 100000)
{
	srand((uint)time(NULL));
	double ts = getTime();
	Config.ShowRet = false;

#if TEST_FILEDATA
	if (!freopen("prime.tn", "rb", stdin)) {
		puts("can not read test data file");
		Config.ShowRet = true;
		return -1;
	}
#else
	if (!freopen("prime.tn", "wb", stdout)) {
		puts("can not write test data file");
		Config.ShowRet = true;
		return -2;
	}
#endif

	const char* sformat1 = "%d T(%d) = %I64d\n";
	const char* sformat2 = "T(%d) = %I64d\n";
	for (int i = 1; i <= t; i++) {
#if TEST_FILEDATA == 0
		uint n = rand() * rand() % (MAXN / 10);
		n = (n + 4) * 2 + 1;
		int64 tn = getTn(n);
		printf(sformat1, i, n, tn);
#else
		int index = 0, n;
		int64 res;
		if (scanf(sformat1, &index, &n, &res) != 3 &&
			scanf(sformat2, &n, &res) != 2) {
			printf("read case %d with wrong data\n", i);
			break;
		}
		int64 tn = getTn(n);
		if (tn != res)
			printf("case %d with wrong data T(%d) = %I64d != getTn(%d) = %I64d\n",
					i, n, res, n, tn);
		if ((i & 63) == 0)
			printf("case pass %d%%\r", i * 100 / t);
#endif
	}
	printf("test case time use %.lf ms\n", getTime() - ts);

	Config.ShowRet = true;
	freopen(CONSOLE, "w", stdout);
	freopen(CONSOLE, "r", stdin);

	return 0;
}

//
static void unitTest( )
{
	const char* tndata[] =
	{
		"T(1000001) = 104528645",
		"T(2000001) = 235945824",
		"T(3000001) = 724870542",
		"T(4000001) = 1210277117",
		"T(5000001) = 1200084673",
		"T(6000001) = 2411437969",
		"T(7000001) = 3285264440",
		"T(8000001) = 2779050013",
		"T(9000001) = 5145558388",
		"T(10000001)= 6154337544"
	};

	for (int i = 0; i < sizeof(tndata) / sizeof(tndata[0]); i++) {
		int n;
		int64 tab;
		if (sscanf (tndata[i], "T(%d) = %I64d", &n, &tab) == 2) {
			int64 cal = getTn(n);
			if (cal != tab)
				printf ("T(%d) = %I64d != %I64d (cal) \n", n, tab, cal);
		}
	}
}

//list gp by the input value start, end, step
static void listGp(const char params[][80])
{
	double ts = getTime();

	int start = 100000000, end = start * 10 / 10 + 1000, step = 2;
	int ni = 1;
	int buf[ ] = {0, start, end, step, 0};

	for (int i = 0; params[i][0] && ni < sizeof(buf) / sizeof(buf[0]); i++) {
		if (isdigit(params[i][0]) || params[i][0] == 'e')
			buf[ni++] = atouint(params[i]);
	}

	start = buf[1], end = buf[2], step = buf[3];

	start += start & 1;
	step += step & 1;

	if (step < 2 || step > 100000000)
		step = 2;
	if (start >= end)
		end = start + end * step - 1;
	if (end > MAXN * BEST_RATION * 2)
		end = MAXN * BEST_RATION * 2 - 20000;

	printf("%d:%d:%d\n", start, (end - start + 2) / 2, step);

	if (Config.FileSave) {
		Config.ShowTime = true;
		freopen("prime.gp", "wb", stdout);
	}

	Config.CreateGpt = true;
	int cnt = 1, gp = 0;
	for (int n = start; n <= end; n += step) {
		if (isdigit(params[5][0]))
			printf("%d ", cnt++);
#if TABLE_GP
		if (n / 2 < sizeof(GPT) / sizeof(GPT[0]) && GPT[n / 2])
			gp = GPT[n / 2] + getSmallGp(n);
		else
			gp = countGp1(n, 0) + getSmallGp(n);
#else
		gp = countGp1(n, 0) + getSmallGp(n);
#endif
		if (Config.ShowTime)
			printf("G(%d) = %d\n", n, gp);
	}

	if (Config.FileSave) {
		freopen(CONSOLE, "w", stdout);
	}
	Config.CreateGpt = false;
	printf("gp case time use %.lf ms\n", getTime() - ts);
}

//list Tn by the input value start, end, step
static void listTn(const char params[][80])
{
	double ts = getTime();

	int start = 10000000, end = start * 10 / 10 + 100, step = 2;
	int ni = 1;
	int buf[ ] = {0, start, end, step, 0};

	for (int i = 0; params[i][0] && ni < sizeof(buf) / sizeof(buf[0]); i++) {
		if (isdigit(params[i][0]) || params[i][0] == 'e')
			buf[ni++] = atouint(params[i]);
	}

	start = buf[1], end = buf[2], step = buf[3];

	start += 1 - (start & 1);
	step += step & 1;

	if (step < 2 || step > 10000000)
		step = 2;
	if (start > end && end < 10000000)
		end = start + end * step - 1;

	printf("%d:%d:%d\n", start, step, end);

	if (Config.FileSave) {
		Config.ShowTime = false;
		freopen("prime.tn", "wb", stdout);
	}

	int cnt = 1;
	for (int n = start; n <= end; n += step) {
		if (isdigit(params[5][0]))
			printf("%d ", cnt++);
		getTn(n);
	}
	if (Config.FileSave)
		fclose(stdout);
	printf("TN case time use %.lf ms\n", getTime() - ts);
}

//print the benchmark form startn to MAXN
static void benchMark(const int startn)
{
	Config.ShowRet = true;
	int start = startn, gap = start - 1;
	for (int i = 0; start * 10 < MAXN; i++) {
		for (int j = 0; j < 9; j++) {
			getTn(j * gap + start);
		}
		start = start / 10 * 100 + 1;
		gap *= 10;
	}
	getTn(start);
}

//
static void printInfo( )
{
	puts("---------------------------------------------------------------");
	puts("Count number of tuples T(n): (p1, p2, p3)\n,\
n = p1 + p2 + p3(p1 <= p2 <= p3), p1, p2, p3\nare all odd prime numbers");
	puts("version 6.2 Copyright (c) by Huang Yuanbing 2009 - 2011 bailuzhou@163.com");

	getCpuInfo();

	printf("[MARCO] : FACTP = %d, LOAD_SEG = %d, TABLE_GP = %d, POPCNT = %d\n",
			FACTP, LOAD_SEG, TABLE_GP, POPCNT);
	printf("[MARCO] : PRIME_DIFF = %d, Work thread = %d\n", PRIME_DIFF, Config.Threads);

	const int mems = (1000 << 10) + sizeof(SrcBitArray) * 2 +
	sizeof(LeftShiftArray) + sizeof(RightShiftArray)
#if TABLE_GP
		+ sizeof(GPT)
#endif
		+ sizeof(WordReverse) + sizeof(WordNumBit0) + sizeof(Pairs) + sizeof(Pattern);

#ifdef _MSC_VER
	printf("Compiled by MS/vc++ %d", _MSC_VER);
#else
	printf("Compiled by GNU/g++ %d.%d.%d",
		__GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#endif

#if _M_AMD64 || __x86_64__
	printf(" on Windows 64 bit");
#endif
	printf(" on %s %s\n", __TIME__, __DATE__);

	printf("total physical memory use ~= %.3lf M\n", mems / (1024 * 1024.0));
	puts("--------------------------------------------------------------\n");
}

//
static int setConfig(char params[][80])
{
	int cmdi = -1;

	for (int i = 0; params[i][0]; i++) {
		char c = params[i][0];
		if (c >= 'a' && c <= 'z')
			c += 'A' - 'a';
		if (isdigit(c) || c == 'E') {
			if (cmdi < 0)
				cmdi = i;
			continue;
		}

		switch(c)
		{
			case 'U':
				Config.TestCase = atoi(params[i] + 1);
				cmdi = i;
				break;
			case 'A':
				Config.EnumMax = !Config.EnumMax;
				break;
			case 'S':
				Config.FileSave = !Config.FileSave;
				break;
			case 'P':
				Config.ShowTime = !Config.ShowTime;
				break;
			case 'D':
				Config.ShowLog = !Config.ShowLog;
				break;
			case 'T':
				Config.Threads = atoi(params[i] + 1);
				if (Config.Threads > MAX_THREADS || Config.Threads < 1)
					Config.Threads = 4;
				break;
			default:
				cmdi = i;
				break;
		}
	}

	return cmdi;
}

//split ccmd to params array
static int splitParams(const char* ccmd, char params[][80])
{
	int nwords = 0;

	for (int i = 0;  ; i++) {
		while (isspace(*ccmd)) {
			ccmd++;
		}
		if (*ccmd == 0 || *ccmd == ';') {
			break;
		}
		char* pc = params[i];
		char c = *ccmd;
		while (isalnum(c) || c == '^' ||
				c == '+' || c == '-' || c == '*') {
			*pc++ = c;
			c = *++ccmd;
		}
		nwords++;
	}

	return nwords;
}

//
static bool excuteCommand(const char* cmd)
{
	while (cmd) {
		char* pcmd = (char*)strchr(cmd, ';');
		char params[10][80] = {0};

		if (splitParams(cmd, params) <= 0) {
			return false;
		}

		int cmdi = setConfig(params);

		if (cmdi < 0) {
			return true;
		}

		const char cmdc = toupper(params[cmdi][0]);

		if (cmdc == 'H') {
			puts(HelpInfo);
			puts(HelpCommand);
		} else if (cmdc == 'B') {
			puts("-----------start benchmark ----------");
			benchMark(1000001);
		} else if (cmdc == 'C') {
#if TABLE_GP
			puts("----------clear g(n) table-----------");
			memset(GPT, 0, sizeof(GPT) );
#endif
		} else if (cmdc == 'U') {
			puts("----------start unit test -----------");
			startTest(Config.TestCase);
		} else if (cmdc == 'L') {
			puts("-----------start list T(n) ----------");
			listTn(params);
		} else if (cmdc == 'G') {
			puts("-----------start list G(n) ----------");
			listGp(params);
		} else if (cmdc == 'Z') {
			return false;
		} else if (isdigit(params[cmdi][0]) ||
					toupper(params[cmdi][0]) == 'E') {
			puts("------------start get T(n) ----------");
			int n = atouint(params[cmdi]);
			getTn(n / 2 * 2 + 1);
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
		printInfo( );
		puts(HelpInfo);
	}

	initTnTable( );

#if TABLE_GP && LOAD_SEG == 0
	if (!loadGpTable(MAXN - 10 * FACTP))
		createGpTable(MAXN - 10 * FACTP);
	//	loadGpTable2(MAXN - 10*FACTP);
#endif

	if (argc > 1 && argv[1][0] == 'u') {
		unitTest();
	}

	for (int i = 1; i < argc; i++) {
		excuteCommand(argv[i]);
	}

	excuteCommand("1e7+1");

	char ccmd[256] = {0};
	while (true) {
		printf("\n[input command] : ");
		if (!gets(ccmd) || !excuteCommand(ccmd))
			break;
	}

	return 0;
}

/*********************************************
T(10000001) = 6154337544, time use 0.58 s
T(20000001) = 13954099494, time use 1.33 s
T(30000001) = 44959615571, time use 4.22 s
T(40000001) = 75674918752, time use 7.22 s
T(50000001) = 75359294527, time use 7.48 s
T(60000001) = 157860860173, time use 15.50 s
T(70000001) = 208780970076, time use 20.98 s
T(80000001) = 177564511381, time use 18.50 s
T(90000001) = 316917615110, time use 33.28 s
T(100000001) = 398220733411, time use 41.92 s

INTEL E2180 2.0G
T(1000001) = 104528645, time use 15 ms
T(2000001) = 235945824, time use 63 ms
T(3000001) = 724870542, time use 140 ms
T(4000001) = 1210277117, time use 219 ms
T(5000001) = 1200084673, time use 219 ms
T(6000001) = 2411437969, time use 437 ms
T(7000001) = 3285264440, time use 578 ms
T(8000001) = 2779050013, time use 485 ms
T(9000001) = 5145558388, time use 890 ms
T(10000001) = 6154337544, time use 1063 ms
T(100000003) = 399782217089, time use 707 s
T(100000005) = 238111692025, time use 450 s

T5500 1.66G
T(100000001) = 398220733411, time use 801 s mod 2310
T(100000001) = 398220733411, time use 912 s mod 210
T(1000000001)= 25759757327768, time use 34371 s mod 210
T(1000000001)= 25759757327768, time use 53079/2 s mod 2310

T(334639305) = 2105252805433, time use 2787.50 s
T(324939615) = 1998576058582, time use 2616.47 s

T(557732175) = ?

Rp(324939567) ＝ 2296838469750
Rp(324939569) ＝ 3445409598643
Rp(324939571) ＝ 3406411282972
Rp(324939573) ＝ 2222554318680
Rp(324939575) ＝ 3180381493504
Rp(324939577) ＝ 3429181266765
Rp(324939579) ＝ 2296741535721

todo:
1. Remove marco MAXN, FACTP, PATTERNS
2.
**********************************************/

