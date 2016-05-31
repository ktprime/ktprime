# include <stdio.h>
# include <memory.h>
# include <stdlib.h>
# include <ctype.h>
# include <string.h>
# include <time.h>
# include <assert.h>
# include <limits.h>

# define WHEEL           30
# define WHEEL_SKIP      0x799b799b

# define L1_CACHE_SIZE   63
# define L1_SIEVE_SEG    4

//max continuous goldbach partition number
# define MAX_GPCOUNT     100000
# define MAX_THREADS     16
# define OMP             0

# define BIT_SCANF       1
# define PRIME_DIFF      0

//use of the SSE4.2/ SSE4a POPCNT instruction for fast bit counting.
#if _MSC_VER > 1300
	# define POPCNT      1
	# include <intrin.h>
#elif (__GNUC__ * 10 + __GNUC_MINOR__ > 44)
	# define POPCNT      1
	# include <popcntintrin.h>
#else
	# define POPCNT      0
#endif

#ifdef _MSC_VER
	#define MEM_ALIGN(n) //__declspec(align(n))
#else
	#define MEM_ALIGN(n) __attribute__ ((aligned(n)))
#endif
//	#pragma pack (16)

#if OMP
	#include <omp.h>
#endif

#if defined _M_AMD64
	# define ASM_X86     0
#elif _MSC_VER
	# define ASM_X86     1
#else
	# define ASM_X86     0
#endif

typedef unsigned char uchar;
typedef unsigned short ushort;
typedef unsigned int uint;

#ifdef _WIN32
	static const char* GPFORMAT = "G(%I64d) = %I64d\n";
	typedef __int64 int64;
	typedef unsigned __int64 uint64;
	# include <windows.h>
	# define CONSOLE "CON"
#else
	static const char* GPFORMAT = "G(%llu) = %llu\n";
	typedef long long int64;
	typedef unsigned long long uint64;
	# include <sys/time.h>
	# include <pthread.h>
	# include <unistd.h>
	# define CONSOLE "/dev/tty"
#endif

//intel 4, amd 5 > 3 > 4
# define BSHIFT 4
# if BSHIFT == 3
	typedef uchar utype;
# elif BSHIFT == 4
	typedef ushort utype;
# elif BSHIFT == 5
	typedef uint utype;
# endif

#if __x86_64__ || _M_AMD64 || __amd64__
	#define X86_64   1
#endif

# if PRIME_DIFF
	# define NEXT_PRIME(p, j) p += Prime[++j]
# else
	# define NEXT_PRIME(p, j) p = Prime[++j]
#endif

#if BIT_SCANF == 0
	#define PRIME_OFFSET(mask)   LeftMostBit1[mask]
	typedef ushort stype;
#else
	#define PRIME_OFFSET(mask)   bitScanForward(mask)
#if X86_64
	typedef uint64 stype;
#else
	typedef uint stype;
#endif
#endif

static struct
{
	int L1Size;
	int L1Maxp;
	int L1Index;
	int L2Size;
}
CpuCache =
{
	(L1_CACHE_SIZE << 10) * WHEEL,
	(L1_CACHE_SIZE << 10) / L1_SIEVE_SEG,
	0,
	256
};

enum
{
	FACTP        = 510510 * 1,
	BLOCKS       = 6,
	FIRST_SIEVE  = 19,
	BLOCK_SIZE   = FACTP * BLOCKS,
	BUFFER_SIZE  = BLOCK_SIZE / WHEEL,
	BUFFER_SIZE8 = (BUFFER_SIZE / 8) - (BUFFER_SIZE / 8) % 32,
};

# define MASK_N(n)         1 << (n & ((1 << BSHIFT) - 1))
# define SET_BIT(a, n)     a[(n) >> BSHIFT] |= MASK_N(n)
# define CLR_BIT(a, n)     a[(n) >> BSHIFT] &= ~MASK_N(n)
# define FLP_BIT(a, n)     a[(n) >> BSHIFT] ^= MASK_N(n)
# define TST_BIT(a, n)     (a[(n) >> BSHIFT] & MASK_N(n))

static const char* CmdInfo = "\
	[A: Advanced thread algorithm]\n\
	[S: Set sieve size S[16 - 512]]\n\
	[G: Algorithm choose G[1 - 3]]\n\
	[T: Threads number T[2 - 64]]\n\
	[M: Monitor progress M[0 - 30]]\n\
	[P: Print time/result/progress P[t/r/g]]\n\
	[F: Save result to file]\n\
	[B: Benchmark start from 1e9 B[1 - 20000]]\n\
	[C: Set L1 cache size C[16-128]]\n\
	[U: Unit test U [count] [cases]]\n\
	[R: Single g(n) R [start] [count]]";

static const char* HelpCommand = "\n\
	example: command/param as follow:\n\
	B, B 10000\n\
	C63 S252\n\
	G1-3, T2-32\n\
	P[t, r, g]\n\
	U 1000, U 1000 10*10\n\
	M15\n\
	R 1e10 100\n\
	A R pr pg 100\n\
	2^31-100 10000, 2e10*3 3e3\n\
	10^12 100, 400000000+100 2000";

enum DATA_RESULT
{
	CON_VTOWHEEL = 1,
	COPY_BYBIT,
	SAVE_PRIMEDIF
};

//global config
static struct
{
	bool PrintRet;
	bool PrintTime;
	bool PrintGppair;
	bool SaveResult;
	bool Advanced;

	//algorithm in 1 - 3
	int Algorithm;
	int Threads;
	int SieveSize;
	int PrintGap;
	int Maxp;
	//goldbach partition mask
	ushort GpMask;
} Config =
{
	false, true, false, false, true,
	3, 4, 4 * L1_CACHE_SIZE * (WHEEL << 10),
	(1 << 16) - 1, 1, 0xffff
};

//small prime buffer up to 1e7
#if PRIME_DIFF
	static uchar Prime[664579 + 1200];
#else
	static uint Prime[664579 + 1200];
#endif

//SievedTpl: cross out the first 7/8th Prime's multiple
static uchar SievedTpl[BUFFER_SIZE];
static uchar SievedTpl8[8][(2 * BUFFER_SIZE) / 8 + 32];

#if BIT_SCANF == 0
//bit 1 left most table
static char LeftMostBit1[1 << 16];
#endif

//15 - 8 sgi | 7 - 4 step |0 - 3 wi;
static ushort WheelSkip[WHEEL][8];

//number of bits 1 binary representation table in Range[0-2^16)
#if POPCNT == 0
static uchar WordNumBit1[1 << 16];
#endif

//WordReverse[i] is equal to the bit reverse of i (i < 2^16)
static ushort WordReverse[1 << 16];

//map 8-bit char to 30 bit integer number
static uint Map16To30[1 << 16];

//map wheel byte
static uchar PatternMask[1 << 16];

//G(i) is goldbach partition of start + 2*i
static uint64 GP[MAX_GPCOUNT + WHEEL + 2];

//max number of goldbach partition of start
static const uint64 Maxn = (uint64)1e14;

//the crossing out bit mod 30, the first
//16 bit of SievedTpl map to
//----------------------------------------
//|01/1|07/1|11/1|13/1|17/1|19/1|23/1|29/1| = 0x1111 1111 = SievedTpl[0]
//----------------------------------------
//|31/1|37/1|41/1|43/1|47/1|49/0|53/1|59/1| = 0x1101 1111 = SievedTpl[1]
//----------------------------------------

static const uchar Pattern[ ] =
{
	1,  7,  11, 13, 17, 19, 23, 29,
	31, 37, 41, 43, 47, 49, 53, 59
};

//index of Pattern map to range[0-29]
//WheelIndex[Pattern[i]] = i;
static const int WheelIndex[ ] =
{
	-1, 0, -1, -1, -1,-1,
	-1, 1, -1, -1, -1, 2,
	-1, 3, -1, -1, -1, 4,
	-1, 5, -1, -1, -1, 6,
	-1,-1, -1, -1, -1, 7
};

static const uchar WheelMask[ ] =
{
	255, 1, 255, 255, 255,255,
	255, 2, 255, 255, 255, 4,
	255, 8, 255, 255, 255, 16,
	255, 32,255, 255, 255, 64,
	255,255, 255, 255,255, 128
};

//current index map to closest Pattern index
//WheelLeng[i] : number of Pattern in Range[0, i]
static const uchar WheelLeng[ ] =
{
	0, 0, 1, 1, 1, 1,
	1, 1, 2, 2, 2, 2,
	3, 3, 4, 4, 4, 4,
	5, 5, 6, 6, 6, 6,
	7, 7, 7, 7, 7, 7
};

static const uint GpMask[] =
{
	0xffffffff, 0x29292929, 0x54545454,
	0xfafafafa, 0x23232323, 0xd4d4d4d4,
	0xedededed, 0x0b0b0b0b, 0xd0d0d0d0,
	0xb7b7b7b7, 0x2b2b2b2b, 0xc4c4c4c4,
	0x5f5f5f5f, 0x2a2a2a2a, 0x94949494,
};

static struct ThreadInfo
{
	int gpcount;
	uint64 minn;
	uint64 first;
	uint64* gp;
	uint64 gpbuf[MAX_GPCOUNT];
} Tdata[MAX_THREADS];

static int setSieveSize(uint sievesize);
static void getGp2(uint64 minn, int gpcount, uint64 gp[]);
static void coreSieve1(uint64 minn, int gpcount, uint64 gp[]);
static void coreSieve2(uint64 minn, int gpcount, uint64 first, uint64 gp[]);
typedef void (*sieve_func)(const uint64 , const uint64 , const uint , const int , uint64 []);

static double getTime( )
{
#ifdef _WIN32
	LARGE_INTEGER freq, count;
	QueryPerformanceFrequency(&freq);
	QueryPerformanceCounter(&count);
	return 1000. * count.QuadPart / (double)freq.QuadPart;
#else
	struct timeval tmVal;
	gettimeofday(&tmVal, NULL);
	return tmVal.tv_sec * 1000. + tmVal.tv_usec / 1000.;
#endif
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
		n /= 2;
	}

	return result;
}

//convert str to uint64 < MAXN
//num1[+-*num2], [num1]Enum2[+-*mum3], num1^mum2[+-*num3]
static uint64 atoint64(const char* str, uint64 defaultn)
{
	uint64 ret = 0;

	while (isspace(*str)) {
		str++;
	}

	if (!isdigit(*str) && !(toupper(*str) == 'E')) {
		return defaultn;
	}

	while (isdigit(*str)) {
		ret = ret * 10 + *str++ - '0';
	}

	if (*str && isdigit(str[1])) {
		if (toupper(*str) == 'E') {
			if (ret == 0) {
				ret = 1;
			}
			ret *= ipow(10, atoi(str + 1));
		} else if (*str == '^') {
			ret = ipow(ret, atoi(str + 1));
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

#if BIT_SCANF
inline static uint bitScanForward(stype n)
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

static uint isqrt(const uint64 x)
{
	uint s = 1;
	for (int i = 1; i < 64; i++) {
		if (0 == ((x - 1) >> i)) {
			s = i - 1;
			break;
		}
	}

	uint64 g0 = (uint64)1 << s;
	uint64 g1 = (g0 + (x >> s)) >> 1;

	while (g1 < g0) {
		g0 = g1;
		g1 = (g0 + (x / g0)) >> 1;
	}

	return (uint)g0;
}

static void cpuid(int cpuinfo[4], int id)
{
#if _MSC_VER > 1200
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

//http://msdn.microsoft.com/en-us/library/hh977022%28v=vs.110%29.aspx
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
	} else {
		CpuCache.L1Size = 32 * (WHEEL << 10);
	}
	CpuCache.L1Maxp = CpuCache.L1Size / (WHEEL * L1_SIEVE_SEG);
	CpuCache.L2Size = cpuinfo[2] >> 16;

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
	return Config.Threads;
#else
	Config.Threads = sysconf(_SC_NPROCESSORS_CONF);
#endif

	return Config.Threads;
}

#if _WIN32
static int WINAPI threadProc(void* ptinfo)
#else
static void* threadProc(void* ptinfo)
#endif
{
	struct ThreadInfo* pThreadInfo = (struct ThreadInfo*)ptinfo;
	if (Config.Advanced)
		coreSieve2(pThreadInfo->minn, pThreadInfo->gpcount, pThreadInfo->first, pThreadInfo->gpbuf);
	else
		coreSieve1(pThreadInfo->minn, pThreadInfo->gpcount, pThreadInfo->gp);
	return 0;
}

static int divideThreadData(int threads, uint64 minn, int gpcount)
{
	int i;

	int tsize = gpcount * 2 / threads;
	tsize += tsize & 1;

	Tdata[0].minn = minn;
	Tdata[0].gp = GP;

	for (i = 1; i < threads; i++) {
		Tdata[i].gp = Tdata[i - 1].gp + tsize / 2;
		Tdata[i].first = Config.SieveSize * i;
	}

	for (i = 0; i < threads; i++) {
		if (!Config.Advanced) {
			Tdata[i].gpcount = tsize / 2;
			if (i > 0)
				Tdata[i].minn = Tdata[i - 1].minn + tsize;
		} else {
			Tdata[i].gpcount = gpcount;
			Tdata[i].minn = minn;
		}
	}

	//the last task must not overflowed
	if (!Config.Advanced && Tdata[threads - 1].minn +
		Tdata[threads - 1].gpcount * 2 != minn + 2 * gpcount) {
		Tdata[threads - 1].gpcount = (minn + 2 * gpcount - Tdata[threads - 1].minn) >> 1;
	}

	return threads;
}

static void startWorkThread(int threads, uint64 minn, int gpcount)
{
	int i;
	divideThreadData(threads, minn, gpcount);

#ifdef _WIN32
	HANDLE thandle[MAX_THREADS];
	DWORD tid[MAX_THREADS];
	for (i = 0; i < threads; i++) {
		thandle[i] = CreateThread(NULL, 0,
			(LPTHREAD_START_ROUTINE)threadProc, (LPVOID)(&Tdata[i]), 0, &tid[i]);
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
		pthread_create(&tid[i], NULL, threadProc, &Tdata[i]);
	}
	for (i = 0; i < threads; i++) {
		pthread_join(tid[i], NULL);
	}
#endif

	if (Config.Advanced && Tdata[0].gpbuf) {
		for (i = 0; i < gpcount; i++) {
			GP[i] = 0;
			for (int t = 0; t < threads; t++)
				GP[i] += Tdata[t].gpbuf[i];
		}
	}
}

//reverse bit of a byte with binary representation
//((c * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
static uchar reverseByte(const uchar c)
{
	uchar n =
		(c & 0x55) << 1 | (c & 0xAA) >> 1;
	n = (n & 0x33) << 2 | (n & 0xCC) >> 2;
	n = (n & 0x0F) << 4 | (n & 0xF0) >> 4;
	return n;
}

//set bitarray1 bit position from s1 to sievesize by step
//set bitarray2 bit position from s2 to sievesize by step
inline static void
set2BitArray2(utype bitarray1[], utype bitarray2[], int s1, int s2, const int step, const int sievesize)
{
	for (; s2 < sievesize; ) {
		SET_BIT(bitarray1, s1); s1 += step;
		SET_BIT(bitarray2, s2); s2 += step;
	}
	if (s1 < sievesize) {
		SET_BIT(bitarray1, s1);
	}
}

#if ASM_X86 == 1
	#define ESP_OFFSET 8
	_declspec(naked)
#endif
inline static void
crossOutFactor2(uchar* ps0, uchar* ps1, const uchar* pend, const ushort wordmask, const int step)
{
#if ASM_X86
	__asm
	{
		push esi
		push edi

		mov edi, dword ptr [esp + ESP_OFFSET +  4] //ps0
		mov esi, dword ptr [esp + ESP_OFFSET +  8] //ps1
		mov edx, dword ptr [esp + ESP_OFFSET + 12] //pend
		mov eax, dword ptr [esp + ESP_OFFSET + 16] //wordmask
		mov ecx, dword ptr [esp + ESP_OFFSET + 20] //step

		jmp LCMP2

LOOP2:
		or byte ptr [edi], al
		or byte ptr [esi], ah
		add edi, ecx
		add esi, ecx
LCMP2:
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

inline static int countBit1(uint64 n)
{
#if POPCNT
	//popcnt instruction : INTEL i7/SSE4.2, AMD Phonem/SSE4A
	#if X86_64
	return _mm_popcnt_u64(n);
	#else
	return _mm_popcnt_u32(n) + _mm_popcnt_u32(n >> 32);
	#endif
#elif __GNUC__
	return __builtin_popcountll(n);
#elif 1
	const uint hig = (uint)(n >> 32), low = (uint)n;
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

//count number of bit 0 in binary representation of array
static int countBit0Array(uint64 bitarray[], const int bitleng)
{
	int bit1s = 0;
	int loops = bitleng >> 6;

	while (loops-- >= 0) {
		bit1s += countBit1(*bitarray++);
	}

	return ((1 + (bitleng >> 6)) << 6) - bit1s;
}

//only fast with -O3 on 32 bit
static int countBit0ArrayOr32(uint bitarray1[], uint bitarray2[], int bitleng)
{
	int bit1s = 0;
	int loops = bitleng >> 5;
	const uint mask7 = 0x77777777;

	while (loops-- >= 0) {
#if 1
		uint n = *bitarray1++ | *bitarray2++;
		n -= ((n >> 1) & 0x55555555);
		n = (n & 0x33333333) + ((n >> 2) & 0x33333333);
		n = (n + (n >> 4)) & 0x0F0F0F0F;
		n += n >> 8;
		n += n >> 16;
		bit1s += (n & 0x0000003F);
#else
		uint n = *bitarray1++ | *bitarray2++;
		uint tmp = (n >> 1) & mask7;
		n -= tmp;
		tmp = (tmp >> 1) & mask7;
		n -= tmp;
		tmp = (tmp >> 1) & mask7;
		n -= tmp;
		n = (n + (n >> 4)) & 0x0F0F0F0F;
		n *= 0x01010101;
		bit1s += (n >> 24);
#endif
	}

	return ((1 + (bitleng >> 5)) << 5) - bit1s;
}

//fast on 64 bit OS
static uint countBit0ArrayOr64(const uint64* data1, const uint64* data2, const int bitleng)
{
	const uint64 m1  = 0x5555555555555555;
	const uint64 m2  = 0x3333333333333333;
	const uint64 m4  = 0x0F0F0F0F0F0F0F0F;
	const uint64 m8  = 0x00FF00FF00FF00FF;
	const uint64 m16 = 0x0000FFFF0000FFFF;
	const uint64 h01 = 0x0101010101010101;

	uint bitCount = 0;
	uint64 count1, count2, half1, acc;
	const uint size = bitleng >> 6;
	const uint limit30 = size - size % 30;

	// 64-bit tree merging (merging3)
	for (int i = 0; i < limit30; i += 30) {
		acc = 0;
		for (uint j = 0; j < 10; j ++) {
			count1  = *data1++ | *data2++;
			count2  = *data1++ | *data2++;
			half1   = *data1++ | *data2++;
			count1 -= (count1 >> 1) & m1;
			count2 -= (count2 >> 1) & m1;
			count1 += half1 & m1;
			count2 += (half1 >> 1) & m1;
			count1 = (count1 & m2) + ((count1 >> 2) & m2);
			count1 += (count2 & m2) + ((count2 >> 2) & m2);
			acc    += (count1 & m4) + ((count1 >> 4) & m4);
		}
		acc = (acc & m8) + ((acc >>  8)  & m8);
		acc = (acc       +  (acc >> 16)) & m16;
		acc =  acc       +  (acc >> 32);
		bitCount += (uint)acc;
	}

	// count the bits of the remaining bytes (MAX 29*8) using
	// "Counting bits set, in parallel" from the "Bit Twiddling Hacks",
	// the code uses wikipedia's 64-bit popcount_3() implementation:
	// http://en.wikipedia.org/wiki/Hamming_weight#Efficient_implementation
	for (uint j = 0; j < size - limit30; j++) {
		uint64 x = *data1++ | *data2++;
		x =  x       - ((x >> 1)  & m1);
		x = (x & m2) + ((x >> 2)  & m2);
		x = (x       +  (x >> 4)) & m4;
		bitCount += (uint)((x * h01) >> 56);
	}

	bitCount += countBit1(*data1 | *data2);
	return ((1 + (bitleng >> 6)) << 6) - bitCount;
}

//total time 6.3 / 9.4 = 67 %
static int countBit0ArrayOrPopcnt(uint64 bitarray1[], uint64 bitarray2[], const int bitleng)
{
	int bit1s = 0;
	int loops = bitleng >> 6;

	while (loops-- >= 0) {
		bit1s += countBit1(*bitarray1++ | *bitarray2++);
	}

	return ((1 + (bitleng >> 6)) << 6) - bit1s;
}

#if SSE2 //vc++ >= 2003
static int countBit0ArrayOrSSE2(const uchar bitarray1[], const uchar* bitarray2, const int bitleng)
{
	const __m128i* pma128 = (__m128i*)bitarray1;
	const __m128i* puma128 = (__m128i*)bitarray2;

	int bit1s = 0;
	int loops = bitleng / 128;

#if _MSC_VER
	__m128i xmm1;
#else
	union
	{
		__m128i m128i;
		uint64 m128i_u64[2];
		uint   m128i_u32[4];
	} xmm1;
#endif
//	assert((uint)pma128 % 16 == 0);

	while (loops-- >= 0) {
		xmm1
#if __GNUC__
		.m128i
#endif
		= _mm_or_si128(_mm_loadu_si128(puma128++), _mm_loadu_si128(pma128++));

#if X86_64
		bit1s +=
			_mm_popcnt_u64(xmm1.m128i_u64[0]) +
			_mm_popcnt_u64(xmm1.m128i_u64[1]);
#else
		bit1s +=
			_mm_popcnt_u32(xmm1.m128i_u32[0]) +
			_mm_popcnt_u32(xmm1.m128i_u32[1]) +
			_mm_popcnt_u32(xmm1.m128i_u32[2]) +
			_mm_popcnt_u32(xmm1.m128i_u32[3]);
#endif
	}

	return (bitleng / 128 + 1) * 128 - bit1s;
}

#elif AVX2

static int countBit0ArrayOrAvx2(const uchar* bitarray1, const uchar* bitarray2, const int bitleng)
{
	int loops = bitleng / 256;
	int bit1s = 0;

	const __m256i *pd1 = (__m256i*) bitarray1;
	const __m256i *pd2 = (__m256i*) bitarray2;

#if _MSC_VER
	__m256i avx2;
#else
	union {
		__m256i m256i;
		uint64 m256i_u64[4];
		uint   m256i_u32[8];
	} avx2;
#endif

	while (loops-- >= 0) {
		avx2
#if __GNUC__
		.m256i
#endif
//		= _mm256_or_si256(_mm256_loadu_si256(pd1++), _mm256_loadu_si256(pd2++));
		= _mm256_or_si256(*pd1++, *pd2++);

#if X86_64
		bit1s +=
			_mm_popcnt_u64(avx2.m256i_u64[0]) +
			_mm_popcnt_u64(avx2.m256i_u64[1]) +
			_mm_popcnt_u64(avx2.m256i_u64[2]) +
			_mm_popcnt_u64(avx2.m256i_u64[3]);
#else
		bit1s +=
			_mm_popcnt_u32(avx2.m256i_u32[0]) + _mm_popcnt_u32(avx2.m256i_u32[1]) +
			_mm_popcnt_u32(avx2.m256i_u32[2]) + _mm_popcnt_u32(avx2.m256i_u32[3]) +
			_mm_popcnt_u32(avx2.m256i_u32[4]) + _mm_popcnt_u32(avx2.m256i_u32[5]) +
			_mm_popcnt_u32(avx2.m256i_u32[6]) + _mm_popcnt_u32(avx2.m256i_u32[7]);
#endif
	}

	return (bitleng / 256 + 1) * 256 - bit1s;
}
#endif

inline static int countBit0ArrayOr(const uchar bitarray1[], const uchar bitarray2[], const int bitleng)
{
#if AVX2
	return countBit0ArrayOrAvx2(bitarray1, bitarray2, bitleng);
#elif SSE2
	return countBit0ArrayOrSSE2(bitarray1, bitarray2, bitleng);
#elif POPCNT
	return countBit0ArrayOrPopcnt((uint64*)bitarray1, (uint64*)bitarray2, bitleng);
#elif X86_64
	return countBit0ArrayOr64((uint64*)bitarray1, (uint64*)bitarray2, bitleng);
#else
	return countBit0ArrayOrPopcnt((uint64*)bitarray1, (uint64*)bitarray2, bitleng);
//	return countBit0ArrayOr32((uint*)bitarray1, (uint*)bitarray2, bitleng);
#endif
}

inline static void countBit0ArrayOr2(const uint64 bitarray1[], const uint64 bitarray2[], const int bitleng, uint64 gp[])
{
	int bits[8] = {0};
	int loops = bitleng >> 6;

	while (loops-- >= 0) {
		const uint64 beg1 = *bitarray1++;
		const uint64 beg2 = *bitarray2++;
		const uint64 next = *bitarray2;
		bits[0] += countBit1(beg1 | (beg2 >> 0));
		bits[1] += countBit1(beg1 | (beg2 >> 1) | (next << 63));
		bits[2] += countBit1(beg1 | (beg2 >> 2) | (next << 62));
		bits[3] += countBit1(beg1 | (beg2 >> 3) | (next << 61));
		bits[4] += countBit1(beg1 | (beg2 >> 4) | (next << 60));
		bits[5] += countBit1(beg1 | (beg2 >> 5) | (next << 59));
		bits[6] += countBit1(beg1 | (beg2 >> 6) | (next << 58));
		bits[7] += countBit1(beg1 | (beg2 >> 7) | (next << 57));
	}

	for (int i = 0; i < 8; i++) {
		gp[0 - i * WHEEL / 2] += ((1 + (bitleng >> 6)) << 6) - bits[i];
	}
}

//the remaining bit of the word where the lastbitpos
//is in is filled with bit 1
static void packQword(uchar bitarray[], const int lastbitpos)
{
	uint64* memstart = (uint64*)bitarray + (lastbitpos >> 6);
	memstart[0] |= ~(((uint64)1u << (lastbitpos % (1 << 6))) - 1);
	memstart[1] = memstart[2] = memstart[3] = (uint64)(~0);
}

//copy bit from srcarray to dstarray by table Map16To30
//it's quite slow
static int convertByteToWheel(ushort srcarray[], const int wordleng, uchar dstarray[])
{
	int lastbit = 32;
	uint* pdst = (uint*)dstarray;
	*pdst = 0;

	for (int i = 0; i <= wordleng; i++) {
		const uint mask = Map16To30[*srcarray++];
		if (lastbit != 32) {
			*pdst++ |= mask << (32 - lastbit);
			*pdst = mask >> lastbit;
			lastbit += 2;
		} else {
			*pdst |= mask;
			lastbit = 2;
		}
	}

	return wordleng;
}

//move bit from [0, bitleng) to [highpos, highpos + bitleng)
static void shiftBitToHigh(uchar bitarray[], const int bitleng, const int highpos)
{
	if (highpos % CHAR_BIT == 0) {
		memmove(bitarray + highpos / CHAR_BIT, bitarray, bitleng / CHAR_BIT + 32);
	} else {
		uint* plowdword = (uint*)bitarray + (highpos + bitleng) / 32;
		const int bitmove = highpos % 16;
		for (int i = bitleng / 32 + 1; i > 0; i--) {
			const uint lowword = *(uint64*)plowdword << bitmove;
			*plowdword -- = lowword;
		}
	}
}

//move bit from [lowpos, lowpos + bitleng) to [0, bitleng)
static void shiftBitToLow(uchar bitarray[], const int bitleng, const int lowpos)
{
	if (lowpos % CHAR_BIT == 0) {
		memcpy(bitarray, bitarray + lowpos / CHAR_BIT, bitleng / CHAR_BIT + 32);
	} else {
		uint* plowdword = (uint*)bitarray + lowpos / CHAR_BIT;
		const int bitmove = lowpos % 16;
		for (int i = bitleng / 32 + 1; i > 0; i--) {
			const uint newword = *(uint64*)plowdword >> bitmove;
			*plowdword++ = newword;
		}
	}
}

//copy from srcarray with bit in [lowpos, lowpos + bitleng) to dstarray in [0, bitleng]
static void shiftBitToLow2(uchar srcarray[], const int bitleng, const int lowpos, uchar dstarray[])
{
	if (lowpos % CHAR_BIT == 0) {
		memcpy(dstarray, srcarray + lowpos / CHAR_BIT, bitleng / CHAR_BIT + 32);
	} else {
		uint* plowdword = (uint*)(srcarray + lowpos / 16);
		uint* pdstdword = (uint*) dstarray;
		const int bitmove = lowpos % 16;

		for (int i = bitleng / 32 + 1; i > 0; i--) {
			*pdstdword++ = *(uint64*)plowdword++ >> bitmove;
		}
	}
}

static void reverseArray(uchar* pbeg, uchar* pend)
{
	while (pbeg < pend) {
		const uchar tmp = *pbeg;
		*pbeg++ = *pend;
		*pend-- = tmp;
	}
}

//reverse word array bitarray with number of byteleng
static void reverseByteArray(uchar bitarray[], const int byteleng)
{
	uint* ps = (uint*)(bitarray);
	uint* pe = (uint*)(bitarray + byteleng - 4);

	while (ps <= pe) {
		const uint tmp = WordReverse[*ps >> 16] | (WordReverse[*ps & 0xffff] << 16);
		*ps++ = WordReverse[*pe >> 16] | (WordReverse[*pe & 0xffff] << 16);
		*pe-- = tmp;
	}

	uchar* pcs = (uchar*)ps, *pce = ((uchar*)pe) + 3;
	while (pcs <= pce) {
		const uchar tmps = WordReverse[*pcs] >> CHAR_BIT;
		const uchar tmpe = WordReverse[*pce] >> CHAR_BIT;
		*pcs++ = tmpe, *pce-- = tmps;
	}
}

//reverse bitarray with bit in [0, bitleng) with bitleng > 0
//swap the kth and (bitleng - kth - 1) bit value of bitarray
static void reverseBitPack(uchar bitarray[], const int bitleng)
{
	const int bitremains = bitleng % CHAR_BIT;
	const int bytes = bitleng / CHAR_BIT;

#if 1
	reverseByteArray(bitarray, (bitleng + CHAR_BIT - bitremains) / CHAR_BIT);
	shiftBitToLow(bitarray, bitleng + 16, CHAR_BIT - bitremains);
#else
	uint buffer[BLOCK_SIZE / 64 + 9];
	uint* pe = (uint*)(bitarray + bytes);
	uint* ps = buffer;

	while ((uchar*)pe >= bitarray) {
		const uint tmp = *(uint64*)--pe >> bitremains;
		*ps++ = (WordReverse[(ushort)tmp]) << 16u | (WordReverse[tmp >> 16]);
	}
	memcpy(bitarray, buffer, bytes + 9);
#endif
}

static void splitToBitArray(uchar srcarray[], const int bitleng, uchar dstarray[][BUFFER_SIZE8])
{
	union
	{
		uint64 mask;
		ushort word[4];
	} utmp;

#if 1
	for (int i = 0; i < 8; i++) {
		//if ((Config.GpMask & (1 << i)) == 0) {
		//		continue;
		//}
		uchar* prc = dstarray[i];
		for (int j = 0; j < bitleng; j += 8) {
			utmp.mask = (*(uint64*)(srcarray + j) >> i) & (0x0101010101010101);
#if BMI2
			*prc++ = _pext_u64(utmp.mask, 0x0101010101010101);//avx2
#else
			*prc++ =
				  (utmp.word[0] | (utmp.word[0] >> 7)) << 0
				| (utmp.word[1] | (utmp.word[1] >> 7)) << 2
				| (utmp.word[2] | (utmp.word[2] >> 7)) << 4
				| (utmp.word[3] | (utmp.word[3] >> 7)) << 6;
#endif
		}
	}
#else
	for (int j = 0; j * 8 < bitleng; j += 1) {
		uint64 bqword = *(uint64*)(srcarray + j * 8);
		for (int i = 0; i < 8; i++) {
			utmp.mask = bqword & 0x0101010101010101;
#if BMI2
			dstarray[i][j] = _pext_u64(bqword, 0x0101010101010101);//avx2
#else
			dstarray[i][j] =
				(utmp.word[0] | (utmp.word[0] >> 7)) << 0 |
				(utmp.word[1] | (utmp.word[1] >> 7)) << 2 |
				(utmp.word[2] | (utmp.word[2] >> 7)) << 4 |
				(utmp.word[3] | (utmp.word[3] >> 7)) << 6;
#endif
			bqword >>= 1;
		}
	}
#endif
}

//get prime from bit buffer
static int savePrime(const uint64 offset, const int wordleng, const ushort* bitarray, uint* dstarray, int init)
{
	int primes = 0;
	static uint total = 0;
	if (offset == 0) {
		total = init;
	}

	for (int bi = 0; bi <= wordleng; bi++) {
		ushort mask = ~bitarray[bi];
		uint64 prime = offset + bi * WHEEL * 2;
		while (mask > 0) {
			const int pi = PRIME_OFFSET(mask);
			uint64 p = prime + Pattern[pi];
			if (dstarray == NULL)
				printf("%u %llu\n", ++total, p);
			else
				dstarray[primes] = uint(p - offset);
			primes++;
			mask &= mask - 1;
		}
	}

	return primes;
}

//get prime from bit buffer
static int dumpGp(uint64 offset, const int wordlengs, const ushort* bitarray, uint gps)
{
	int primes = 0;
	static uint total = 1;
	total = gps;
	if (offset == 0) {
		total = gps;
	}

	offset += 1;
	for (int bi = 0; bi <= wordlengs; bi++) {
		ushort mask = ~bitarray[bi];
		while (mask > 0) {
			printf("%u %llu\n", total++, offset + PRIME_OFFSET(mask) * 2);
			mask &= mask - 1;
			primes++;
		}
		offset += 32;
	}

	return primes;
}

//init the current segment bitarray and pack sievesize
static int setSieveTpl(const uint64 start, const uint sievesize, uchar bitarray[])
{
	int bitleng = sievesize;

	const int si = (start % BLOCK_SIZE) / WHEEL;
	bitleng += (int)(start % WHEEL);
	bitleng = bitleng / WHEEL * 8 + WheelLeng[bitleng % WHEEL];
	const int bytes = bitleng / CHAR_BIT + 1;

	if (si + bytes < sizeof(SievedTpl)) {
		memcpy(bitarray, SievedTpl + si, bytes);
	} else {
		memcpy(bitarray, SievedTpl + si, sizeof(SievedTpl) - si);
		memcpy(bitarray + sizeof(SievedTpl) - si, SievedTpl, bytes + si - sizeof(SievedTpl));
	}

	if (start < WHEEL) {
		bitarray[0] = 0x1;
	}
	//pack the first byte
	bitarray[0] |= (1 << WheelLeng[start % WHEEL]) - 1;

	//pack the last qword bugs
	packQword(bitarray, bitleng);
	if (bitleng % CHAR_BIT) {
		bitleng += CHAR_BIT - bitleng % CHAR_BIT;
	}

	return bitleng / CHAR_BIT;
}

static int getPackLen(uint64 start, const uint sievesize)
{
	int bitleng = sievesize + (int)(start % WHEEL);
	bitleng = bitleng / WHEEL * 8 + WheelLeng[bitleng % WHEEL];
	const int ei = bitleng % CHAR_BIT;

	if (ei != 0) {
		bitleng += CHAR_BIT - ei;
	}
	return bitleng / CHAR_BIT;
}

//init the current segment bitarray and pack sievesize
static int setSieveTpl2(uint64 start, const uint sievesize, uchar bitarray[][BUFFER_SIZE8])
{
	int bitleng = sievesize + (int)(start % WHEEL);
	bitleng = bitleng / WHEEL * 8 + WheelLeng[bitleng % WHEEL];
	const int ei = bitleng % CHAR_BIT;
	const int si = (start % BLOCK_SIZE) / WHEEL;
	const int byteleng = bitleng / CHAR_BIT;

	for (int i = 0; i < 8; i++) {
		shiftBitToLow2(SievedTpl8[i] + si / CHAR_BIT, byteleng + 9, si % CHAR_BIT, bitarray[i]);
		if (Pattern[i] < start % WHEEL)
			SET_BIT(bitarray[i], 0);
		if (ei > i) {
			packQword(bitarray[i], byteleng + 1);
		} else {
			packQword(bitarray[i], byteleng);
		}
	}

	//
	for (int j = 0; j < 8; j++) {
		if (start <= Pattern[j] && start + sievesize > Pattern[j]) {
			CLR_BIT(bitarray[j], 0);
			SET_BIT(bitarray[0], 0);
		}
	}

	if (ei != 0) {
		bitleng += CHAR_BIT - ei;
	}

	return bitleng / CHAR_BIT;
}

inline static void
crossOutFactor4(uchar* pbeg[], const uchar* pend, const uint mask, const uint step)
{
	union
	{
		uint dmask;
		uchar bmask[4];
	} umask;
	umask.dmask = mask;

	uchar* ps0 = pbeg[0], *ps1 = pbeg[1];
	uchar* ps2 = pbeg[2], *ps3 = pbeg[3];

	while (ps3 <= pend) {
#if 0
		*ps0 |= umask.bmask[0]; ps0 += step;
		*ps1 |= umask.bmask[1]; ps1 += step;
		*ps2 |= umask.bmask[2]; ps2 += step;
		*ps3 |= umask.bmask[3]; ps3 += step;
#else
		*ps0 |= mask >> 0; ps0 += step;
		*ps1 |= mask >> 8; ps1 += step;
		*ps2 |= mask >> 16; ps2 += step;
		*ps3 |= mask >> 24; ps3 += step;
#endif
	}

	if (ps0 <= pend)
		*ps0 |= mask >> 0;
	if (ps1 <= pend)
		*ps1 |= mask >> 8;
	if (ps2 <= pend)
		*ps2 |= mask >> 16;
}

#define CHECK_OR(n) if (ps##n <= pend) *ps##n |= 1 << n
inline static void
crossOutFactor8(uchar* pbeg[], const uchar* pend, const uint step)
{
	uchar* ps0 = pbeg[0], *ps1 = pbeg[1];
	uchar* ps2 = pbeg[2], *ps3 = pbeg[3];

	while (ps3 <= pend) {
		*ps0 |= 1 << 0, ps0 += step;
		*ps1 |= 1 << 1, ps1 += step;
		*ps2 |= 1 << 2, ps2 += step;
		*ps3 |= 1 << 3, ps3 += step;
	}
	//not safe for step > L1Size
	CHECK_OR(0);
	CHECK_OR(1);
	CHECK_OR(2);

	uchar* ps4 = pbeg[4], *ps5 = pbeg[5];
	uchar* ps6 = pbeg[6], *ps7 = pbeg[7];
	while (ps7 <= pend) {
		*ps4 |= 1 << 4, ps4 += step;
		*ps5 |= 1 << 5, ps5 += step;
		*ps6 |= 1 << 6, ps6 += step;
		*ps7 |= 1 << 7, ps7 += step;
	}
	CHECK_OR(4);
	CHECK_OR(5);
	CHECK_OR(6);
}

inline static void
sieveGp(uchar bitarray[], const uchar* pend, const uint p, uint offset, uint multiples)
{
	uint mask0 = 0;
#if 1
	uchar* ps0 = NULL;
	for (int m = 0; m < 8; m++) {
		const uchar wordmask = 1 << WheelIndex[offset % WHEEL];
		if (Config.GpMask & wordmask) {
			uchar* ps1 = bitarray + offset / WHEEL;
			if (mask0 == 0) {
				mask0 = wordmask;
				ps0 = ps1;
			} else {
				crossOutFactor2(ps0, ps1, pend, wordmask << CHAR_BIT | mask0, p);
				mask0 = 0;
			}
		}
		offset += (multiples % 4) * 2 * p; multiples /= 4;
	}
	if (mask0) {
		while (ps0 <= pend) {
			*ps0 |= mask0; ps0 += p;
		}
	}
#else
	uchar* pbeg[4] = {0};
	for (int m = 0, pn = 0; m < 8; m++) {
		const uchar wordmask = 1 << WheelIndex[offset % WHEEL];
		if (Config.GpMask & wordmask) {
			mask0 |= wordmask << (pn * 8);
			pbeg[pn++] = bitarray + offset / WHEEL;
		}
		offset += (multiples % 4) * 2 * p; multiples /= 4;
	}
	if (pbeg[3] == 0) {
		pbeg[3] = pbeg[2];
		mask0 |= (mask0 >> 16) << 24;
	}
	crossOutFactor4(pbeg, pend, mask0, p);
#endif
}

static void
sieveSmall0(uchar bitarray[], const uchar* pend, const uint p, uint offset, uint multiples)
{
	uchar* pbeg[8];
	for (int i = 0; i < 8; i++) {
		pbeg[WheelIndex[offset % WHEEL]] = bitarray + offset / WHEEL;
		offset += (multiples % 4) * 2 * p; multiples /= 4;
	}
	crossOutFactor8(pbeg, pend, p);
}

static void
sieveSmall1(uchar bitarray[], const uchar* pend, const uint p, uint offset, uint multiples)
{
	for (int k = 4; k > 0; k--) {
		uchar* ps0 = bitarray + offset / WHEEL;
		ushort wordmask = WheelMask[offset % WHEEL];
		offset += (multiples % 4) * 2 * p; multiples /= 4;

		uchar* ps1 = bitarray + offset / WHEEL;
		wordmask |= WheelMask[offset % WHEEL] << 8;
		offset += (multiples % 4) * 2 * p; multiples /= 4;
		crossOutFactor2(ps0, ps1, pend, wordmask, p);
	}
}

static void
sieveSmall2(uchar dstarray[][BUFFER_SIZE8], const uint sleng, const uint p, uint offset, uint multiples)
{
	for (int k = 4; k > 0; k--) {
		const int s1 = offset / WHEEL;
		const uchar* ps1 = dstarray[WheelIndex[offset % WHEEL]];
		offset += (multiples % 4) * 2 * p; multiples /= 4;

		const int s2 = offset / WHEEL;
		const uchar* ps2 = dstarray[WheelIndex[offset % WHEEL]];
		offset += (multiples % 4) * 2 * p; multiples /= 4;
		set2BitArray2((utype*)ps1, (utype*)ps2, s1, s2, p, sleng);
	}
}

static void
sieveSmall3(uchar dstarray[][BUFFER_SIZE8], const uint sleng, const uint p, uint offset, uint multiples)
{
	for (int k = 0; k < 8; k++) {
		const uint mapi = WheelIndex[offset % WHEEL];
		if (Config.GpMask & (1 << mapi)) {
			utype* ps = (utype*)dstarray[mapi];
			for (uint s = offset / WHEEL; s < sleng; s += p) {
				SET_BIT(ps, s);
			}
		}
		offset += (multiples % 4) * 2 * p; multiples /= 4;
	}
}

static void doSieve2(uchar dstarray[][BUFFER_SIZE8], const uint64 start, const uint sievesize)
{
	uint sqrtp = isqrt(start + sievesize) + 1;
	if ((start + sievesize) < ((uint64)sqrtp) * sqrtp) {
		sqrtp = isqrt(start + sievesize) + 2;
	}

	const int sleng = sievesize / WHEEL * 1 + 8;

	uint j = 8 + FIRST_SIEVE / 23, p = FIRST_SIEVE;
	for (; p < sqrtp; NEXT_PRIME(p, j)) {
		uint offset = p - (uint)(start % p);
		if (start <= p)
			offset = p * p - (uint)start;

		const uint wi = WheelSkip[offset % WHEEL][WheelIndex[p % WHEEL]];
		offset += (wi & 15) * p;
		const uint multiples = WHEEL_SKIP >> ((wi >> 8) * 2);

		if (Config.GpMask != 0xffff) {
			sieveSmall3(dstarray, sleng, p, offset, multiples);
		} else {
			sieveSmall2(dstarray, sleng, p, offset, multiples);
		}
	}
}

//core code of this algorithm
//sieve prime multiples in [start, start + sievesize)
static int segmentedSieve2(const uint64 start, const uint sievesize, uchar dstarray[][BUFFER_SIZE8])
{
	const int byteleng = setSieveTpl2(start, sievesize, dstarray);
	doSieve2(dstarray, start - start % WHEEL, sievesize);

	return byteleng;
}

static void doSieve(uchar bitarray[], const uint64 start, const uint sievesize, const uint minp, uint sqrtp)
{
	if ((start + sievesize) < ((uint64)sqrtp) * sqrtp) {
		sqrtp = isqrt(start + sievesize) + 2;
	}

	const uint bestp = sievesize / WHEEL / 2;
	const uchar* pend = bitarray + sievesize / WHEEL;

	uint j = 8 + FIRST_SIEVE / 23, p = FIRST_SIEVE;
	//TODO:
	while (p < minp) {
		NEXT_PRIME(p, j);
	}

	for (; p < sqrtp; NEXT_PRIME(p, j)) {
		//the min number which satisfy : (start + offset) % p == 0A
		uint offset = p - (uint)(start % p);
		if (start <= p)
			offset = p * p - (uint)start;

		const uint wi = WheelSkip[offset % WHEEL][WheelIndex[p % WHEEL]];
		offset += (wi & 15) * p;
		int sgi = (wi >> 8) * 2;
		uint multiples = (WHEEL_SKIP >> sgi) | (WHEEL_SKIP << (32 - sgi));

		if (p > bestp) {
			while (offset <= sievesize) {
				bitarray[offset / WHEEL] |= WheelMask[offset % WHEEL];
				offset += (multiples % 4) * 2 * p;
				multiples = (multiples >> 2) | (multiples << 30);
			}
		} else {
			//only fast on intel core ix cpu
			//if (0 && WordNumBit1[Config.GpMask] < 5)
			//sieveGp(bitarray, pend, p, offset, multiples);
			//else
			sieveSmall1(bitarray, pend, p, offset, multiples);
		}
	}
}

//core code of this algorithm
//sieve prime multiples in [start, start + sievesize)
static int segmentedSieve(uint64 start, uint sievesize, uchar dstarray[])
{
	uchar bitarray[2 * BUFFER_SIZE + 32];
	const uint maxp = isqrt(start + sievesize) + 1;
	const uint byteleng = setSieveTpl(start, sievesize, bitarray);
	sievesize += start % WHEEL; start -= start % WHEEL;

	for (uint sieveindex = 0, segsize = CpuCache.L1Size; sieveindex < sievesize; sieveindex += segsize) {
		if (segsize + sieveindex > sievesize)
			segsize = sievesize - sieveindex;
		doSieve(bitarray + sieveindex / WHEEL, start + sieveindex, segsize, FIRST_SIEVE, CpuCache.L1Maxp);
	}
	doSieve(bitarray, start, sievesize, CpuCache.L1Maxp, maxp);

	int primes = byteleng;
	const int cmd = (dstarray == NULL) ? 0 : dstarray[0];
	if (cmd == 0) {
		primes = countBit0Array((uint64*)bitarray, byteleng * CHAR_BIT);
	} else if (cmd == CON_VTOWHEEL) {
		primes = convertByteToWheel((ushort*)bitarray, byteleng / 2 + 2, dstarray);
	} else if (cmd == COPY_BYBIT) {
		memcpy(dstarray, bitarray, byteleng + 16);
	} else if (cmd == SAVE_PRIMEDIF) {
		primes = savePrime(start, byteleng / 2, (ushort*)bitarray, (uint*)dstarray, 0);
	}

	return primes;
}

static int setSieveSize(uint sievesize)
{
	if (sievesize < 6000 && sievesize > 15) {
		sievesize *= (WHEEL << 10);
	} else if (sievesize < 16 * (WHEEL << 10)) {
		sievesize = L1_CACHE_SIZE * (WHEEL << 10);
	}

	if (sievesize > BLOCK_SIZE - 2 * MAX_GPCOUNT) {
		sievesize = BLOCK_SIZE - 3 * MAX_GPCOUNT;
	} else if (sievesize < 16 * (WHEEL << 10)) {
		sievesize = L1_CACHE_SIZE * (WHEEL << 10);
	}

	sievesize -= sievesize % (32 * WHEEL);

	Config.SieveSize = sievesize;
	return sievesize;
}

//the sieve of Eratosthenes implementation by bit packing
//all prime less than 2^16 will be saved in prime buffer List
//Prime[0] is the first sieve prime, Prime[i] is the difference
//of the adjacent prime, Prime[i] = Prime[i] - Prime[i - 1];
static void eratoSieve(const uint maxp)
{
	int primes = 1;

	utype bitarray[1020100 >> (BSHIFT + 1)];
	memset(bitarray, 0, (100 + maxp) >> 4);

#if PRIME_DIFF
	Prime[0] = BLOCK_SIZE % 19 == 0 ? 23 : 29;
	Prime[1] = 2;
	uint lastprime = 2;
#endif

	for (uint p = 3; p <= maxp; p += 2) {
		if (!TST_BIT(bitarray, p / 2)) {
#if PRIME_DIFF
			Prime[++primes] = p - lastprime;
			lastprime = p;
#else
			Prime[++primes] = p;
#endif

			if (p > 30000)
				continue;
			for (uint j = p * p / 2; j <= maxp / 2; j += p)
				SET_BIT(bitarray, j);
		}
	}

#if PRIME_DIFF
	//pack the last two byte for safety
	Prime[primes + 2] = Prime[primes + 1] = 255;
#else
	Prime[primes + 1] = 1 << 29;
#endif

	Config.Maxp = maxp;
}

//The first presieved template
//sieve the first 7th prime multiples
static void initSieveTpl( )
{
	memset(SievedTpl, 0, sizeof(SievedTpl));
	const int sievetab[ ] = {7, 11, 13, 17, 19, 23, 29};

	for (int i = 0; BLOCK_SIZE % sievetab[i] == 0; i++) {
		int start = sievetab[i], p2 = 2 * sievetab[i];
		for (; start < FACTP; start += p2) {
			int wi = WheelIndex[start % WHEEL];
			if (wi >= 0)
				SievedTpl[start / WHEEL] |= (1 << (wi & 7));
		}
	}

	for (int j = 1; j < BLOCKS; j++)
		memcpy(SievedTpl + j * FACTP / WHEEL, SievedTpl, FACTP / WHEEL);
	for (int k = 0; k < 8; k++)
		memset(SievedTpl8[k], ~0, sizeof(SievedTpl8[k]));

	const int ssize = sizeof(SievedTpl);
	for (int m = 0; m < 2 * ssize; m++) {
		uchar mask = 0;
		if (m < ssize)
			mask = ~SievedTpl[m];
		else
			mask = ~SievedTpl[m - ssize];
		for (int pi = 0; mask > 0; pi++) {
			if (mask & 1)
				SievedTpl8[pi][m >> 3] &= ~(1 << (m & 7));
			mask >>= 1;
		}
	}
}

static void initWheelSkip()
{
	const uint nextMultiple[8] =
	{
		0x74561230, 0x67325401,
		0x51734062, 0x42170653,
		0x35607124, 0x26043715,
		0x10452376, 0x03216547
	};

	const uchar primewheel[ ] =
	{
		1, 7, 11, 13, 17, 19, 23, 29
	};

	for (int i = 0; i < WHEEL; i += 1) {
		for (int k = 0; k < 8; k++) {
			int step = 0, offset = i;
			if (i % 2 == 0) {
				step += 1; offset += primewheel[k];
			}

			int wi = WheelIndex[offset % WHEEL];
			while (wi < 0) {
				offset += primewheel[k] * 2;
				step += 2;
				wi = WheelIndex[offset % WHEEL];
			}
			WheelSkip[i][k] = step | (wi << 4) | ((nextMultiple[wi] >> (k * 4)) & 15) << 8;
		}
	}
}

//init 4 bit table: WordNumBit1, WordReverse, LeftMostBitIndex, Map16To30
static void initBitTable( )
{
	//1. init WordNumBit1 table in 0-2^16
	int nbitsize = 1 << 16;
	int i;
#if POPCNT == 0
	WordNumBit1[0] = 0;
	for (i = 1; i < nbitsize; i++)
		WordNumBit1[i] = WordNumBit1[i >> 1] + (i & 1);
#endif

	//2. init bit reverseByte of word table
	uchar bytereverse[256] = {0};
	nbitsize = sizeof(WordReverse) / sizeof(WordReverse[0]);
	for (i = 1; i < (1 << 8); i++)
		bytereverse[i] = reverseByte((uchar)i);
	for (i = 1; i < nbitsize; i++)
		WordReverse[i] = bytereverse[i >> 8] | (bytereverse[i & 255] << 8);

	//3. init Map16To30 table
	nbitsize = sizeof(Map16To30) / sizeof(Map16To30[0]);
	for (i = 0; i < nbitsize; i++) {
		uint mask = (uint)(~0 - (0x11 << 30));
		for (int j = 0; j < 16; j++) {
			if ((i & (1 << j)) == 0)
				mask &= ~(1 << (Pattern[j] / 2));
		}
		Map16To30[i] = mask;
	}

#if BIT_SCANF == 0
	//4. init LeftMostBit1 table
	for (int m = 2; m < (1 << 16); m += 2) {
		LeftMostBit1[m + 0] = LeftMostBit1[m >> 1] + 1;
		LeftMostBit1[m + 1] = 0;
	}
#endif
}

//init Prime, SievedTpl and WordNumBit1 table
static void initGp( )
{
	eratoSieve(100130);
	getSystemInfo();
	getCpuInfo();
	initSieveTpl( );
	initBitTable( );
	initWheelSkip( );
	setSieveSize(L1_CACHE_SIZE * 1);
}

static int startTestGp(const int checkloops, int maxn)
{
	if (!freopen("gp.txt", "rb", stdin)) {
		puts("can not read test data file");
		freopen(CONSOLE, "r", stdin);
		return -1;
	}

	const char* gpformat1 = GPFORMAT;
#if _WIN32
	const char* gpformat2 = "%I64d:%d:%I64d\n";
#else
	const char* gpformat2 = "%llu:%d:%llu\n";
#endif

	int maxgpcout = 0;
	uint64 start = 1 << 30;
	uint64 filedata[MAX_GPCOUNT + 32];
	if (scanf(gpformat2, &start, &maxgpcout, &filedata[0]) != 3 || filedata[0] != 2) {
		puts("wrong gp data type start:maxgpcout:step in the first line");
		freopen(CONSOLE, "r", stdin);
		return -2;
	}
	if (maxgpcout > MAX_GPCOUNT)
		maxgpcout = MAX_GPCOUNT;
	for (int j = 0, failred = 0; j < maxgpcout && failred < 10; j++) {
		uint64 tmp;
		if (scanf(gpformat1, &tmp, filedata + j) != 2) {
			printf("wrong gp data format %s\n", gpformat1);
			maxgpcout = j--;
			failred++;
		}
	}

	if (maxgpcout < maxn && maxgpcout > 0)
		maxn = maxgpcout;

	Config.PrintRet = Config.PrintTime = false;
	Config.PrintGap = 0;

	printf("Test gp data %llu:%d, input cases %d, maxn %d\n", start, maxgpcout, checkloops, maxn);

	uint64 minn = start;
	srand(time(0));
	for (int i = 0; i < checkloops; i++) {
		int itemfails = 0;
		int gcount2 = rand() % maxn + 1;
		if (i + gcount2 > maxgpcout) {
			minn += 2;
			continue;
		}

		setSieveSize(rand() % 64 + 32);
		Config.Algorithm = rand() % 3 + 1;
		Config.Threads = rand() % 4 + 1;
		Config.Advanced = !Config.Advanced;

		getGp2(minn, gcount2, GP);
		for (int j = 0; j < gcount2; j++) {
			if (filedata[j + i] != GP[j]) {
				itemfails++;
				printf("case %d wrong data G(%llu) = %llu != %llu (file)\n", i + 1, minn + j * 2, GP[j], filedata[i + j]);
			}
		}

		printf("loop %d : [%llu, %d] ", i + 1, minn, gcount2);

		if (itemfails > 0)
			printf("%d cases fail; executeCmd( c%d a%d t%d)\n",
				itemfails, Config.SieveSize, Config.Algorithm, Config.Threads);
		else
			printf(" g(%d) t(%d) pass\r", Config.Algorithm, Config.Threads);
		minn += 2;
	}
	putchar('\n');

	//restore the default configuration
	Config.PrintRet = Config.PrintTime = true;
	Config.PrintGap = 0;

	freopen(CONSOLE, "r", stdin);

	return 0;
}

//split command to cmdparams array
static int splitCmdParams(const char* ccmd, char cmdparams[][40])
{
	int nwords = 0;

	for (int i = 0; ; i++) {
		while (isspace(*ccmd)) {
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

static int addOneGp35(const uint64 maxn)
{
	int gps = 2, gmask = 0;
	for (uint j = 2, p = 3; p <= maxn / p && gps > 0; NEXT_PRIME(p, j)) {
		const uint rid = maxn % p;
		if (rid == 0 && !(gmask & 0x1)) {
			gmask |= 0x1;
			gps--;
		} else if (rid == 2 && !(gmask & 0x2)) {
			gmask |= 0x2;
			gps--;
		}
	}
	if (maxn == 3 || maxn == 5) {
		gps = 1;
		gmask |= 0x2;
	}

	if (gps > 0 && Config.PrintGppair) {
		if (gps == 2)
			puts("1 3\n2 5");
		else if (gmask & 0x02)
			puts("1 3");
		else
			puts("1 5");
	}

	return gps;
}

//count goldbach partition [minn, minn + 2 * gpcount - 2]
//with one of the smallest prime is 3 or 5
static void addMultiGp35(const uint64 minn, const int gpcount, uint64 gp[])
{
	utype bitarray[(MAX_GPCOUNT >> BSHIFT) + 64];

	const int64 start = minn - 6;
	const int mod30bits = start % WHEEL / 2;

	bitarray[0] = CON_VTOWHEEL;
	int bytes = segmentedSieve(start, gpcount * 2 + 2, (uchar*)bitarray);

	//5, minn - 5
	if (!TST_BIT(bitarray, mod30bits) && start + 1 >= 5)
		gp[0] ++;
	for (int i = 1; i <= gpcount; i++) {
		if (!TST_BIT(bitarray, (i + mod30bits))) {
			const uint64 x = start + (i + mod30bits) * 2 + 1;
			if (x >= 5)
				gp[i - 0] ++; // 5 + x
			if (x >= 3)
				gp[i - 1] ++; // 3 + x
		}
	}
}

static void addGpInMiddle(const uint64 minn, const int gpcount, uint64 gp[])
{
	const uint64 half = minn / 2;
	const uint64 start = half + half % 2;
	const int offset = (start - start % WHEEL) - half;
	uint prime[MAX_GPCOUNT / 6];

	prime[0] = SAVE_PRIMEDIF;

	int pn = segmentedSieve(start, gpcount * 2 - 1, (uchar*)prime);

	for (int i = 0; i < pn; i++) {
		for (int j = i; j < pn; j++) {
			int gpi = (prime[i] + prime[j]) / 2 + offset;
			if (gpi < gpcount) {
				gp[gpi] ++;
			} else {
				break;
			}
		}
	}
}

/*** ------- (start)---> (start + sievesize) ---(start2)----->----(start2 + sievesize)-----
                   \      /                  /|\                  /|\
                    \    /                    |                    |
                     \  /                     |                    |
                      \/      -----------\ _  |                    |
                      /\      -----------/    |                    |
                     /  \                     |                    |
                    /    \                    |                    |
                   /      \                   |                    |
                  /        \                 \|/                  \|/
// ----(n-start-leng)->--(n - start) ---(n - start2)---<---(n - start2 - sievesize)----
**/
//start + 1, start + sievesize + 1
static int segmentedGp01(const uint64 start1, const uint64 start2, const uint sievesize)
{
	//assert(start + sievesize <= start2 + 1 && 0 == start % WHEEL);
	uchar bitarray1[BLOCK_SIZE / (CHAR_BIT * 2)] MEM_ALIGN(32);
	uchar bitarray2[BLOCK_SIZE / (CHAR_BIT * 2)] MEM_ALIGN(32);
	int gps = 0;

	bitarray1[0] = bitarray2[0] = CON_VTOWHEEL;
	segmentedSieve(start1, sievesize, bitarray1);
	segmentedSieve(start2, sievesize, bitarray2);
	//10% time use
	reverseBitPack(bitarray2, (sievesize + start2 % WHEEL) >> 1);

	//add smallest prime 3 and 5
	if (start1 == 0 && (bitarray2[0] & 0x6) != 0x6) {
		gps = addOneGp35(start2 + sievesize - 4);
	}

	for (int i = 0; i * (CHAR_BIT * 2) <= sievesize; i += 8)
		*(uint64*)(bitarray1 + i) |= *(uint64*)(bitarray2 + i);

	if (Config.PrintGppair) {
		gps += dumpGp(start1, sievesize >> 5, (ushort*)bitarray1, gps + 1);
	} else {
		gps += countBit0Array((uint64*)bitarray1, sievesize >> 1);
	}

	return gps;
}

//more fast than other
static int segmentedGp02(const uint64 start1, const uint64 start2, const uint sievesize)
{
	//assert(start + sievesize <= start2 + 1 && 0 == start % WHEEL);
	uchar bitarray1[BUFFER_SIZE + 32] MEM_ALIGN(32);
	uchar bitarray2[BUFFER_SIZE + 32] MEM_ALIGN(32);
	int gps = 0;

	bitarray1[0] = bitarray2[0] = COPY_BYBIT;
	int byteleng1 = segmentedSieve(start1, sievesize, bitarray1);
	int byteleng2 = segmentedSieve(start2, sievesize, bitarray2);

	//12% time, convert prime wheel pattern to gp pattern
	for (int i = byteleng2 - 1, j = 0; i >= 0; j += 2) {
		const uint gmask = *(uint*)(bitarray2 + (i -= 2));
		*(ushort*)(bitarray1 + j) |=
			(PatternMask[(ushort)(gmask >> CHAR_BIT)]) |
			(PatternMask[(ushort)gmask] << CHAR_BIT);
	}

	if (start1 == 0) {
		//add smallest prime 3 and 5
		gps = addOneGp35(start2 + sievesize - 4);
	}

	if (Config.PrintGppair) {
		gps += savePrime(start1, byteleng1 / 2, (ushort*)bitarray1, NULL, gps);
	} else {
		gps += countBit0Array((uint64*)bitarray1, byteleng1 * CHAR_BIT);
	}

	return gps;
}

//the first algorithm to segmented goldbach partition
static void segmentedGp1(const uint64 start1, const uint64 start2, const uint sievesize, const int gpcount, uint64 gp[])
{
	//	assert(start + sievesize <= start2 + 1 && 0 == start % WHEEL);
	uchar bitarray1[BLOCK_SIZE / (CHAR_BIT * 2)] MEM_ALIGN(32);
	uchar bitarray2[BLOCK_SIZE / (CHAR_BIT * 2)] MEM_ALIGN(32);

	bitarray1[0] = bitarray2[0] = CON_VTOWHEEL;
	//3, 5 pos is set bit 1
	segmentedSieve(start1, sievesize, bitarray1);
	segmentedSieve(start2, sievesize + gpcount * 2 + WHEEL, bitarray2);

	int bitleng = (sievesize + start2 % WHEEL + gpcount * 2 - 2) >> 1;
	reverseBitPack(bitarray2, bitleng);

	if (gpcount == 1) {
		gp[0] += countBit0ArrayOrPopcnt((uint64*)bitarray1, (uint64*)(bitarray2), sievesize >> 1);
		return;
	}

	for (int movei = 1; movei <= (1 << 3); movei++) {
		int gpi = gpcount - movei;
		if (gpi < 0)
			break;
		for (int i = 0; gpi >= 0; i++) {
			const int bit0s = countBit0ArrayOr(bitarray1, bitarray2 + i, sievesize >> 1);
#if OMP
			#pragma omp atmoic
#endif
			{
				gp[gpi] += bit0s;
			}
			gpi -= 1 << 3;
		}
		if (movei < (1 << 3)) {
			shiftBitToLow(bitarray2, bitleng, 1);
		}
	}
}

//the seconds algorithm to segmented goldbach partition
//10% fast than the third algorithm for gpcount < 1000
static void segmentedGp2(const uint64 start1, const uint64 start2, const uint sievesize, const int gpcount, uint64 gp[])
{
	//assert(start + sievesize <= start2 + 1 && 0 == start % WHEEL);
	uchar tmparray[BUFFER_SIZE + 32] MEM_ALIGN(32);
	uchar bitarray1[8][BUFFER_SIZE8] MEM_ALIGN(32);
	uchar bitarray2[8][BUFFER_SIZE8] MEM_ALIGN(32);

	//sieve the first bit array
	tmparray[0] = COPY_BYBIT;
	const int begbit = segmentedSieve(start1, sievesize, tmparray);

	//time use 20%
	//assert(((size_t)&bitarray1[0][0]) % 16 == 0);
	//assert(((uint)&bitarray2[0][0]) % 16 == 0);
	splitToBitArray(tmparray, begbit, bitarray1);
	for (int k = 0; k < 8; k++)
		packQword(bitarray1[k], begbit);

	//sieve the seconds bit array
	const int leng2 = sievesize + (gpcount - 1) * 2 - 1;
	tmparray[0] = COPY_BYBIT;
	const int maxbit = segmentedSieve(start2, leng2, tmparray);

	reverseArray(tmparray, tmparray + maxbit - 1);
	splitToBitArray(tmparray, maxbit, bitarray2);

	int offset = (start1 + start2 + leng2) % WHEEL;
	if (0 == offset)
		offset = WHEEL;
	offset = 2 * gpcount - 2 - offset;

	for (int lwsi = 0; lwsi < 512; lwsi++) {
		const int li = lwsi >> 6, wi = (lwsi >> 3) % 8, si = lwsi % 8;
		int j = 0;
		int gpi = (Pattern[si] + Pattern[wi] - li * WHEEL + offset) / 2;
		if (gpi >= gpcount) {
			gpi -= 4 * WHEEL;
			j = 1;
		}
		for ( ; gpi >= 0; j++) {
			const int bit0s = countBit0ArrayOr(bitarray1[si], bitarray2[wi] + j, begbit);
#if OMP
			#pragma omp atmoic
#endif
			{
				gp[gpi] += bit0s;
			}
			gpi -= 4 * WHEEL;
		}

		if (si == 7 && li < 7)
			shiftBitToLow(bitarray2[wi], maxbit, 1);
	}
}

//the third algorithm to segmented goldbach partition
static void segmentedGp3(const uint64 start1, const uint64 start2, const uint sievesize, const int gpcount, uint64 gp[])
{
	//assert(start + sievesize <= start2 + 1);// & BUFFER_SIZE > sievesize);
	uchar bitarray1[8][BUFFER_SIZE8] MEM_ALIGN(32);
	uchar bitarray2[8][BUFFER_SIZE8] MEM_ALIGN(32);

	//assert(((size_t)&bitarray1[0]) % 16 == 0);
	//assert(((uint)&bitarray2[0]) % 16 == 0);

	const int begbit = segmentedSieve2(start1, sievesize, bitarray1);
	const int leng2 = sievesize - 1 + (gpcount - 1) * 2;
	int packlen = 0;

#if 0
	for (int j = 0; j < 240 && start2 > j; j += 8) {
		const int pack = getPackLen(start2 - j, leng2 + j);
		if (pack % CHAR_BIT == 0) {
			packlen = j;
			break;
		}
	}
#endif

	const int maxbit = segmentedSieve2(start2 - packlen, leng2 + packlen, bitarray2);
	for (int i = 0; i < 8; i++) {
		reverseBitPack(bitarray2[i], maxbit);
	}

	int offset = (start1 + start2 + leng2) % WHEEL;
	if (0 == offset)
		offset = WHEEL;
	offset = 2 * gpcount - 2 - offset;

	for (int lwsi = 0; lwsi < 512; lwsi++) {
		const int li = lwsi >> 6, wi = (lwsi >> 3) % 8, si = lwsi % 8;
		int gpi = (Pattern[si] + Pattern[wi] - li * WHEEL + offset) / 2;
		int j = 0;
		if (gpi >= gpcount) {
			gpi -= 4 * WHEEL;
			j = 1;
		}
		for ( ; gpi >= 0; j++) {
			const int bit0s = countBit0ArrayOr(bitarray1[si], bitarray2[wi] + j, begbit);
#if OMP
			#pragma omp atmoic
#endif
			{
				gp[gpi] += bit0s;
			}
			gpi -= 4 * WHEEL;
		}

		if (si == 7 && li < 7)
			shiftBitToLow(bitarray2[wi], maxbit, 1);
	}
}

static void initPatternMask(const uint wheelmod)
{
	Config.GpMask = 0;

	uchar gppair[16][2], pairs = 0;
	for (int i = 0; i < 8; i++) {
		for (int j = 0; j < 8; j++) {
			const uint psum = Pattern[i] + Pattern[j];
			if (psum % WHEEL == wheelmod) {
				gppair[pairs][0] = i;
				gppair[pairs][1] = j + 8;
				if (psum > WHEEL)
					gppair[pairs][1] = j;
				Config.GpMask |= (1 << i) | (1 << j);
				pairs++;
				break;
			}
		}
	}

	//assert(Config.GpMask == (uchar)GpMask[wheelmod / 2]);

	memset(PatternMask, (uchar)(-1), sizeof(PatternMask));
	for (int k = 0; k < (1 << 16); k++) {
		for (int j = 0; j < pairs; j++) {
			if (!(k & (1 << gppair[j][1]))) {
				//set bit 0
				PatternMask[k] &= ~(1 << gppair[j][0]);
			}
		}
	}
}

/*** ------- (start) ----- (start+sievesize) -----(start)----->---(start + sievesize)------
                 \        /                     /|\                /|\
                  \      /                       |                  |
                   \    /                        |                  |
                    \  /                         |                  |
                     \/        ---------\        |                  |
                     /\        ---------/        |                  |
                    /  \                         |                  |
                   /    \                        |                  |
                  /      \                       |                  |
                 /        \                     \|/                \|/
--(n-start-leng)--(n-start+gpcount*2) --(n-start+gpcount*2)--<---(n-start-leng)-----
<------------------------------------------------------------------------------
loop move gpcount bit
**/
static void getGp1(const uint64 minn, const uint gpcount, const uint gstep, uint64 gp[])
{
	if (Config.SaveResult) {
		Config.PrintGap = 0;
		Config.PrintRet = true;
		Config.PrintTime = false;
		freopen("gp.txt", "wb+", stdout);
		printf("%llu:%u:%d\n", minn, gpcount, 2);
	}

	assert(minn >= 6 && minn % 2 == gstep % 2);

	if (minn / Config.Maxp > Config.Maxp) {
		eratoSieve(isqrt(minn) + MAX_GPCOUNT);
	}

	uint64 begin = minn;

	for (int i = 0; i < gpcount; i++) {

		const double ts = getTime();
		const uint64 half = begin / 2;
		uint64 sgp = 0;
		int step = Config.SieveSize;
		initPatternMask(begin % WHEEL);

		if (step > half)
			step = (int)half;
#if (OMP)
		omp_set_num_threads(Config.Threads);
		#pragma omp parallel for reduction(+:sgp) \
		schedule(dynamic) if (minn > 10000000)
#endif
		for (uint64 start = 0; start < half; start += step) {
			if (start > half - step)
				step = half - start;
			if ((start & Config.PrintGap) == 0 && Config.PrintGap > 0)
				printf("finish ... %02d%%\r", (int)(100 * start / half));
			if (Config.Algorithm > 1)
				sgp += segmentedGp02(start, begin - start - step, step + 1);
			else
				sgp += segmentedGp01(start, begin - start - step, step + 1);
		}

		if (Config.PrintRet) {
			printf("G(%lld) = %u", begin, (uint)sgp);
			if (Config.PrintTime && begin > 10000000)
				printf(", time use %.3f sec", (getTime() - ts) / 1000);
			putchar('\n');
		}

		if (gp && i < MAX_GPCOUNT)
			gp[i] = sgp;
		begin += gstep;
	}

	if (Config.SaveResult) {
		Config.PrintGap = (1 << 16) - 1;
		Config.SaveResult = false;
		freopen(CONSOLE, "w", stdout);
	}
	Config.GpMask = 0xffff;
}

//thread proc
static void coreSieve2(uint64 minn, const int gpcount, const uint64 first, uint64 gp[])
{
	const uint64 half = minn / 2;
	const int gstep = Config.Threads * Config.SieveSize;
	uint sievesize = Config.SieveSize;

	double ts = getTime();
	if (minn > Maxn + 2 * MAX_GPCOUNT)
		minn = Maxn;

	memset(gp, 0, sizeof(gp[0]) * gpcount);
	sieve_func func[] = {segmentedGp1, segmentedGp2, segmentedGp3};
	sieve_func segmentedGp = func[Config.Algorithm - 1];

	for (uint64 start = first; start < half; start += gstep) {
		if (start > half - sievesize || half < sievesize)
			sievesize = uint(half - start);
		if ((start & Config.PrintGap) == 0 && Config.PrintGap > 0) {
			int dots = 100 * start / half;
			printf("\r %.2f sec, %02d%% ", (getTime() - ts) / 10 / dots, dots);
			for (int i = 0; i * 2 < dots; i++)
				putchar('*');
		}
		segmentedGp(start, minn - start - sievesize, sievesize + 1, gpcount, gp);
	}

	if (first == 0) {
		addGpInMiddle(minn, gpcount, gp);
	}
}

static void coreSieve1(uint64 minn, const int gpcount, uint64 gp[])
{
	const uint64 half = minn / 2;
	uint sievesize = Config.SieveSize;
	if (sievesize > half)
		sievesize = half;
	if (minn > Maxn + 2 * MAX_GPCOUNT)
		minn = Maxn;

	double ts = getTime();
	memset(gp, 0, sizeof(gp[0]) * gpcount);

	sieve_func func[] = {segmentedGp1, segmentedGp2, segmentedGp3};
	sieve_func segmentedGp = func[Config.Algorithm - 1];

#if OMP
	omp_set_num_threads(Config.Threads);
	#pragma omp parallel for schedule(dynamic) if (minn > 10000000)
#endif
	for (uint64 start = 0; start < half; start += sievesize) {
		if (start > half - sievesize)
			sievesize = uint(half - start);
		if ((start & Config.PrintGap) == 0 && Config.PrintGap > 0) {
			int dots = (100 * start / half);
			printf("\r %.2f sec, %02d%% ", (getTime() - ts) / 10 / dots, dots);
			for (int i = 0; i * 2 < dots; i++)
				putchar('#');
		}
		segmentedGp(start, minn - start - sievesize, sievesize + 1, gpcount, gp);
	}

	if (gpcount > 1)
		addGpInMiddle(minn, gpcount, gp);
}

//the second algorithm to calculate
//GP [minn, minn + 2 * gpcount - 2]
static void getGp2(uint64 minn, int gpcount, uint64 gp[])
{
	assert(minn % 2 == 0 && gpcount > 0);

	double ts = getTime();
	if (gpcount > MAX_GPCOUNT)
		gpcount = MAX_GPCOUNT;

	if (minn / Config.Maxp > Config.Maxp) {
		eratoSieve(isqrt(minn) + MAX_GPCOUNT);
	}

	Config.GpMask = GpMask[minn % WHEEL / 2];
	for (int i = 1; i < gpcount && Config.GpMask != 0xffff; i++)
		Config.GpMask |= GpMask[(minn + 2 * i) % WHEEL / 2];

	if (Config.SaveResult) {
		Config.PrintGap = 0;
		Config.PrintRet = true;
		freopen("gp.txt", "wb+", stdout);
		printf("%llu:%d:%d\n", minn, gpcount, 2);
	}

	if ((minn > 10000000 || gpcount > 1000) && OMP == 0) {
		if (Config.Threads > MAX_THREADS) {
			Config.Threads = 8;
		}
		if (Config.Threads > gpcount && !Config.Advanced) {
			Config.Threads = 1;
		}
		startWorkThread(Config.Threads, minn, gpcount);
	} else {
		coreSieve1(minn, gpcount, gp);
	}

	Config.GpMask = 0xffff;
	addMultiGp35(minn, gpcount, gp);

	ts = getTime() - ts;

	if (Config.PrintRet) {
		putchar('\n');
		for (int j = 0; j < gpcount; j++) {
			printf(GPFORMAT, minn + j * 2, gp[j]);
		}
	}

	if (Config.PrintTime) {
		printf("\r\ncount = %d, algorithm = %d, threads = %d", gpcount, Config.Algorithm, Config.Threads);
		printf(", time use %.3lf sec\n", ts / 1000);
	}

	if (Config.SaveResult) {
		Config.PrintGap = (1 << 16) - 1;
		freopen(CONSOLE, "w", stdout);
		Config.SaveResult = false;
	}
}

static void printInfo( )
{
	puts("---------------------------------------------------------------");
	puts("Count Multi Goldbach Partition in [6, 1e14], version 2.0\n");
	puts("Copyright @ by Huang Yuanbing 2011 - 2018 bailuzhou@163.com\n"
	"Code:<https://github.com/ktprime/ktprime/blob/master/FastGn.cpp>\n"
	"CXXFLAG:g++ -march=native -mpopcnt -funroll-loops -O3 -s -pipe");

#ifdef _MSC_VER
	printf("Compiled by MS/vc++ %d", _MSC_VER);
#else
	printf("Compiled by GNU/g++ %s", __VERSION__);
#endif

#if X86_64
	printf(" x86-64 ");
#endif
	printf(" on %s %s\n", __TIME__, __DATE__);

	puts("---------------------------------------------------------------");
	puts("---------------------------------------------------------------");
	printf("MACRO: ASM_X86 = %d, DATA_TYPE = %d, POPCNT = %d\n",
		ASM_X86, BSHIFT, POPCNT);

	printf("Work threads = %d, Segment Cache ~= %d kb, Algorithm = %d\n",
		Config.Threads, (Config.SieveSize / WHEEL) >> 10, Config.Algorithm);

	puts("---------------------------------------------------------------\n");
}

static void doCompile()
{
	char exename[256];
	strcpy(exename, __FILE__);
	char* pdot = strchr(exename, '.');
	if (pdot) {
		strcpy(pdot, "_.exe");
	}

	const char* const cxxflag =
#if _MSC_VER
		"cl /O2 /Oi /Ot /Oy /GT /GL %s %s";
#elif X86_64
		"g++ -m64 -mpopcnt -march=native -funroll-loops -O2 -s -pipe %s -o %s";
#else
		"g++ -m32 -mpopcnt -march=native -funroll-loops -O3 -s -pipe %s -o %s";
#endif

	char compileLine[256] = {0};
	sprintf(compileLine, cxxflag, __FILE__, exename);
	puts(compileLine);
	system(compileLine);
}

//set config and return the first cmd index
static int parseConfig(const char cmdparams[][40])
{
	int cmdi = -1;

	for (int i = 0; cmdparams[i][0]; i++) {
		char c = toupper(cmdparams[i][0]);

		if (isdigit(c) || c == 'E') {
			if (cmdi < 0)
				cmdi = i;
			continue;
		}

		int tmp = atoi(cmdparams[i] + 1);
		switch (c)
		{
			case 'S':
				setSieveSize(tmp);
				break;
			case 'T':
				if (tmp < MAX_THREADS && tmp > 0) {
					Config.Threads = tmp;
				}
				break;
			case 'G':
				if (tmp < 4 && tmp > 0)
					Config.Algorithm = tmp;
				break;
			case 'A':
				Config.Advanced = !Config.Advanced;
				break;
			case 'C':
				if (tmp > 128 || tmp < 8) {
					tmp = 63;
				}
				CpuCache.L1Size = (tmp << 10) * WHEEL;
				CpuCache.L1Maxp = (tmp << 10) / L1_SIEVE_SEG;
			case 'M':
				if (tmp >= 0 && tmp <= 30) {
					Config.PrintGap = (1 << tmp) - 1;
				}
				break;
			case 'P':
				c = cmdparams[i][1];
				if (c == 't')
					Config.PrintTime = !Config.PrintTime;
				else if (c == 'r')
					Config.PrintRet = !Config.PrintRet;
				else if (c == 'g')
					Config.PrintGppair = !Config.PrintGppair;
				break;
			case 'F':
				Config.SaveResult = true;
				break;
			case 'H':
				puts(CmdInfo);
				puts(HelpCommand);
				break;
			default:
				cmdi = i;
		}
	}

	return cmdi;
}

//execute command cmd:[H, B, U, R, N, S]
static bool executeCmd(const char* cmd)
{
	while (cmd) {

		//each command split by ';'
		char* pcmd = (char*) strchr(cmd, ';');
		char cmdparams[8][40] = {0};

		if (splitCmdParams(cmd, cmdparams) <= 0) {
			return false;
		}

		int cmdi = parseConfig(cmdparams);
		if (cmdi < 0)
			return true;

		char cmdc = toupper(cmdparams[cmdi][0]);

		if (cmdc == 'B') {
			uint64 minn = atoint64("e9", 0);
			int gpcount = (int)atoint64(cmdparams[cmdi + 1], 10000);
			getGp2(minn, gpcount, GP);
		} else if (cmdc == 'U') {
			int cases = (int)atoint64(cmdparams[cmdi + 1], 1000);
			int gpcount = (int)atoint64(cmdparams[cmdi + 2], 0);
			if (gpcount > 0) {
				startTestGp(cases, gpcount);
			}
		} else if (cmdc == 'R') {
			uint64 minn = atoint64(cmdparams[cmdi + 1], 1000000000);
			int gpcount = (int)atoint64(cmdparams[cmdi + 2], 10);
			int step = (int)atoint64(cmdparams[cmdi + 3], 2);
			getGp1(minn, gpcount, step, 0);
		} else {
			char d = cmdparams[cmdi][0];
			if (isdigit(d) || toupper(d) == 'E') {
				uint64 minn = atoint64(cmdparams[cmdi], 1000000000);
				int gpcount = (int)atoint64(cmdparams[cmdi + 1], 1000);
				getGp2(minn, gpcount, GP);
			}
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
		puts(CmdInfo);
	}

	initGp( );

	for (int i = 1; i < argc; i++) {
		const char c = argv[i][0];
		if (c == 'm') {
			doCompile( );
		} else if (c == 'u') {
			srand(time(0));
			executeCmd("f g1 r 30 e4; u e4 10; u e3 20");
			executeCmd("f g2 r 3e5 e4; u 10000 100; u 1000 1000");
			for (int i = 10; i <= 32; i++) {
				char cmd[60] = {0};
				sprintf(cmd, "f 2^%d 10000; u %d %d;", i, 1000, rand() % 1000 + 10);
				executeCmd(cmd);
			}
			executeCmd("g3 f 1e10 300; u 100 200");
		} else {
			executeCmd(argv[i]);
		}
	}

//	executeCmd("t1 e9 e3 g1; e9 e3 g2; e9 e3 g3");
	executeCmd("t4 e9 e4 g2; e9 e4 g3");

	char ccmd[255] = {0};
	while (true) {
		printf("\n>> ");
		if (!gets(ccmd) || !executeCmd(ccmd))
			break;
	}

	return 0;
}

/******************************************************************
copyright (C) 2009 - 2018 by Huang Yuan bing
mail to: bailuzhou@163.com
free use for non-commercial purposes

G[10^9] to G[10^9 + 2e4]
2.66G Intel I5 560M 8.30 x86
2.26G Intel I3 350M 7.90 x86
2.26G Intel I3 350M 5.68 x64
3.20G Intel I5 3470 3.80 x86
3.20G Intel I5 3470 1.62 x64

G[10^9] to G[10^9 + 2e3]
popcnt64 : 1620
popcnt32 : 2063
tree3    : 2600
table16  : 3560

http://graphics.stanford.edu/~seander/bithacks.html
*******************************************************************

n    pi(n)         g(n)
e08  5761455       291400
e09  50847534      2274205
e10  455052511     18200488
e11  4118054813    149091160
e12  37607912018   1243722370
e13  346065536839  10533150855
e14  3204941750802 90350630388
e15                783538341852
http://code.google.com/p/pcxprj/wiki/CompilerAndLinkerSwitchGuide
-fprofile-generate
-fprofile-use

g++ -O3 -funroll-loops -mpopcnt -s -pipe -march=native -fomit-frame-pointer -lpthread FastGn.cpp -o FastGn.exe

build by vc++
	cl /O2 FastGp.cpp

************************************************************************/

