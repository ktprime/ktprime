#include <algorithm>
#include <numeric>
#include <queue>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <climits>

#define  MAXN  10000*10000
#define  MAXK  100000

using namespace std;

template<class T>
class maxheap
{
public:
	maxheap(int maxk) {
		size = 0;
		a[0] = INT_MAX;
		for (int i = 0; i <= maxk; i++)
			a[maxk + i] = INT_MIN;
	}

	T top() const { return a[1]; }

	void push(const T x)
	{
		int c = ++size;
		int p = size / 2;

		while (x > a[p] /*&& p >= 1*/) {
			a[c] = a[p];
			c = p;
			p /= 2;
		}
		a[c] = x;
	}

	void pop()
	{
		const T x = a[size--];
		int p = 1, c = 1;

		while (x < a[c]/* && c <= size &&*/) {
			a[p] = a[c];
			p = c;
			c *= 2;
			if (a[c + 1] > a[c])
				c++;
		}
		a[p] = x;
	}

//private:

	T a[MAXK * 2 + 2];//1..n
	int size;
};

static int buff[MAXN + 2];

void rand_swap(int a[], int n, int k)
{
	const int step = n / k;
//	std::sort(a, a + k);
	for (int i = 1; i < k; i ++) {
#if RD
		int h = rand() % k, t = i * step + rand() % step;
		if (a[h] > a[t])
			std::swap(a[h], a[t]);
		else if (a[h + 1] > a[t + 1])
			std::swap(a[h + 1], a[t + 1]);
		assert(t < n);
#endif
	}
}

void reset(int a[], int n)
{
	memcpy(a, buff, n * sizeof(a[0]));
	rand_swap(a, n, 10000);
}

#if __linux__
#include <sys/resource.h>
#elif _WIN32
//#include <windows.h>
#include <intrin.h>
#endif

static clock_t getTime()
{
#if 0
	FILETIME ptime[4] = {0};
	GetThreadTimes(GetCurrentThread(), &ptime[0], &ptime[1], &ptime[2], &ptime[3]);
	return (ptime[2].dwLowDateTime + ptime[3].dwLowDateTime) / 10000;
	//return clock();
#elif __linux__
	struct rusage rup;
	getrusage(RUSAGE_SELF, &rup);
	long sec  = rup.ru_utime.tv_sec  + rup.ru_stime.tv_sec;
	long usec = rup.ru_utime.tv_usec + rup.ru_stime.tv_usec;
	return sec * 1000 + usec / 1000;
#else
	return clock();
#endif
}

void stl_sort(int a[], int n, const int k)
{
	reset(a, n);
	auto ts = getTime();

	std::sort(a, a + n);
	printf("stl sort       %4ld ms, a[%d] = %d\n\n", getTime() - ts, k, a[k - 1]);
}

void stl_nth(int a[], int n, const int k)
{
	reset(a, n);
	auto ts = getTime();

	std::nth_element(a, a + k, a + n);
	int maxe = *std::max_element(a, a + k);
	int sum =  accumulate(a, a + k, 0);
	printf("stl nth_element  %4ld ms, a[%d] = %d, sum = %d\n", getTime() - ts, k, maxe, sum);
}

void stl_makeheap(int a[], int n, const int k)
{
	reset(a, n);
	auto ts = getTime();

	std::make_heap(a, a + k);

	int maxe = a[0];
	for (int i = k; i < n; ++i) {
		if (a[i] < maxe) {
			pop_heap(a, a + k);
			a[k - 1] = a[i];
			push_heap(a, a + k);
			maxe = a[0];
		}
	}

	int sum = accumulate(a, a + k, 0);
	printf("stl make_heap   %4ld ms, a[%d] = %d, sum = %d\n", getTime() - ts, k, maxe, sum);
}

void stl_priqueue(int a[], int n, const int k)
{
	reset(a, n);
	auto ts = getTime();

	priority_queue<int> pri_queue;
	for (int i = 0; i < k; i++) {
		pri_queue.push(a[i]);
	}

	int maxe = pri_queue.top();
	for (int i = k; i < n; i++) {
		if (a[i] < maxe) {
			pri_queue.pop();
			pri_queue.push(a[i]);
			maxe = pri_queue.top();
		}
	}

	printf("stl priority_queue %ld ms, a[%d] = %d\n", getTime() - ts, k, maxe);
}

void max_heap(int a[], int n, const int k)
{
	reset(a, n);
	auto ts = getTime();

	maxheap<int> my_heap(k);
	for (int i = 0; i < k; i++) {
		my_heap.push(a[i]);
	}

	int maxe = my_heap.top();
	for (int i = k; i < n; i++) {
		if (a[i] < maxe) {
			my_heap.pop();
			my_heap.push(a[i]);
			maxe = my_heap.top();
		}
	}

	int sum = accumulate(my_heap.a + 1, my_heap.a + k + 1, 0);
	printf("my max_heap     %4ld ms, a[%d] = %d, sum = %d\n", getTime() - ts, k, maxe, sum);
}

void bucket_sort(int a[], int n, const int k)
{
	reset(a, n);
	auto ts = getTime();

	const unsigned int segment = 20;
	const unsigned int segsize = 1 << (32 - segment);
	int bucket[segsize] = {0}; //32k
	for (int i = 0; i < n; i++) {
		bucket[a[i] >> segment] ++;
	}

	int maxe = 1 << segment;
	for (int i = 0, total = 0; i < segsize; i++) {
		total += bucket[i];
		if (total >= k) {
			maxe = (i + 1) << segment;
			break;
		}
	}

	int j = 0;
	for (int i = 0; i < n; i++)	{
		if (a[i] < maxe) {
			a[j++] = a[i];
			if (j >= k * 2) {
				std::sort(a, a + j);
				maxe = a[(j = k) - 1];
			}
		}
	}

	std::sort(a, a + j);
	int sum = accumulate(a, a + k, 0);
	printf("bucket_sort      %4ld ms, a[%d] = %d, sum = %d\n", getTime() - ts, j, a[k - 1], sum);
}

void merge_sort(int a[], int n, const int k)
{
	reset(a, n);

	auto ts = getTime();

	std::sort(a, a + k);

	int* best_a = a + k;
	int maxe = a[k - 1];
	int bestn = 0;
	for (int i = k; i < n; i++) {
		if (a[i] < maxe) {
			best_a[bestn++] = a[i];
			if (bestn >= k / 16 + 10) {
				std::sort(best_a, best_a + bestn);
				std::inplace_merge(a, best_a, best_a + bestn);
				//std::sort(a, best_a + bestn);
				maxe = a[k - 1];
				bestn = 0;
			}
		}
	}

	std::sort(a, a + k + bestn);
	maxe = a[k - 1];
	int sum = accumulate(a, a + k, 0);
	printf("sort-merge      %4ld ms, a[%d] = %d, sum = %d\n", getTime() - ts, k, maxe, sum);
}

void merge_array(int a[], int b[], const int k)
{
	std::sort(b, b + k);
	if (a[0] >= b[k - 1]) {
		memcpy(a, b, sizeof(int) * k);
		return;
	}

	int i = k - 1, s = 128;
	int j = k - 1 + s;
	while (i > s && a[i - s] > b[j - i]) {
		i -= s;
	}

	for (; i >= 0; i--) {
		j = k - 1 - i;
		//a[0, i]  b[0, j - 1]
		if (a[i] <= b[j]) {
#ifndef ME
			int m = k - 1; j --;
			while (j >= 0) {
				if (a[i] <= b[j])
					a[m --] = b[j--];
				else
					a[m --] = a[i--];
			}
#else
			memcpy(a + i + 1, b, sizeof(int) * j);
			std::inplace_merge(a, a + i + 1, a + k);
#endif
			break;
		}
	}
}

void merge_sort2(int a[], int n, const int k)
{
	reset(a, n);

	auto ts = getTime();
	std::sort(a , a + k);

	int* best_a = a + k;
	int bestn = 0;
	int maxe = a[k - 1];
	for (int i = k; i < n; i++) {
		if (a[i] < maxe) {
			best_a[bestn++] = a[i];
			if (bestn == k) {
				merge_array(a, best_a, k);
				maxe = a[k - 1];
				bestn = 0;
			}
		}
	}

	std::sort(a, best_a + bestn);
	maxe = a[k - 1];
	int sum = accumulate(a, a + k, 0);
	printf("sort-merge2     %4ld ms, a[%d] = %d, sum = %d\n", getTime() - ts, k, maxe, sum);
}

static void cpuidInfo(int cpuinfo[4], int id)
{
#if _MSC_VER >= 1600 //2010
	__cpuidex(cpuinfo, id, 0);
#elif _MSC_VER >= 1400 //2005
	__cpuid(cpuinfo, id);
#elif __GNUC__ || __TINYC__
	__asm__ ( "cpuid\n" :"=a"(cpuinfo[0]),"=b"(cpuinfo[1]),"=c"(cpuinfo[2]),"=d"(cpuinfo[3]) :"a"(id));
#elif ASM_X86
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
#endif
}

static int getIntelL3info(int *l3size /*, uint *assoc, uint *linesize*/)
{
	int regs[4];
	int i;

	cpuidInfo(regs, 0); /* Maximum Input Value */
	int max_leaf = regs[0];
	if (max_leaf < 2) {
		return -1; /* no way to find L3 cache info */
	}

	cpuidInfo(regs, 1); /* Additional Information */
	int family = (regs[0] >> 8) & 0xF;
	int model = (regs[0] >> 4) & 0xF;

	cpuidInfo(regs, 2); /* Cache and TLB Information */

	regs[0] &= 0xFFFFFF00; /* least significant byte of EAX is invalid */
	for (i = 0; i < 4; i++) {
		if (regs[i] < 0) { /* invalid if most significant bit set */
			regs[i] = 0;
		}
	}

	unsigned char *descriptors = (unsigned char *) regs;

	const int mb = 1024;

	#define RETINFO(s, a, l) *l3size = (s); /* *assoc = (a); *linesize = (l);*/ return 0

	int use_leaf_4 = 0;
	for (i = 0; i < 32; i++) {
		switch(descriptors[i]) {
			case 0x22: RETINFO(512 , 4, 64);
			case 0x23: RETINFO(1 * mb, 8, 64);
			case 0x25: RETINFO(2 * mb, 8, 64);
			case 0x29: RETINFO(4 * mb, 8, 64);
			case 0x40: RETINFO(0, 0, 0); /* no L3 cache */
			case 0x46: RETINFO(4 * mb, 4, 64);
			case 0x47: RETINFO(8 * mb, 8, 64);
			case 0x49:
					 if (family == 0x0F && model == 0x06) {
						RETINFO(4 * mb, 16, 64);
					 }
					break;
			case 0x4A: RETINFO(6 * mb, 12, 64);
			case 0x4B: RETINFO(8 * mb, 16, 64);
			case 0x4C: RETINFO(12 * mb, 12, 64);
			case 0x4D: RETINFO(16 * mb, 16, 64);
			case 0xD0: RETINFO(512 , 4, 64);
			case 0xD1: RETINFO(1 * mb, 4, 64);
			case 0xD2: RETINFO(2 * mb, 4, 64);
			case 0xD6: RETINFO(1 * mb, 8, 64);
			case 0xD7: RETINFO(2 * mb, 8, 64);
			case 0xD8: RETINFO(4 * mb, 8, 64);
			case 0xDC: RETINFO(1 * mb + 512, 12, 64);
			case 0xDD: RETINFO(3 * mb, 12, 64);
			case 0xDE: RETINFO(6 * mb, 12, 64);
			case 0xE2: RETINFO(2 * mb, 16, 64);
			case 0xE3: RETINFO(4 * mb, 16, 64);
			case 0xE4: RETINFO(8 * mb, 16, 64);
			case 0xEA: RETINFO(12 * mb, 24, 64);
			case 0xEB: RETINFO(18 * mb, 24, 64);
			case 0xEC: RETINFO(24 * mb, 24, 64);
			case 0xFF:
					use_leaf_4 = 1;
					break;
		}
	}

	if (!use_leaf_4 || max_leaf < 4) {
		return -1; /* failed, no L3 info found */
	}

	i = 0;
	while(1) {
#if _MSC_VER >= 1400
		__cpuidex(regs, 4, i); /* Deterministic Cache Parameters */
#else
		cpuidInfo(regs, 4); /* Deterministic Cache Parameters */
#endif
		if ((regs[0] & 0x1F) == 0) {
			RETINFO(0, 0, 0); /* no L3 cache */
		}
		if (((regs[0] >> 5) & 0x7) == 3) {
			int lsize = (regs[1] & 0xFFF) + 1;
			int partitions = ((regs[1] >> 12) & 0x3FF) + 1;
			int ways = ((regs[1] >> 22) & 0x3FF) + 1;
			int sets = regs[2] + 1;
			RETINFO(ways * partitions * lsize * sets, ways, lsize);
		}
		i++;
	}

	return 0;
}

static int getCpuInfo()
{
	char vendor[0x40] = {0};
	int (*pTmp)[4] = (int(*)[4])vendor;
	cpuidInfo(*pTmp ++, 0x80000002);
	cpuidInfo(*pTmp ++, 0x80000003);
	cpuidInfo(*pTmp ++, 0x80000004);

	for (int i = 0; vendor[i]; i ++) {
		if (vendor[i] != ' ' || vendor[i + 1] != ' ')
			putchar(vendor[i]);
	}

	int cpuinfo1[4]; cpuidInfo(cpuinfo1, 0x80000005); //work for amd cpu
	int cpuinfo2[4]; cpuidInfo(cpuinfo2, 0x80000006);

	int l3Size = 0;
	if (strstr(vendor, "Intel") != NULL && !strstr(vendor, "AMD")) {
		getIntelL3info(&l3Size);
	}
	else {
		l3Size = (cpuinfo2[3] >> 18) * 512;
	}

	printf(" Cpu Cache L1Dsize = %d kb, L2Size/L3Size = %d/%d kb\n", cpuinfo1[2] >> 24, cpuinfo2[2] >> 16, l3Size);
	return 0;
}

static int arr[MAXN];
int main(int argc, char* argv[])
{
	int maxn = MAXN, k = 1000, type = 0;
	srand(time(NULL));

	getCpuInfo();
	printf("cmd:topk k(<=%d) n(<=%d) type(0 rand,1 decrease,2 increase,3 wavy)\n\n", MAXK, MAXN);

//	if (argc > 1) { type = atoi(argv[1]); }
	if (argc > 1) { k  = atoi(argv[1]);   if (k <= 0  || k > maxn || k > MAXK) k = MAXK;   }
	if (argc > 2) { maxn = atoi(argv[2]); if (maxn <= 0 || maxn > MAXN) maxn = MAXN; }
	assert(k  <= MAXK && k < MAXN / 2);
	assert(maxn <= MAXN);

	for (int j = 3; j >= 0; j --) {
		int s = rand();
		type = j;
		for (int i = 0; i < maxn; i++) {
			int r = rand() * rand() + s * rand() + i;
			if (type == 1)
				r = -s + i;
			else if (type == 2)
				r = maxn - i - s;
			else if (type == 3)
				r = (s + i) * (i + j);
			if (r < 0)
				r = 0 - r;
			buff[i] = r;
		}

		printf("maxn = %d, topk = %d, type = %d\n", maxn, k, type);
//		stl_nth(arr, maxn, k);
//		bucket_sort(arr, maxn, k);
		max_heap(arr, maxn, k);
//		stl_priqueue(arr, maxn, k);
		stl_makeheap(arr, maxn, k);
		merge_sort2(arr, maxn, k);
//		merge_sort(arr, maxn, k);
		putchar('\n'); putchar('\n');
	}

	return 0;
}

