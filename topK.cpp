#include <algorithm>
#include <numeric>
#include <queue>
#include <cstdio>
//#include <cstring>
//#include <cstdlib>
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
	maxheap(int maxk)
	{
		size = 0;
		a = (T*)malloc(sizeof(T) * (2*maxk + 2));
		a[0] = INT_MAX;
		for (int i = 0; i <= maxk; i++)
			a[maxk + i] = INT_MIN;
	}

	~maxheap() { free(a); }
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

	T *a;
	int size;
};

static int buff[MAXN + 2];

void rand_swap(int a[], int n, int k)
{
	const int step = n / k;
//	std::sort(a, a + k);
	for (int i = 1; i < k; i ++) {
#if 0
		int h = rand() % k, t = i * step + rand() % step;
		if (a[h] > a[t])
			std::swap(a[h], a[t]);
		else if (a[i] > a[t])
			std::swap(a[i], a[t]);
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
	clock_t ts = getTime();

	std::sort(a, a + n);
	printf("stl sort       %4ld ms, a[%d] = %d\n\n", getTime() - ts, k, a[k - 1]);
}

void stl_nth(int a[], int n, const int k)
{
	reset(a, n);
	clock_t ts = getTime();

	std::nth_element(a, a + k, a + n);
	int maxe = *std::max_element(a, a + k);
	int sum =  accumulate(a, a + k, 0);
	printf("stl nth_element  %4ld ms, a[%d] = %d, sum = %d\n", getTime() - ts, k, maxe, sum);
}

void stl_makeheap(int a[], int n, const int k)
{
	reset(a, n);
	clock_t ts = getTime();

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
	clock_t ts = getTime();

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

	printf("stl pri_queue   %4ld ms, a[%d] = %d\n", getTime() - ts, k, maxe);
}

void max_heap(int a[], int n, const int k)
{
	reset(a, n);
	clock_t ts = getTime();

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
	clock_t ts = getTime();

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

	clock_t ts = getTime();

	std::sort(a, a + k);

	int* best_a = a + k;
	int maxe = a[k - 1], bestn = 0;
	for (int i = k; i < n; i++) {
		if (a[i] < maxe) {
			best_a[bestn++] = a[i];
			if (bestn >= k / 8) {
				std::sort(best_a, best_a + bestn);
//				if (a[k - 1] <= best_a[0] && bestn == k) {
//					memcpy(a, best_a, sizeof(a[0]) * bestn);
//					a += k; best_a = a + k;
//				} else
				std::inplace_merge(a, best_a, best_a + bestn);
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

	int i = k / 100 * 52 + 16, s = 16;
	int j = k - 1 + s;

	if (a[k / 2] > b[k / 2])
		i = k / 2;
	else if (j < i || a[i - s] <= b[j - i])
		i = k - 1;
	while (i > s && a[i - s] > b[j - i])
		i -= s;
	for (; i >= 0; i--) {
		j = k - 1 - i;
		if (a[i] > b[j])
			continue ;

		//merge a[0, i]/b[0, j - 1]
#ifndef ME
		int m = k - 1; j --;
		while (j >= 0)
			a[m --] = a[i] <= b[j] ? b[j --] : a[i --];
#else
		memcpy(a + i + 1, b, sizeof(int) * j);
		std::inplace_merge(a, a + i + 1, a + k);
#endif
		return;
	}
}

void merge_sort2(int a[], int n, const int k)
{
	reset(a, n);

	clock_t ts = getTime();
	std::sort(a , a + k);
	a[-1] = INT_MIN;

	int* best_a = a + k;
	int maxe = a[k - 1], bestn = 0;
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

static void printInfo()
{
	const char* sepator =
		"------------------------------------------------------------------------------------------------------------";
	puts(sepator);
	puts("Copyright (C) by 2010-2018 Huang Yuanbing 22738078@qq.com/bailuzhou@163.com\n");

	char buff[500];
	char* info = buff;
#ifdef __clang__
	info += sprintf(info, "clang %s", __clang_version__); //vc/gcc/llvm
#if __llvm__
	info += sprintf(info, " on llvm/");
#endif
#endif

#if _MSC_VER
	info += sprintf(info, "Compiled by vc++ %d", _MSC_VER);
#elif __GNUC__
	info += sprintf(info, "Compiled by gcc %d.%d.%d", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#elif __INTEL_COMPILER
	info += sprintf(info, "Compiled by intel c++ %d", __INTEL_COMPILER);
#elif __TINYC__
	info += sprintf(info, "Compiled by tcc %d", __TINYC__);
#endif

#if __cplusplus
	info += sprintf(info, " __cplusplus = %d", (int)__cplusplus);
#endif

#if __x86_64__ || __amd64__ || _M_X64 || __amd64 || __x86_64
	info += sprintf(info, " x86-64");
#elif __i386__ | _M_IX86 | _X86_ | __i386
	info += sprintf(info, " x86");
#elif __arm64__
	info += sprintf(info, " arm64");
#elif __arm__
	info += sprintf(info, " arm");
#else
	info += sprintf(info, " unknow");
#endif

	puts(buff);
	puts(sepator);
}

int main(int argc, char* argv[])
{
	int maxn = MAXN, k = 1000, type = 0;
	srand(time(NULL));
	printf("\ncmd:topk k(<=%d) n(<=%d) type[0 rand,1 decrease,2 increase,3 wavy, 4 dup]\n\n", MAXK, MAXN);
	printInfo();

	if (argc > 1) { k = atoi(argv[1]); }
	if (argc > 2) { maxn = atoi(argv[2]); if (maxn <= 0 || maxn > MAXN) maxn = MAXN; }

	int* arr = (int *)malloc(sizeof(int) * maxn);
	assert(k < MAXN / 2 && maxn <= MAXN);

	for (int j = 0; j <= 4; j ++) {
		int s = rand(), r = 0;
		type = j;
		for (int i = 0; i < maxn; i++) {
			if (i % RAND_MAX == 0)
				srand(time(NULL));
			if (type == 0)
				r = (rand() << 16) + rand(); //
			else if (type == 1)
				r =  i;
			else if (type == 2)
				r = maxn - i;
			else if (type == 3)
				r = ((s + i) * rand()) % (MAXN + 1);
			else if (type == 4)
				r = (s - i) * i;
			if (r < 0)
				r = -r;
			buff[i] = r;
		}

		printf("maxn = %d, topk = %d, type = %d\n", maxn, k, type);
#if __cplusplus
		stl_nth(arr, maxn, k);
//		bucket_sort(arr, maxn, k);
		max_heap(arr, maxn, k);
		stl_priqueue(arr, maxn, k);
		stl_makeheap(arr, maxn, k);
#endif
		merge_sort2(arr, maxn, k);
		merge_sort(arr, maxn, k);
		putchar('\n'); putchar('\n');
	}

	return 0;
}

// https://www.jdoodle.com/online-compiler-c++
// https://repl.it/repls/ProductivePowerlessSlope
// https://www.tutorialspoint.com/online_cpp_compiler.php 8
// http://rextester.com/l/cpp_online_compiler_gcc 8
// https://www.onlinegdb.com/ 9
// https://www.ideone.com/3rCdob
