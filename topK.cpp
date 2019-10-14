//https://www.geeksforgeeks.org/k-largestor-smallest-elements-in-an-array/
//https://github.com/siddontang/leetcode/blob/master/median-of-two-sorted-arrays.cpp
//https://github.com/MaskRay/LeetCode/blob/master/median-of-two-sorted-arrays.cc
//
#include <algorithm>
#include <numeric>
#include <queue>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <climits>
#include <random>

#if I32
typedef int stype;
const stype max_v = INT_MAX;
#elif U32
typedef unsigned int stype;
stype max_v = UINT_MAX;
#else
typedef int64_t stype;
stype max_v = INT64_MAX;
#endif

#ifndef AK
static int AK = 4;
#endif

template<class T>
class max_heap
{
public:
	max_heap(int k)
	{
		size = 0;

#if POP_PUSH || CHECK
		a = (T*)malloc(sizeof(T) * (2*k + 3));
		for (int i = 0; i <= k + 1; i++)
			a[k + i] = ~max_v;
#else
		a = (T*)malloc(sizeof(T) * (k + 3));
		a[k + 1] = a[k + 2] = ~max_v;
#endif
		a[0] = max_v;
	}

	~max_heap() { free(a); }
	T top() const { return a[1]; }

	void push(const T x)
	{
		size_t c = ++size;
		size_t p = size / 2;

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
		size_t p = 1, c = 1;

		while (x < a[c]/* && c <= size &&*/) {
			a[p] = a[c];
			p = c;
			c *= 2;
			//if (a[c + 1] > a[c]) c++;
			c += (a[c + 1] > a[c]);
		}
		a[p] = x;
	}

	T pop_push(const T v)
	{
		size_t p = 1, l = 2, r = 3;

#if !CHECK
		while (l <= size) {
#else
		while (true) {
#endif
			const T c = a[l] >= a[r] ? l : r;
			if (v >= a[c])
				break;

			a[p] = a[c], p = c;

			l = c * 2;
			r = l + 1;
		}

		a[p] = v;
		return a[1];
	}

	//private:

	T *a;
	int size;
};

void rand_swap(stype a[], int n, int k)
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

void reset(stype a[], int n, int k)
{
	memcpy(a, a + n, k * sizeof(a[0]));
	rand_swap(a, n, 10000);
}

void check(const stype a[], int n, int k)
{
	if (a[2*n] == -1) {
		return;
	}

	for (int i = 0; i < k; i++) {
		if (a[i] != a[n * 2 + i]) {
			printf("%d %ld != %ld\n", i, (int64_t)a[i], (int64_t)a[n * 2 + i]);
			break;
		}
	}
}

#if __linux__
#include <sys/resource.h>
#elif _WIN32
#include <windows.h>
#endif

static clock_t getTime()
{
#if _WIN32
	FILETIME ptime[4] = {0};
	GetThreadTimes(GetCurrentThread(), &ptime[0], &ptime[1], &ptime[2], &ptime[3]);
	return (ptime[2].dwLowDateTime + ptime[3].dwLowDateTime) / 10000;
	//return clock();
#elif __linux__ || __unix__
	struct rusage rup;
	getrusage(RUSAGE_SELF, &rup);
	long sec = rup.ru_utime.tv_sec + rup.ru_stime.tv_sec;
	long usec = rup.ru_utime.tv_usec + rup.ru_stime.tv_usec;
	return sec * 1000 + usec / 1000;
#elif __linux__ || __unix__
	return clock() / 1000;
#else
	return clock();
#endif
}

#if __cplusplus
void stl_sort(stype a[], int n, const int k)
{
	reset(a, n, k);
	clock_t ts = getTime();

	std::sort(a, a + n);
	printf("  stl sort       %5ld ms, a[%d] = %ld\n\n", getTime() - ts, k, (int64_t)a[k - 1]);
}

void stl_nth(stype a[], int n, const int k)
{
	reset(a, n, n);
	clock_t ts = getTime();

	std::nth_element(a, a + k, a + n);
	std::sort(a, a + k);
	//	std::partial_sort(a, a + k, a + n);

	stype maxe = a[k - 1];
	int64_t sum = std::accumulate(a, a + k, 0);
	printf("  stl nth_element %5ld ms, a[%d] = %ld, sum = %ld\n", getTime() - ts, k, (int64_t)maxe, sum);

	memcpy(a + n * 2, a, k * sizeof(a[0]));
	check(a, n, k);
}

void stl_priqueue(stype a[], int n, const int k)
{
	reset(a, n, k);
	clock_t ts = getTime();

	std::priority_queue<stype> pri_queue;
	for (int m = 0; m < k; m++) {
		pri_queue.push(a[m]);
	}

	stype maxe = pri_queue.top();
	for (int i = n + k; i < n * 2; i++) {
		if (a[i] < maxe) {
			pri_queue.pop();
			pri_queue.push(a[i]);
			maxe = pri_queue.top();
		}
	}

	ts = getTime() - ts;
	for (int j = 0; j < k; j ++) {
		a[j] = pri_queue.top();
		pri_queue.pop();
	}

	int64_t sum = std::accumulate(a, a + k, 0);
	printf("  stl pri_queue   %5ld ms, a[%d] = %ld, sum = %ld\n", ts, k, (int64_t)maxe, sum);

	std::sort(a, a + k);
	check(a, n, k);
}

void stl_makeheap(stype a[], int n, const int k)
{
	reset(a, n, k);
	clock_t ts = getTime();

	std::make_heap(a, a + k);

	stype maxe = a[0];
	for (int i = n + k; i < n * 2; i++) {
		if (a[i] < maxe) {
			std::pop_heap(a, a + k);
			a[k - 1] = a[i];
			std::push_heap(a, a + k);
			maxe = a[0];
		}
	}

	int64_t sum = std::accumulate(a, a + k, 0);
	printf("  stl make_heap   %5ld ms, a[%d] = %ld, sum = %ld\n", getTime() - ts, k, (int64_t)maxe, sum);

	std::sort(a, a + k);
	check(a, n, k);
}

void mymax_heap(stype a[], int n, const int k)
{
	reset(a, n, k);
	clock_t ts = getTime();

	max_heap<stype> my_heap(k);
	for (int i = 0; i < k; i++) {
		my_heap.push(a[i]);
	}

	stype maxe = my_heap.top();
	for (int i = n + k; i < n * 2; i++) {
		if (a[i] < maxe) {
#ifdef POP_PUSH
			my_heap.pop(); my_heap.push(a[i]);
			maxe = my_heap.top();
#else
			maxe = my_heap.pop_push(a[i]);
#endif
		}
	}

	int64_t sum = std::accumulate(my_heap.a + 1, my_heap.a + k + 1, 0);
	printf("  ktprime max_heap%5ld ms, a[%d] = %ld, sum = %ld\n", getTime() - ts, k, (int64_t)maxe, sum);

	memcpy(a, my_heap.a + 1 , k * sizeof(a[0]));
	std::sort(a, a + k);
	check(a, n, k);
}
#endif

void bucket_sort(unsigned int a[], int n, const int k)
{
	reset((stype*)a, n, n);

	clock_t ts = getTime();

#if I32
	for (int i = 0; i < n; i++) a[i] += INT_MAX;
#endif

	const unsigned int segment = 20;
	const unsigned int segsize = 1 << (32 - segment);
	unsigned int bucket[segsize] = {0}, mink = 0; //32k
	for (int m = 0; m < n; m++) {
		const int bindex = a[m] >> segment;
		bucket[bindex] ++;
		if (bindex < 8 && mink ++ >= k) { // why
			break;
		}
	}

	unsigned int maxe = 1 << segment;
	int i = 0, j = 0, bsize = 0;
	for (i = 0, bsize = 0; i < segsize; i++) {
		bsize += bucket[i];
		if (bsize >= k) {
			maxe = (i + 1) << segment;
			break;
		}
	}

	for (i = 0; i < n; i++) {
		if (a[i] < maxe) {
			a[j++] = a[i];
			if (j >= k * 2) {
				std::sort(a, a + j);
				maxe = a[(j = k) - 1];
			}
		}
	}

	std::sort(a, a + j);
#if I32
	for (int i = 0; i < j; i++) a[i] -= INT_MAX;
#endif

	int64_t sum = std::accumulate((stype*)a, (stype*)a + k, 0);
	ts = getTime() - ts;

	printf("  bucket_sort     %5ld ms, a[%d] = %d, sum = %ld\n", ts, j, (int)a[k - 1], sum);
	check((stype*)a, n, k);
}

void swap_array(stype a[], const int k)
{
	for (int i = 0; i < k / 2; i ++) {
		std::swap(a[i], a[k - i - 1]);
	}
}

static stype findKth(stype a[], int m, stype b[], int n, int k)
{
	if(m > n) {
		return findKth(b, n, a, m, k);
	}
	if(m == 0) {
		return b[k - 1];
	}
	if(k == 1) {
		return a[0]< b[0] ? a[0] : b[0];
	}

	int pa = k / 2 < m ? k / 2 : m;
	int pb = k - pa;
	if(a[pa - 1] < b[pb - 1]) {
		return findKth(a + pa, m - pa, b, n, k - pa);
	} else if(a[pa - 1] > b[pb - 1]) {
		return findKth(a, m, b + pb, n - pb, k - pb);
	} else {
		return a[pa - 1];
	}
}

static stype findMedianSortedArrays2(stype A[], int m, stype B[], int n)
{
	int t = m + n;
	if(t & 0x01) {
		return findKth(A, m, B, n, t / 2 + 1);
	} else {
		return (findKth(A, m, B, n, t / 2) + findKth(A, m, B, n, t / 2 + 1)) / 2;
	}
}

/*****
find kth element's positon of a of two sorted array
*/
static int find_kth(stype a[], const int m, stype b[], const int n, const int kth)
{
	int i = 0, j = 0, k = kth;
	while (k > 0) {
		const int p = (k - 1) >> 1;
		if (j + p >= n || i + p < m && a[i + p] < b[j + p])
			i += p + 1;
		else
			j += p + 1;
		k -= p + 1;
	}

	assert(i + j == kth);
	int s = j >= n || i < m && a[i] < b[j] ? i++ : j++;
	return i - 1;
//	if ((m + n) & 1)
//		return s;
//	return j >= n || i < m && a[i] < b[j] ? i : n - j;
}

static void merge_array(stype a[], stype b[], const int n, const int m)
{
	std::sort(b, b + m);
	if (a[0] >= b[m - 1]) {
		if (m >= n)
			memcpy(a, b + m - n, sizeof(a[0]) * n);
		else {
			memmove(a + m, a, sizeof(a[0]) * (n - m));
			memcpy(a, b, sizeof(a[0]) * m);
		}
		return;
	}

	int i = find_kth(a, n, b, m, n - 1);
	int j = n - 2 - i;
	int t = n - 1;

	//merge a[0, i]/b[0, j] into a[0, n - 1]
	while (j >= 0) {
		a[t --] = a[i] <= b[j] ? b[j --] : a[i --];
	}
}

const int l2_cpu_size = 1024, l1_cpu_size = 32;
void merge_sort2(stype a[], int n, const int k)
{
	reset(a, n, k);
	clock_t ts = getTime();

	std::sort(a, a + k);
	a[-1] = a[-2] = ~max_v;

	stype* ax_a = a + k;
	stype maxe = a[k - 1];
	int auxn = k / AK, axn = 0;
	if (auxn * sizeof(a[0]) > l2_cpu_size * 1024)
		auxn = l2_cpu_size / sizeof(a[0]);
	if (auxn % 2 == 1) auxn ++;

	for (int i = n + k; i < n + n; i++) {
		if (a[i] < maxe) {
			ax_a[axn++] = a[i];
			if (axn == auxn) {
				merge_array(a, ax_a, k, axn);
				maxe = a[k - 1];
				axn = 0;
			}
		}
	}

    std::sort(a, ax_a + axn);
//	std::sort(ax_a, ax_a + axn); std::inplace_merge(a, ax_a, ax_a + axn);
//	std::partial_sort(a, ax_a, ax_a + axn);
//	merge_array(a, ax_a, k, axn);

	maxe = a[k - 1];
	int64_t sum = std::accumulate(a, a + k, 0);
	printf("  sort_merge2     %5ld ms, a[%d] = %ld, sum = %ld\n", getTime() - ts, k, (int64_t)maxe, sum);
	check(a, n, k);
}

void merge_sort(stype a[], int n, const int k)
{
	reset(a, n, k);
	clock_t ts = getTime();

	std::sort(a, a + k);
	stype* ax_a = a + k;
	stype maxe = a[k - 1];
	int auxn = k / AK + 10, axn = 0;

	if (auxn * sizeof(a[0]) > l2_cpu_size * 1024)
		auxn = l2_cpu_size / sizeof(a[0]);
	if (auxn % 2 == 1)
		auxn ++;

	for (int i = n + k; i < n * 2; i++) {
		if (a[i] >= maxe) {
			continue ;
		}

		ax_a[axn++] = a[i];
		if (axn == auxn) {
			//if (std::is_sorted(ax_a, ax_a + axn, std::greater<stype>()))
			//swap_array(ax_a, axn);
			//if (!std::is_sorted(ax_a, ax_a + axn))
			std::sort(ax_a, ax_a + axn);
			if (a[0] >= ax_a[axn - 1]) {
				if (axn >= k)
					memcpy(a, ax_a + axn - k, sizeof(a[0]) * k);
				else {
					memmove(a + axn, a, sizeof(a[0]) * (k - axn));
					memcpy(a, ax_a, sizeof(a[0]) * axn);
				}
			} else {
				std::inplace_merge(a, ax_a, ax_a + axn);
			}
			maxe = a[k - 1];
			axn = 0;
		}
	}

	std::sort(a, ax_a + axn);
//	std::sort(ax_a, ax_a + axn); std::inplace_merge(a, ax_a, ax_a + axn);

	maxe = a[k - 1];
	int64_t sum = std::accumulate(a, a + k, 0);
	printf("  sort_merge      %5ld ms, a[%d] = %ld, sum = %ld\n", getTime() - ts, k, (int64_t)maxe, sum);
	check(a, n, k);
}

static void printInfo()
{
	const char* sepator =
		"------------------------------------------------------------------------------------------------------------";
	puts(sepator);
	//	puts("Copyright (C) by 2018-2020 Huang Yuanbing 22738078@qq.com/bailuzhou@163.com\n");

	char cbuff[256];
	char* info = cbuff;
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
#endif

#if __cplusplus
	info += sprintf(info, " __cplusplus = %d", static_cast<int>(__cplusplus));
#endif

#if __x86_64__ || __amd64__ || _M_X64 || __amd64 || __x86_64
	info += sprintf(info, " x86-64");
#elif __i386__ || _M_IX86 || _X86_ || __i386
	info += sprintf(info, " x86");
#elif __arm64__
	info += sprintf(info, " arm64");
#elif __arm__
	info += sprintf(info, " arm");
#else
	info += sprintf(info, " unknow");
#endif

	puts(cbuff);
	puts(sepator);
}

int main(int argc, char* argv[])
{
	srand(time(nullptr));
	int const MAXN = 10000*10000;
	printf("\ncmd:topk -k(<=%d) -n(<=%d) -r(2-16) [type = 0-1 rand, 2-3 wavy, 4-5 rand, 6 decrease, 7 incre]\n\n", 1000000, MAXN);
	printInfo();

	int n = MAXN, k = MAXN / 1000, type = 0;
	for (int i = 1; i < argc; i++)
	{
		char c = argv[i][0];
		int disgi = (c == '-') ? 1 : 0;
		c = argv[i][disgi];
		int r = atoi(argv[i] + disgi + 1);

		if (c == 'k') {
			if (r > 0) { k = r; }
		} else if (c == 'n') {
			if (r >= -100 && r < 0) { n = MAXN / (-r); }
			else if (r <= 100 && r > 0) { n = MAXN * r; }
			else if (r > 0) { n = r; }
		} else if (c == 'r' && r > 0) {
			AK = r;
		}
	}

	if (k > n) {
		k = n;
	}

	stype* arr = static_cast<stype *>(malloc(sizeof(stype) * (n * 2 + k + 8))) + 4;
	stype* buff = arr + n;
	buff[n] = -1;

	std::default_random_engine e(time(nullptr));
	std::mt19937_64 rng; rng.seed(time(nullptr));
	std::uniform_int_distribution<int64_t> u(0, INT64_MAX);
	std::normal_distribution<> d(1 << 30, 1 << 16);
//	std::exponential_distribution<> p(0.1);

	printf("n = %d, topk = %d, r = %d\n", n, k, AK);
	for (int j = 0; j <= 7; j ++) {
		int64_t r = 0;
		type = j;
		for (size_t i = 1; i <= n; i++) {
			if (type == 0) {
				r = e();
			} else if (type == 1) {
				r = u(e);
			} else if (type == 2) {
				r = d(rng);
			} else if (type == 3) {
				r = rng() * i - i % (1 << 16);
			} else if (type == 4) {
				r = i + (512 - e() % 1024) * (i % 8);
			} else if (type == 5) {
				r = e() - (int)(i * i) % 1024;
			} else if (type == 6) {
				r = n - i - rng() % 4; // i + 1
			} else if (type == 7) {
				r = i + e() % 256;
			}

			buff[i] = (stype)r + 1;
		}

		printf("type = %d\n", type);
		stl_nth(arr, n, k);

		if (sizeof(arr[0]) <= sizeof(int)) bucket_sort((unsigned int*)arr, n, k);

#if __cplusplus
		stl_priqueue(arr, n, k);
		stl_makeheap(arr, n, k);
		mymax_heap(arr, n, k);
#endif
		merge_sort(arr, n, k);
		merge_sort2(arr, n, k);
		putchar('\n');
	}

	free(arr - 4);
	return 0;
}

// https://www.jdoodle.com/online-compiler-c++
// https://repl.it/repls/ProductivePowerlessSlope
// https://www.tutorialspoint.com/online_cpp_compiler.php
// http://rextester.com/l/cpp_online_compiler_gcc
// https://www.onlinegdb.com/
// https://www.ideone.com/3rCdob

