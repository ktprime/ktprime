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
#if _WIN32
	typedef int64_t stype;
#else
	typedef long stype;
#endif
	stype max_v = INT64_MAX;
#endif

static int AK = 4;
static int AK2 = 1;

typedef unsigned int uint;

template<class T>
class max_heap
{
public:
	max_heap(int k)
	{
		size = 0;

#if POP_PUSH || NC
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
		uint p = 1, l = 2, r = 3;

#if !NC
		while (l <= size) {
#else
		while (true) {
#endif
			const T c = a[l] >= a[r] ? l : r;
			if (v >= a[c])
				break;

			l = c * 2; r = l + 1;
			a[p] = a[c]; p = c;
		}

		a[p] = v;
		return a[1];
	}

	//private:

	T *a;
	size_t size;
};

void rand_swap(stype a[], int n, int k)
{
#if 0
	const int step = n / k;
	//	std::sort(a, a + k);
	for (int i = 1; i < k; i ++) {

		int h = rand() % k, t = i * step + rand() % step;
		if (a[h] > a[t])
			std::swap(a[h], a[t]);
		else if (a[i] > a[t])
			std::swap(a[i], a[t]);
	}
#endif
}

void reset(stype a[], int n, int k)
{
	memcpy(a, a + n, k * sizeof(a[0]));
	rand_swap(a, n, 10000);
}

void check(const stype a[], int n, int k)
{
	if (a[2*n] == (stype)(-1)) {
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

static clock_t NowMs()
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
#if 0
void stl_sort(stype a[], int n, const int k)
{
	reset(a, n, k);
	clock_t ts = NowMs();

	std::sort(a, a + n);
	printf("  stl sort       %5ld ms, a[%d] = %ld\n\n", NowMs() - ts, k, (int64_t)a[k - 1]);
}
#endif

void stl_nth(stype a[], int n, const int k)
{
	reset(a, n, n);
	clock_t ts = NowMs();

#if 1
	std::nth_element(a, a + k, a + n);
	std::sort(a, a + k);
#else
	std::partial_sort(a, a + k, a + n);
#endif

	stype maxe = a[k - 1];
	int64_t sum = std::accumulate(a, a + k, 0);
	printf("  stl nth_element %5ld ms, a[%d] = %ld, sum = %ld\n", NowMs() - ts, k, (int64_t)maxe, sum);

	memcpy(a + n * 2, a, k * sizeof(a[0]));
	check(a, n, k);
}

void stl_priqueue(stype a[], int n, const int k)
{
//	reset(a, n, k);
	clock_t ts = NowMs();

	std::priority_queue<stype> pri_queue;
	for (int m = 0; m < k; m++) {
		pri_queue.push(a[m + n]);
	}

	stype maxe = pri_queue.top();
	for (int i = n + k; i < n * 2; i++) {
		if (a[i] < maxe) {
			pri_queue.pop();
			pri_queue.push(a[i]);
			maxe = pri_queue.top();
		}
	}

	ts = NowMs() - ts;
	for (int j = 0; j < k; j ++) {
		a[j] = pri_queue.top();
		pri_queue.pop();
	}

	int64_t sum = std::accumulate(a, a + k, 0);
	std::sort(a, a + k);
	printf("  stl pri_queue   %5ld ms, a[%d] = %ld, sum = %ld\n", ts, k, (int64_t)maxe, sum);

	check(a, n, k);
}

void stl_makeheap(stype a[], int n, const int k)
{
	reset(a, n, k);
	clock_t ts = NowMs();

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
	std::sort(a, a + k);

	printf("  stl make_heap   %5ld ms, a[%d] = %ld, sum = %ld\n", NowMs() - ts, k, (int64_t)maxe, sum);

	check(a, n, k);
}

void mymax_heap(stype a[], int n, const int k)
{
	clock_t ts = NowMs();

	max_heap<stype> my_heap(k);
#if 0
	memcpy(my_heap.a + 1, a + n, k * sizeof(a[0]));
	std::sort(my_heap.a + 1, my_heap.a + 1 + k, std::greater<stype>());
	my_heap.size = k;
#else
	for (int i = 0; i < k; i++)
		my_heap.push(a[n + i]);
#endif

	stype maxe = my_heap.top();
	for (int i = n + k; i < n * 2; i++) {
		if (a[i] < maxe) {
#if POP_PUSH
			my_heap.pop(); my_heap.push(a[i]);
			maxe = my_heap.top();
#else
			maxe = my_heap.pop_push(a[i]);
#endif
		}
	}

	int64_t sum = std::accumulate(my_heap.a + 1, my_heap.a + k + 1, 0);
	memcpy(a, my_heap.a + 1, k * sizeof(a[0]));
	std::sort(a, a + k);

	printf("  ktprime max_heap%5ld ms, a[%d] = %ld, sum = %ld\n", NowMs() - ts, k, (int64_t)maxe, sum);

	check(a, n, k);
}
#endif

void bucket_sort(uint a[], uint n, uint k)
{
	reset((stype*)a, n, n);

	clock_t ts = NowMs();

#if I32
	//bug for zero
	for (int i = 0; i < n; i++) a[i] += INT_MAX;
#endif

	constexpr uint segbits = 20; //16 - 22 why
	constexpr uint segsize = 1 << (32 - segbits);
	constexpr uint min_bucket = segsize >> (6 + 1);

	//set bucket fill in L1 cache
	uint bucket[segsize] = { 0 }; //
	uint m = 0;
	for (uint finds = k * 3; m < n; m += 1) {
		const uint bindex = a[m] >> segbits;
		bucket[bindex] ++;
#if DF
		////try find a small range
		if (bindex < min_bucket && finds -- == 0) {
#ifndef DP
			printf("ration = %d%%, m = %u, min_bucket = %u\n", m * 100 / n, m, min_bucket);
#endif
			break;
		}
#endif
	}

#if 0
	if (m >= n / 4) {
		//do a full find for some bad case
		for (; m < n; m++) {
			const uint bindex = a[m] >> segbits;
			bucket[bindex] ++;
		}
	}
#endif

	uint maxe = UINT_MAX;
	for (uint i = 0, sums = 0; i < segsize; i++) {
		sums += bucket[i];
		if (sums >= k) {
			maxe = ((i + 1) << segbits);
			break;
		}
	}

	uint sum_size = 0, min_size = k < 100 ? 4 * k + 32 : k * 3 / 2;
	for (uint i = 0; i < n; i++) {
		if (a[i] >= maxe)
			continue ;

		a[sum_size++] = a[i];
		if (sum_size > min_size) {
			std::sort(a, a + sum_size);
			maxe = a[(sum_size = k) - 1];
		}
	}

	std::sort(a, a + sum_size);

#if I32
	for (int i = 0; i < sum_size; i++) a[i] -= INT_MAX;
#endif

	int64_t sum = std::accumulate(a, a + k, 0);

	printf("  bucket_sort     %5ld ms, a[%d] = %d, sum = %ld\n", NowMs() - ts, k, (int)a[k - 1], sum);
	check((stype*)a, n, k);
}

static void swap_array(stype a[], const int k)
{
	for (int i = 0; i < k / 2; i ++) {
		std::swap(a[i], a[k - i - 1]);
	}
}

/*****
find kth element's positon of two sorted array a/b
*/
static int find_kth(stype a[], const int m, stype b[], const int n, const int kth)
{
	int i = 0, j = 0, k = kth;
	while (k > 0) {
		const int p = (k - 1) >> 1;
		if (j + p >= n || (i + p < m && a[i + p] < b[j + p]))
			i += p + 1;
		else
			j += p + 1;
		k -= p + 1;
	}

	assert(i + j == kth);

//	while (a[i - 1] > b[j]) i--, j++;
//	while (a[i + 1] < b[j]) i ++, j--;

	return (j >= n || (i < m && a[i] <= b[j])) ? i : i - 1;
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

	int t = n - 1;
	int i = find_kth(a, n, b, m, t);
	int j = t - 1 - i;

	//merge a[0, i]/b[0, j] into a[0, n - 1]
	//TODO: check i overflow
	while (j >= 0) {
		a[t --] = a[i] <= b[j] ? b[j --] : a[i --];
	}
}

const int l2_cpu_size = 1024 * 1024, l1_cpu_size = 32 * 1024;
void merge_sort(stype a[], int n, const int k)
{
	reset(a, n, k);
	clock_t ts = NowMs();

	std::sort(a, a + k);
	a[-1] = ~max_v;

	stype* ax_a = a + k;
	stype maxe = a[k - 1];
	int bestn = k / AK2, ax_n = 0;
	if (bestn > l2_cpu_size / sizeof(a[0]))
		bestn = l2_cpu_size / sizeof(a[0]);
	if (bestn < 100)
		bestn = k * 4 + 101;

	for (int i = n + k; i < n + n; i++) {
		if (a[i] >= maxe) {
			continue ;
		}

		ax_a[ax_n++] = a[i];
		if (ax_n == bestn) {
			merge_array(a, ax_a, k, ax_n);
			maxe = a[k - 1];
			ax_n = 0;
		}
	}

//	std::sort(a, ax_a + ax_n);
//	std::sort(ax_a, ax_a + ax_n); std::inplace_merge(a, ax_a, ax_a + ax_n);
//	std::partial_sort(a, ax_a, ax_a + ax_n);
	merge_array(a, ax_a, k, ax_n);

	maxe = a[k - 1];
	int64_t sum = std::accumulate(a, a + k, 0);
	printf("  merge_sort      %5ld ms, a[%d] = %ld, sum = %ld\n", NowMs() - ts, k, (int64_t)maxe, sum);
	check(a, n, k);
}

void merge_inplace(stype a[], int n, const int k)
{
	reset(a, n, k);
	clock_t ts = NowMs();

	std::sort(a, a + k);
	stype* ax_a = a + k;
	stype maxe = a[k - 1];
	int bestn = k / AK, ax_n = 0;

	if (bestn > l2_cpu_size / sizeof(a[0]))
		bestn = l2_cpu_size / sizeof(a[0]);
	if (bestn < 100)
		bestn = k * 4 + 100;

	for (int i = n + k; i < n * 2; i++) {
		if (a[i] >= maxe) {
			continue ;
		}
		ax_a[ax_n++] = a[i];
		if (ax_n < bestn) {
			continue;
		}
#if 0
		if (std::is_sorted(ax_a, ax_a + ax_n, std::greater<stype>()))
			std::reverse(ax_a, ax_a + ax_n);
		else if (!std::is_sorted(ax_a, ax_a + ax_n))
#endif
			std::sort(ax_a, ax_a + ax_n);

		if (a[0] >= ax_a[ax_n - 1]) {
			if (ax_n >= k)
				memcpy(a, a + ax_n, sizeof(a[0]) * k);
			else {
				memmove(a + ax_n, a, sizeof(a[0]) * (k - ax_n));
				memcpy(a, ax_a, sizeof(a[0]) * ax_n);
				//swap to a big range
				if (k > 100)
					bestn = k;
			}
		} else {
			//TODO optimize inplace_merge without merge data after ax_a
			std::inplace_merge(a, ax_a, ax_a + ax_n);
		}
		maxe = a[k - 1];
		ax_n = 0;
	}

	std::sort(a, ax_a + ax_n);
//	std::sort(ax_a, ax_a + ax_n); std::inplace_merge(a, ax_a, ax_a + ax_n);

	maxe = a[k - 1];
	int64_t sum = std::accumulate(a, a + k, 0);
	printf("  merge_inplace   %5ld ms, a[%d] = %ld, sum = %ld\n", NowMs() - ts, k, (int64_t)maxe, sum);
	check(a, n, k);
}

static void printInfo()
{
	const char* sepator =
		"------------------------------------------------------------------------------------------------------------";
	puts(sepator);
	//	puts("Copyright (C) by 2018-2020 Huang Yuanbing bailuzhou at 163.com\n");

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
	const int max_n = 1024*1024*32;

	printf("\ncmd:topk -k(<=%d) -n(<=%d) -r(1-16) m(1-16) s(shuffle) [type = 0-1 rand, 2-3 wavy, 4-5 rand, 6 decrease, 7 incre]\n\n", 1000000, max_n);
	printInfo();

	int n = max_n, k = max_n / 100, type = 0, shuff = 0;
	for (int i = 1; i < argc; i++)
	{
		char c = argv[i][0];
		int disgi = (c == '-') ? 1 : 0;
		c = argv[i][disgi];
		int r = atoi(argv[i] + disgi + 1);

		if (c == 'k') {
			if (r > 0) { k = r; }
			else if (r == 0) k = rand() % 100000;
		} else if (c == 'n') {
			if (r >= -100 && r < 0) { n = max_n / (-r); }
			else if (r <= 100 && r > 0) { n = max_n * r; }
			else if (r > 0) { n = r; }
		} else if (c == 'r' && r > 0) {
			AK = r;
		} else if (c == 'm' && r > 0) {
			AK2 = r;
		} else if (c == 's') {
			shuff = r;
		}
	}

	if (k > n) {
		k = n;
	}

	stype* arr = static_cast<stype *>(malloc(sizeof(stype) * (n * 2 + k + 8))) + 4;
	stype* buff = arr + n;
	buff[n] = (stype)(-1);

	std::random_device rd;
	std::default_random_engine e(rd());
	std::mt19937_64 rng; rng.seed(time(nullptr));

	std::uniform_int_distribution<int64_t> u(0, INT64_MAX / 2);
	std::normal_distribution<> d(1 << (20 + rand() % 10), 1 << 16);
//	std::exponential_distribution<> p(0.1);

	printf("n = %d, topk = %d, r1 = %d, r2 = %d, shuff = %d\n", n, k, AK, AK2, shuff);
	for (int j = 0; j <= 7; j ++) {
		int64_t r = 0;
		type = j;
		for (size_t i = 1; i <= (size_t)n; i++) {
			if (type == 0) {
				r = rd() + e();
			} else if (type == 2) {
				r = rng();
			} else if (type == 1) {
				r = u(rd) + (1 << 20);
			} else if (type == 3) {
				r = d(rng);
			} else if (type == 4) {
				r = i * i + (512 - rng() % 1024) * (i % 8);
			} else if (type == 5) {
				r = rng() - (int)(i * i) % 1024;
			} else if (type == 6) {
				r = n - i - rng() % 4; // i + 1
			} else if (type == 7) {
				r = i + rng() % 256;
			}
			buff[i] = (stype)r + 1;
		}

		if (shuff) {
			std::shuffle(buff + 1, buff + n, rng);
		}

		printf("type = %d\n", type);
		stl_nth(arr, n, k);

		if (sizeof(arr[0]) <= sizeof(int))
			bucket_sort((uint*)arr, n, k);

#if __cplusplus
		stl_priqueue(arr, n, k);
		stl_makeheap(arr, n, k);
		mymax_heap(arr, n, k);
#endif
		merge_inplace(arr, n, k);
		merge_sort(arr, n, k);
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

