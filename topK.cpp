/***
maxn = 100000000, topk = 10000, type = 0
stl nth_element  1144 ms
radix_sort        200 ms
stl priority_queue 108 ms
my max_heap      152 ms
stl make_heap    108 ms
sort-merge2      112 ms
sort-merge       108 ms

maxn = 100000000, topk = 10000, type = 1
stl nth_element   128 ms
radix_sort        400 ms
stl priority_queue 92 ms
my max_heap      140 ms
stl make_heap     92 ms
sort-merge2       96 ms
sort-merge        96 ms

maxn = 100000000, topk = 10000, type = 3
stl nth_element   368 ms
radix_sort        216 ms
stl priority_queue 116 ms
my max_heap      160 ms
stl make_heap    112 ms
sort-merge2      116 ms
sort-merge       112 ms

maxn = 100000000, topk = 10000, type = 2
stl nth_element   128 ms
radix_sort        400 ms
stl priority_queue 96 ms
my max_heap      140 ms
stl make_heap     96 ms
sort-merge2       96 ms
sort-merge       136 ms
*/


#include <algorithm>
#include <numeric>
#include <queue>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cassert>

#define  MAXN  10000*10000
#define  MAXK  10000

using namespace std;

template<class T>
class maxheap
{
public:
	maxheap() { size = 0; a[0] = 1 << 30; }
	T top()   { return a[1]; }

	void push(const T& x)
	{
		int c = ++size;
		int p = size / 2;

		while (x > a[p] && p >= 1) {
			a[c] = a[p];
			c = p;
			p /= 2;
		}
		a[c] = x;
	}

	void pop()
	{
		const T& x = a[size--];
		int p = 1;
		int c = 1;

		while (c <= size && x < a[c]) {
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

void reset(int a[], int n)
{
	memcpy(a, buff, n * sizeof(a[0]));
}

void rand_swap(int a[], int n, int k)
{
	const int step = n / k;
//	std::sort(a, a + k);
	for (int i = 1; i < k; i ++) {
#if 1
		int h = rand() % k, t = i * step + rand() % step;
		if (a[h] > a[t])
			std::swap(a[h], a[t]);
		else if (a[h + 1] > a[t + 1])
			std::swap(a[h + 1], a[t + 1]);
		assert(t < n);
#endif
	}
}

#if __linux__
#include <sys/resource.h>
#elif _WIN32
//#include <windows.h>
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
	return clock() / 1000;
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
	rand_swap(a, n, k);
	auto ts = getTime();

	std::nth_element(a, a + k, a + n);
	int maxe = *std::max_element(a, a + k);
	int sum =  accumulate(a, a + k, 0);
	printf("stl nth_element  %4ld ms, a[%d] = %d, sum = %d\n", getTime() - ts, k, maxe, sum);
}

void stl_makeheap(int a[], int n, const int k)
{
	reset(a, n);
	rand_swap(a, n, k);
	auto ts = getTime();

	std::make_heap(a, a + k);

	int maxe = a[0];
	for (int i = k; i < n; ++i) {
		if (a[i] < maxe) {
			a[k] = a[i];
			push_heap(a, a + k + 1);
			pop_heap(a,  a + k + 1);
			maxe = a[0];
		}
	}

	int sum = accumulate(a, a + k, 0);
	printf("stl make_heap   %4ld ms, a[%d] = %d, sum = %d\n", getTime() - ts, k, maxe, sum);
}

void stl_priqueue(int a[], int n, const int k)
{
	reset(a, n);
	rand_swap(a, n, k);
	auto ts = getTime();

	priority_queue<int> pri_queue;
	for (int i = 0; i < k; i++) {
		pri_queue.push(a[i]);
	}

	int maxe = pri_queue.top();
	for (int i = k; i < n; i++) {
		if (a[i] < maxe) {
			pri_queue.push(a[i]);
			pri_queue.pop();
			maxe = pri_queue.top();
		}
	}

	printf("stl priority_queue %ld ms, a[%d] = %d\n", getTime() - ts, k, maxe);
}

void max_heap(int a[], int n, const int k)
{
	reset(a, n);
	rand_swap(a, n, k);
	auto ts = getTime();

	maxheap<int> my_heap;
	for (int i = 0; i < k; i++) {
		my_heap.push(a[i]);
	}

	int maxe = my_heap.top();
	for (int i = k; i < n; i++) {
		if (a[i] < maxe) {
			my_heap.push(a[i]);
			my_heap.pop();
			maxe = my_heap.top();
		}
	}

	int sum = accumulate(my_heap.a + 1, my_heap.a + k + 1, 0);
	printf("my max_heap     %4ld ms, a[%d] = %d, sum = %d\n", getTime() - ts, k, maxe, sum);
}

void radix_sort(int a[], int n, const int k)
{
	reset(a, n);
	rand_swap(a, n, k);
	auto ts = getTime();

	const unsigned int segment = 18;
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
				j = k;
				maxe = a[j - 1];
			}
		}
	}

	std::sort(a, a + j);
	int sum = accumulate(a, a + k, 0);
	printf("radix_sort       %4ld ms, a[%d] = %d, sum = %d\n", getTime() - ts, k, a[k - 1], sum);
}

void merge_sort(int a[], int n, const int k)
{
	reset(a, n);
	rand_swap(a, n, k);

	auto ts = getTime();

	std::sort(a, a + k);

	int* best_a = a + k;
	int maxe = a[k - 1];
	int bestn = 0;
	for (int i = k; i < n; i++) {
		if (a[i] < maxe) {
			best_a[bestn++] = a[i];
			if (bestn == k / 16) {
				std::sort(best_a, best_a + bestn);
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

	int i = k * 6 / 10, s = 64;
	int n = k - 1 + s;
	while (i > s && a[i - s] > b[n - i]) {
		i -= s;
	}

	for (; i >= 0; i--) {
		int j = k - 1 - i;
		if (a[i] <= b[j]) {
#ifndef ME
			int m = k - 1;
			while (j > 0) {
				if (a[i] > b[j])
					a[m --] = a[i--];
				else
					a[m --] = b[--j];
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
	rand_swap(a, n, k);

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

static int arr[MAXN];
int main(int argc, char* argv[])
{
	int maxn = MAXN, k = MAXK, type = 0;
	srand(time(NULL));

	printf("cmd:topk type(0(rand),1(decrease),2(increase),3(wavy)) k(<=%d) n(<=%d)\n\n", MAXK, MAXN);

	if (argc > 1) { type = atoi(argv[1]); }
	if (argc > 2) { k  = atoi(argv[2]);   if (k <= 0  || k > maxn || k > MAXK) k = MAXK;   }
	if (argc > 3) { maxn = atoi(argv[3]); if (maxn <= 0 || maxn > MAXN) maxn = MAXN; }
	assert(k  <= MAXK && k < MAXN / 2);
	assert(maxn <= MAXN);

	int s = rand();
	for (int j = 2; j > 0; j --) {
		for (int i = 0; i < maxn; i++) {
			int r = rand() * rand() + rand() + s;
			if (type == 1)
				r = s + i;
			else if (type == 2)
				r = maxn - i - s;
			else if (type == 3)
				r = (s + i) * (i + j);
			if (r < 0)
				r = 0 - r;
			buff[i] = r;
		}

		printf("maxn = %d, topk = %d, type = %d\n\n", maxn, k, type);
		for (int i = 0; i < maxn; i++) {
			if (buff[i] < 0) {
				buff[i] = 0 - buff[i];
			}
		}

		stl_nth(arr, maxn, k);
		radix_sort(arr, maxn, k);
		stl_priqueue(arr, maxn, k);
		max_heap(arr, maxn, k);
		stl_makeheap(arr, maxn, k);
		merge_sort2(arr, maxn, k);
		merge_sort(arr, maxn, k);

		putchar('\n');
	}

	return 0;
}

