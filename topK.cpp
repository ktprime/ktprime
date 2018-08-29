#include <algorithm>
#include <numeric>
#include <queue>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cassert>

#define  MAXN  10000*10000
#define  MAXK   10000

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
	for (int i = k - 1; i >= 0; i--) {
#if 0
		int h = i, t = n - k + i;
		if (a[h] > a[t])
			std::swap(a[h], a[t]);
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
	return clock();
#endif
}

void stl_sort(int a[], int n, int k)
{
	reset(a, n);
	auto ts = getTime();

	std::sort(a, a + n);
	printf("stl sort       %4ld ms, a[%d] = %d\n\n", getTime() - ts, k, a[k - 1]);
}

void stl_nth(int a[], int n, int k)
{
	reset(a, n);
	auto ts = getTime();

	rand_swap(a, n, k);
	std::nth_element(a, a + k, a + n);
	int maxe = *std::max_element(a, a + k);
	int sum =  accumulate(a, a + k, 0);
	printf("stl nth_element  %4ld ms, a[%d] = %d, sum = %d\n", getTime() - ts, k, maxe, sum);
}

void stl_makeheap(int a[], int n, int k)
{
	reset(a, n);
	auto ts = getTime();

	rand_swap(a, n, k);
	std::make_heap(a, a + k);

	int maxe = a[0];
	for (int i = k; i < n; ++i) {
		if (maxe > a[i]) {
			a[k] = a[i];
			push_heap(a, a + k + 1);
			pop_heap(a,  a + k + 1);
			maxe = a[0];
		}
	}

	int sum = accumulate(a, a + k, 0);
	printf("stl make_heap   %4ld ms, a[%d] = %d, sum = %d\n", getTime() - ts, k, maxe, sum);
}

void stl_priqueue(int a[], int n, int k)
{
	reset(a, n);
	auto ts = getTime();

	rand_swap(a, n, k);
	priority_queue<int> pri_queue;
	for (int i = 0; i < k; i++) {
		pri_queue.push(a[i]);
	}

	int maxe = pri_queue.top();
	for (int i = k; i < n; i++) {
		if (maxe > a[i]) {
			pri_queue.push(a[i]);
			pri_queue.pop();
			maxe = pri_queue.top();
		}
	}

	printf("stl priority_queue %ld ms, a[%d] = %d\n", getTime() - ts, k, maxe);
}

void max_heap(int a[], int n, int k)
{
	reset(a, n);
	auto ts = getTime();

	rand_swap(a, n, k);
	maxheap<int> my_heap;
	for (int i = 0; i < k; i++) {
		my_heap.push(a[i]);
	}

	int maxe = my_heap.top();
	for (int i = k; i < n; i++) {
		if (maxe > a[i]) {
			my_heap.push(a[i]);
			my_heap.pop();
			maxe = my_heap.top();
		}
	}

	int sum = accumulate(my_heap.a + 1, my_heap.a + k + 1, 0);
	printf("my max_heap     %4ld ms, a[%d] = %d, sum = %d\n", getTime() - ts, k, maxe, sum);
}

void radix_sort(int a[], int n, int k)
{
	reset(a, n);
	auto ts = getTime();

	rand_swap(a, n, k);
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
		if (maxe > a[i]) {
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
	printf("radix_sort       %4ld ms, a[%d] = %d, sum = %d\n", getTime() - ts, j, a[k - 1], sum);
}

void merge_sort(int a[], int n, const int k)
{
	reset(a, n);
	auto ts = getTime();

	rand_swap(a, n, k);
	std::sort(a, a + k);

	int* best_a = a + k;
	int maxe = a[k - 1];
	int bestn = 0;
	for (int i = k; i < n; i++) {
		if (maxe > a[i]) {
			best_a[bestn++] = a[i];
			if (bestn >= k / 16) {
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
	printf("sort-merge      %4ld ms,  a[%d] = %d, sum = %d\n", getTime() - ts, k, maxe, sum);
}

void merge_array(int a[], int b[], int n)
{
	std::sort(b, b + n);
	if (a[0] >= b[n - 1]) {
		memcpy(a, b, sizeof(int) * n);
		return;
	}

	for (int i = n - 1; i >= 0; i--) {
		int j = n - 1 - i;
		if (a[i] <= b[j]) {
			memcpy(a + i, b, sizeof(int) * (j + 1));
			/**
			while (a[i - 1] <= a[n - 1]) {
				n--;
			}
			while (a[0] <= b[0]) {
				a++, n--;
			} */
			std::inplace_merge(a, a + i, a + n);
			break;
		}
	}
}

void merge_sort2(int a[], int n, const int k)
{
	reset(a, n);
	auto ts = getTime();

	rand_swap(a, n, k);
	std::sort(a , a + k);

	int* best_a = a + k;
	int bestn = 0;
	int maxe = a[k - 1];
	for (int i = k; i < n; i++) {
		if (maxe > a[i]) {
			best_a[bestn++] = a[i];
			if (bestn >= k) {
				merge_array(a, best_a, bestn);
//				std::sort(a, a + k + bestn);
				maxe = a[k - 1];
				bestn = 0;
			}
		}
	}

	std::sort(a, a + k + bestn);
	maxe = a[k - 1];
	int sum = accumulate(a, a + k, 0);
	printf("sort-merge2     %4ld ms, a[%d] = %d, sum = %d\n", getTime() - ts, k, maxe, sum);
}

static int arr[MAXN];
int main(int argc, char* argv[])
{
	int maxn = MAXN, k = MAXK, type = 0; // 0 1 2 rand + -
	srand(time(NULL));

	if (argc > 1) { type = atoi(argv[1]); }
	if (argc > 2) { maxn = atoi(argv[2]); if (maxn <= 0 || maxn > MAXN) maxn = MAXN; }
	if (argc > 3) { k  = atoi(argv[3]); if (k <= 0  || k > maxn || k > MAXK) k = MAXK;   }
	assert(k  <= MAXK && k < MAXN / 2);
	assert(maxn <= MAXN);

	for (int j = 3; j > 0; j --) {
	for (int i = 0; i < maxn; i++) {
		int r = rand() * rand() + rand() * rand();
		if (type == 1)
			r = i;
		else if (type == 2)
			r = (1 << 30) - i;
		else if (type == 3)
			r = i * i;
		if (r < 0)
			r = 0 - r;
		buff[i] = r;
		//printf("%d %d %d\n", i, r, buff[i]);
	}

	printf("maxn = %d, topk = %d, type = %d\n\n", maxn, k, type);
	for (int i = 0; i < maxn; i++) {
		if (buff[i] < 0) {
			buff[i] = 0 - buff[i];
		}
	}

	//	stl_sort(arr, maxn, k);
		merge_sort2(arr, maxn, k);
		stl_nth(arr, maxn, k);
		radix_sort(arr, maxn, k);
		max_heap(arr, maxn, k);
		stl_priqueue(arr, maxn, k);
		stl_makeheap(arr, maxn, k);
		merge_sort(arr, maxn, k);

		putchar('\n');
	}

	return 0;
}

