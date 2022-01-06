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
#include "d-ary.h"

#define NC   1
//#define PDQS 1
//#define SKA_SORT 1
//#define MIN_MAX_HEAP 1

#ifdef __has_include
    #if __has_include("pdqsort.h")
        #include "pdqsort.h" //https://github.com/orlp/pdqsort/blob/master/pdqsort.h
        #define PDQS 1
    #endif
    #if __has_include("minmax_and_dary_heap.hpp")
        #include "minmax_and_dary_heap.hpp"
        #define MIN_MAX_HEAP 1
    #endif
#else

#if MIN_MAX_HEAP
#include "minmax_and_dary_heap.hpp"
#endif

#if PDQS
#include "pdqsort.h" //https://github.com/orlp/pdqsort/blob/master/pdqsort.h
#endif

#if SKA_SORT
#include "ska_sort.hpp" //https://github.com/skarupke/ska_sort/blob/master/ska_sort.hpp
#endif

#endif

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

// likely/unlikely
#if (__GNUC__ >= 4 || __clang__)
#    define EHASH_LIKELY(condition)   __builtin_expect(condition, 1)
#    define EHASH_UNLIKELY(condition) __builtin_expect(condition, 0)
#else
#    define EHASH_LIKELY(condition)   condition
#    define EHASH_UNLIKELY(condition) condition
#endif

static int AKI = 4;
static int AK2 = 4;
static int pdqs = 2;

static int64_t gsum = 0;
static stype  gmaxe = 0;

typedef unsigned int uint;

template<typename T>
class max_heap
{
public:
    max_heap(int size)
    {
        _size = 0;
        _a = new T[size * 2 + 2];

#ifndef NC
        for (int i = 3; i <= size + 2; i++)
            _a[size + i] = ~max_v;
#endif
        _a[size + 1] = _a[size + 2] = ~max_v;
        _a[0] = max_v;
    }

    void make_heap(const T* a, int size)
    {
#if 0
        memcpy(_a + 1, a, sizeof(a[0]) * size);
        _size = size;
        std::make_heap(_a + 1, _a + _size + 1);
#else
        for (int i = 0; i < size; i++)
            push(a[i]);
#endif
    }

    ~max_heap()   { delete []_a; }
    T top() const { return _a[1];}

    void push(const T& v)
    {
        uint c = ++_size;
        uint p = _size / 2;

        while (v > _a[p] /*&& p >= 1*/) {
            _a[c] = _a[p];
            c = p;
            p /= 2;
        }
        _a[c] = v;
    }

    void pop()
    {
        const T& v = _a[_size--];
        uint p = 1;
        uint c = _size > 2 ? (_a[3] > _a[2]) + 2 : _size;

        while (v < _a[c]) {
            _a[p] = _a[c];
            p = c;
            c *= 2;
#ifdef NC
            if (c > _size)
                break;
#endif
            c += _a[c + 1] > _a[c];
        }
        _a[p] = v;
    }

    //https://en.wikipedia.org/wiki/Heapsort
    //https://stackoverflow.com/questions/39095000/if-siftdown-is-better-than-siftup-why-do-we-have-it
    T pop_push(const T v)
    {
        uint p = 1, l = 2, r = 3;

        while (true) {
            const uint c = _a[l] >= _a[r] ? l : r;
            if (v >= _a[c])
                break;

            _a[p] = _a[c];
            p = c;
            l = c * 2 + 0;
#ifdef NC
            if (l > _size)
                break;
#endif
            r = c * 2 + 1;
        }

        _a[p] = v;
        return _a[1];
    }

    T pop_push2(const T v)
    {
        _a[++_size] = v;
        pop();
        return _a[1];
    }
    //private:

    T *_a;
    uint _size;
};

template<typename T>
inline void QSORT(T* a, T* e)
{
#if PDQS
    if (pdqs == 2)
        pdqsort_branchless(a, e);
#if SKA_SORT
    else if (pdqs == 1)
        ska_sort(a, e);
#endif
    else
        std::sort(a, e);
#elif SKA_SORT
    if (pdqs > 0)
        ska_sort(a, e);
    else
        std::sort(a, e);
#else
    std::sort(a, e);
#endif
}

static void rand_swap(stype a[], const int n, int k)
{
#if 0
    const int step = n / k;
    //    QSORT(a, a + k);
    for (int i = 1; i < k; i ++) {

        int h = rand() % k, t = i * step + rand() % step;
        if (a[h] > a[t])
            std::swap(a[h], a[t]);
        else if (a[i] > a[t])
            std::swap(a[i], a[t]);
    }
#endif
}

static void reset_ank(stype a[], const int n, int k)
{
    memcpy(a, a + n, k * sizeof(a[0]));
    rand_swap(a, n, 10000);
}

static void check_sum(int64_t sum, stype maxe, const char* func)
{
    if (sum != gsum || maxe != gmaxe) {
        if (gsum == 0) {
            gsum = sum;
            gmaxe = maxe;
        }
        else
            printf("%s %ld %ld\n", func, (long)sum, (long)maxe);
    }
}

static void check_ak(const stype a[], const int n, int k)
{
    if (a[2*n] == (stype)(-1)) {
        return;
    }

    for (int i = 0; i < k; i++) {
        if (a[i] != a[n * 2 + i]) {
            printf("%d %ld != %ld\n", i, (long)a[i], (long)a[n * 2 + i]);
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
#if __linux__ || __unix__
    struct rusage rup;
    getrusage(RUSAGE_SELF, &rup);
    long sec = rup.ru_utime.tv_sec + rup.ru_stime.tv_sec;
    long usec = rup.ru_utime.tv_usec + rup.ru_stime.tv_usec;
    return sec * 1000 + usec / 1000;
#elif _WIN32
    LARGE_INTEGER freq = {0};
//    if (freq.QuadPart == 0)
    QueryPerformanceFrequency(&freq);

    LARGE_INTEGER nowus = {0};
    QueryPerformanceCounter(&nowus);
    return (nowus.QuadPart * 1000) / (freq.QuadPart);
#else
    return clock();
#endif
}

#if __cplusplus
#if 0
static void stl_sort(stype a[], const int n, const int k)
{
    reset_ank(a, n, k);
    clock_t ts = NowMs();

    QSORT(a, a + n);
    printf("  stl sort       %5ld ms, a[%d] = %ld\n\n", NowMs() - ts, k, (int64_t)a[k - 1]);
}
#endif

static void stl_nth(stype a[], const int n, const int k)
{
    reset_ank(a, n, n);
    clock_t ts = NowMs();

#if 1
    std::nth_element(a, a + k, a + n);
    QSORT(a, a + k);
#else
    std::partial_sort(a, a + k, a + n);
#endif

    stype maxe = a[k - 1];
    printf("  stl nth_element %5ld ms\n", NowMs() - ts);

    memcpy(a + n * 2, a, k * sizeof(a[0]));
    check_ak(a, n, k);
    int64_t sum = std::accumulate(a, a + k, 0);
    check_sum(sum, maxe, __FUNCTION__);
}

static void stl_priqueue(stype a[], const int n, const int k)
{
//    reset_ank(a, n, k);
    clock_t ts = NowMs();

    std::priority_queue<stype> pri_queue;
    for (int m = 0; m < k; m++)
        pri_queue.push(a[m + n]);

    stype maxe = pri_queue.top();
    for (int i = n + k; i < n * 2; i++) {
        if (a[i] < maxe) {
            pri_queue.pop();
            pri_queue.push(a[i]);
            maxe = pri_queue.top();
        }
    }

    ts = NowMs() - ts;
    for (int j = 1; j <= k; j ++) {
        a[k - j] = pri_queue.top();
        pri_queue.pop();
    }

    printf("  stl pri_queue   %5ld ms\n", ts);

    check_ak(a, n, k);
    int64_t sum = std::accumulate(a, a + k, 0);
    check_sum(sum, maxe, __FUNCTION__);
}

static void kary_heap(stype a[], const int n, const int k)
{
    clock_t ts = NowMs();
    d_ary<stype, 4> kary_heap;
    for (int m = 0; m < k; m++)
        kary_heap.push(a[m + n]);

    stype maxe = kary_heap.top();
    for (int i = n + k; i < n * 2; i++) {
        if (a[i] < maxe) {
            kary_heap.pop();
            kary_heap.push(a[i]);
            maxe = kary_heap.top();
        }
    }

    ts = NowMs() - ts;
    for (int j = 1; j <= k; j ++) {
        a[k - j] = kary_heap.top();
        kary_heap.pop();
    }

    printf("  kary4 heap      %5ld ms\n", ts);

    check_ak(a, n, k);
    int64_t sum = std::accumulate(a, a + k, 0);
    check_sum(sum, maxe, __FUNCTION__);
}

#if MIN_MAX_HEAP
static void min_max_heap(stype a[], const int n, const int k)
{
    reset_ank(a, n, k);
    clock_t ts = NowMs();

    make_minmax_heap(a, a + k, std::greater<stype>());
    auto maxe = a[0];

    for (int i = n + k; i < n * 2; i++) {
        if (a[i] < maxe) {
            pop_minmax_heap_min(a, a + k, std::greater<stype>());
            a[k - 1] = a[i];
            push_minmax_heap(a, a + k, std::greater<stype>());
            maxe = a[0];
        }
    }

    printf("  min_max_heap    %5ld ms\n", NowMs() - ts);
    std::sort(a, a + k);

    check_ak(a, n, k);
    int64_t sum = std::accumulate(a, a + k, 0);
    check_sum(sum, maxe, __FUNCTION__);
}

static void dary_heap(stype a[], const int n, const int k)
{
    reset_ank(a, n, k);
    clock_t ts = NowMs();

    make_dary_heap<4>(a, a + k);
    auto maxe = a[0];

    for (int i = n + k; i < n * 2; i++) {
        if (a[i] < maxe) {
            pop_dary_heap<4>(a, a + k);
            a[k - 1] = a[i];
            push_dary_heap<4>(a, a + k);
            maxe = a[0];
        }
    }

    printf("  dary_heap4      %5ld ms\n", NowMs() - ts);
    std::sort(a, a + k);

    check_ak(a, n, k);
    int64_t sum = std::accumulate(a, a + k, 0);
    check_sum(sum, maxe, __FUNCTION__);
}
#endif

static void stl_makeheap(stype a[], const int n, const int k)
{
    reset_ank(a, n, k);
    clock_t ts = NowMs();

    std::make_heap(a, a + k);
    stype maxe = a[0];

#if 0
    uint lmax = 1, rmax = 2;
    while (lmax < k && rmax < k) {
        auto ll = lmax * 2 + 1;
        if (ll <= k)
            lmax = ll;
        else if (ll - 1 <= k)
            lmax = ll - 1;
        else
            lmax = k;

        auto rr = rmax * 2 + 1;
        if (rr <= k)
            rmax = rr;
        else if (rr - 1 <= k)
            rmax = rr - 1;
        else
            rmax = k;
    }
#endif

    for (int i = n + k; i < n * 2; i++) {
        if (a[i] < maxe) {
#ifndef PPOP_PUSH
            std::pop_heap(a, a + k);
            a[k - 1] = a[i];
            std::push_heap(a, a + k);
            maxe = a[0];
#elif 1
            //a[0] = a[i]; std::pop_push(a, a + k); maxe = a[0];
            a[k] = a[i]; std::pop_heap(a, a + k + 1); maxe = a[0];
#else
            if (a[1] >= a[2]) {
                maxe = a[0] = a[1];
                a[1] = a[i];
                std::sift_down(a + 1, a + lmax);
            } else {
                maxe = a[0] = a[2];
                a[2] = a[i];
                std::sift_down(a + 2, a + rmax);
            }
#endif
        }
    }

    printf("  stl make_heap   %5ld ms\n", NowMs() - ts);
    std::sort_heap(a, a + k);

    check_ak(a, n, k);
    int64_t sum = std::accumulate(a, a + k, 0);
    check_sum(sum, maxe, __FUNCTION__);
}

static void ktprime_heap(stype a[], const int n, const int k)
{
    clock_t ts = NowMs();
    max_heap<stype> my_heap(k);

    my_heap.make_heap(a + n, k);
//    for (int i = 0; i < k; i++)
//        my_heap.push(a[n + i]);

    stype maxe = my_heap.top();
    for (int i = n + k; i < n * 2; i++) {
        if (a[i] < maxe) {
#if POP_PUSH
            my_heap.pop(); my_heap.push(a[i]);
            maxe = my_heap.top();
#else
            maxe = my_heap.pop_push2(a[i]);
#endif
        }
    }

    memcpy(a, my_heap._a + 1, k * sizeof(a[0]));

    printf("  ktprime heap    %5ld ms\n", NowMs() - ts);
    std::sort_heap(a, a + k);

    check_ak(a, n, k);
    int64_t sum = std::accumulate(a, a + k, 0);
    check_sum(sum, maxe, __FUNCTION__);
}
#endif

static void bucket_sort(uint a[], uint n, uint k)
{
    reset_ank((stype*)a, n, n);
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
            continue;

        a[sum_size++] = a[i];
        if (sum_size > min_size) {
            QSORT(a, a + sum_size);
            maxe = a[(sum_size = k) - 1];
        }
    }

    QSORT(a, a + sum_size);

#if I32
    for (int i = 0; i < sum_size; i++) a[i] -= INT_MAX;
#endif


    printf("  bucket_sort     %5ld ms\n", NowMs() - ts);
    check_ak((stype*)a, n, k);
    int64_t sum = std::accumulate(a, a + k, 0);
    check_sum(sum, maxe, __FUNCTION__);
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

//    assert(i + j == kth);

//    while (a[i - 1] > b[j]) i--, j++;
//    while (a[i + 1] < b[j]) i ++, j--;

    return (j >= n || (i < m && a[i] <= b[j])) ? i : i - 1;
}

static void merge_array(stype a[], stype b[], const int n, const int m)
{
    QSORT(b, b + m);
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
    //TODO: check_ak i overflow
    while (j >= 0) {
        a[t --] = b[j] >= a[i] ? b[j --] : a[i --];
    }
}

const int l2_cpu_size = 1024 * 1024, l1_cpu_size = 32 * 1024;
static void merge_sort(stype a[], const int n, const int k)
{
    reset_ank(a, n, k);
    clock_t ts = NowMs();

    QSORT(a, a + k);
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
            continue;
        }

        ax_a[ax_n++] = a[i];
        if (ax_n == bestn) {
            merge_array(a, ax_a, k, ax_n);
            maxe = a[k - 1];
            ax_n = 0;
        }
    }

//    QSORT(a, ax_a + ax_n);
//    QSORT(ax_a, ax_a + ax_n); std::inplace_merge(a, ax_a, ax_a + ax_n);
//    std::partial_sort(a, ax_a, ax_a + ax_n);
    merge_array(a, ax_a, k, ax_n);

    printf("  merge_sort      %5ld ms\n", NowMs() - ts);
    check_ak(a, n, k);
    int64_t sum = std::accumulate(a, a + k, 0);
    maxe = a[k - 1];
    check_sum(sum, maxe, __FUNCTION__);
}

static void merge_inplace(stype a[], const int n, const int k)
{
    reset_ank(a, n, k);
    clock_t ts = NowMs();

    QSORT(a, a + k);
    stype* ax_a = a + k;
    stype maxe = a[k - 1];
    int bestn = k / AKI, ax_n = 0;

    if (bestn > l2_cpu_size / sizeof(a[0]))
        bestn = l2_cpu_size / sizeof(a[0]);
    if (bestn < 100)
        bestn = k * 4 + 100;

    for (int i = n + k; i < n * 2; i++) {
        if (a[i] >= maxe) {
            continue;
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
            QSORT(ax_a, ax_a + ax_n);

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

    QSORT(a, ax_a + ax_n);
//    QSORT(ax_a, ax_a + ax_n); std::inplace_merge(a, ax_a, ax_a + ax_n);

    printf("  merge_inplace   %5ld ms\n", NowMs() - ts);

    check_ak(a, n, k);
    int64_t sum = std::accumulate(a, a + k, 0);
    maxe = a[k - 1];
    check_sum(sum, maxe, __FUNCTION__);
}

static void printInfo()
{
    const char* sepator =
        "------------------------------------------------------------------------------------------------------------";
    puts(sepator);
    //    puts("Copyright (C) by 2018-2021 Huang Yuanbing bailuzhou at 163.com\n");

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
#elif __arm64__ || __aarch64__
    info += sprintf(info, " arm64");
#elif __arm__
    info += sprintf(info, " arm");
#else
    info += sprintf(info, " unknow");
#endif

    puts(cbuff);
    puts(sepator);
}

static void ktprime_pushpop2(stype a[], const int n, const int k)
{
    max_heap<stype> my_heap(k);

    for (int i = 0; i < k; i++)
        my_heap.push(a[n + i]);

    clock_t ts = NowMs();
    int64_t sum = 0;
    for (int i = n + k; i < n * 2; i++) {
        sum += my_heap.top() * i;
        my_heap.pop_push2(a[i]);
    }

    stype maxe = my_heap.top();
    printf("\t%s  %5ld ms\n", __FUNCTION__,  NowMs() - ts);
    check_sum(sum, maxe, __FUNCTION__);
}

static void ktprime_pushpop(stype a[], const int n, const int k)
{
    max_heap<stype> my_heap(k);

    for (int i = 0; i < k; i++)
        my_heap.push(a[n + i]);

    clock_t ts = NowMs();
    int64_t sum = 0;
    for (int i = n + k; i < n * 2; i++) {
        sum += my_heap.top() * i;
        my_heap.pop_push(a[i]);
    }

    stype maxe = my_heap.top();
    printf("\t%s   %5ld ms\n", __FUNCTION__,  NowMs() - ts);
    check_sum(sum, maxe, __FUNCTION__);
}

static void ktprime_heap3(stype a[], const int n, const int k)
{
    clock_t ts = NowMs();
    max_heap<stype> my_heap1(k);
    max_heap<stype> my_heap2(k);

    for (int i = 0; i < k; i++) {
        my_heap1.push(a[n + i]);
        my_heap2.push(a[n + i]);
    }

    int64_t sum = 0;
    for (int i = n + k; i < n * 2; i++) {
        assert(my_heap1.top() == my_heap2.top());

        sum += my_heap1.top();
        my_heap1.pop(), my_heap1.push(a[i]);
        my_heap2.pop_push(a[i]);
        if (i % (1024 * 32) == 0) {
            assert( std::accumulate(my_heap1._a + 1, my_heap1._a + k + 1, 0) ==
                    std::accumulate(my_heap2._a + 1, my_heap2._a + k + 1, 0));
        }
    }

    stype maxe = my_heap1.top();
    printf("\tktprime heap     %5ld ms\n", NowMs() - ts);
    check_sum(sum, maxe, __FUNCTION__);
}

static void kary_pushpop(stype a[], const int n, const int k)
{
    d_ary<stype, 4> kary_heap;
    for (int m = 0; m < k; m++)
        kary_heap.push(a[m + n]);

    clock_t ts = NowMs();
    int64_t sum = 0;
    for (int i = n + k; i < n * 2; i++) {
        sum += kary_heap.top() * i;
        kary_heap.pop_push(a[i]);
    }
    stype maxe = kary_heap.top();

    printf("\tkary4_pushpop     %5ld ms\n", NowMs() - ts);
    check_sum(sum, maxe, __FUNCTION__);
}

#if MIN_MAX_HEAP
static void dary_pushpop(stype a[], const int n, const int k)
{
    reset_ank(a, n, k);
    make_dary_heap<4>(a, a + k);

    clock_t ts = NowMs();
    int64_t sum = 0;
    for (int i = n + k; i < n * 2; i++) {
        sum += a[0] * i;
        a[k] = a[i];
        pop_dary_heap<4>(a, a + k + 1);
    }
    stype maxe = a[0];
    printf("\tdary4_pushpop     %5ld ms\n", NowMs() - ts);
    check_sum(sum, maxe, __FUNCTION__);
}

static void minmax_pushpop(stype a[], const int n, const int k)
{
    reset_ank(a, n, k);
    make_minmax_heap(a, a + k, std::greater<stype>());

    clock_t ts = NowMs();
    int64_t sum = 0;
    for (int i = n + k; i < n * 2; i++) {
        sum += a[0] * i;
        a[k] = a[i];
        pop_minmax_heap_min(a, a + k + 1, std::greater<stype>());
    }
    stype maxe = a[0];
    printf("\t%s    %5ld ms\n", __FUNCTION__, NowMs() - ts);
    check_sum(sum, maxe, __FUNCTION__);
}
#endif

static void stlpri_pushpop(stype a[], const int n, const int k)
{
    //    reset_ank(a, n, k);
    std::priority_queue<stype> pri_queue;
    for (int m = 0; m < k; m++)
        pri_queue.push(a[m + n]);

    clock_t ts = NowMs();
    int64_t sum = 0;
    for (int i = n + k; i < n * 2; i++) {
        sum += pri_queue.top() * i;
        pri_queue.pop();
        pri_queue.push(a[i]);
    }
    stype maxe = pri_queue.top();

    printf("\t%s    %5ld ms\n", __FUNCTION__, NowMs() - ts);
    check_sum(sum, maxe, __FUNCTION__);
}

static void stlheap_pushpop2(stype a[], const int n, const int k)
{
    reset_ank(a, n, k);
    std::make_heap(a, a + k);

    clock_t ts = NowMs();
    int64_t sum = 0;
    for (int i = n + k; i < n * 2; i++) {
        sum += a[0] * i;
        std::pop_heap(a, a + k);
        a[k - 1] = a[i];
        std::push_heap(a, a + k);
    }
    stype maxe = a[0];

    printf("\t%s   %5ld ms\n", __FUNCTION__, NowMs() - ts);
    check_sum(sum, maxe, __FUNCTION__);
}

static void stlheap_pushpop(stype a[], const int n, const int k)
{
    reset_ank(a, n, k);
    std::make_heap(a, a + k);

    clock_t ts = NowMs();
    int64_t sum = 0;
    for (int i = n + k; i < n * 2; i++) {
        sum += a[0] * i;
#ifndef PPOP_PUSH
        a[k] = a[i]; std::pop_heap(a, a + k + 1);
#else
        a[0] = a[i]; std::pop_push(a, a + k);
#endif
    }
    stype maxe = a[0];

    printf("\t%s    %5ld ms\n", __FUNCTION__, NowMs() - ts);
    check_sum(sum, maxe, __FUNCTION__);
}

void heap_push_pop(stype a[], const int n, const int k)
{
    clock_t ts = 0;
    gsum = gmaxe = 0;
    const auto k2 = k * 2;

    {
        ts = NowMs();long mask = 0;
        std::priority_queue<stype> pri_queue;
        for (int m = 0; m < k2; m++)
            pri_queue.push(a[m + n]);
        printf("\t pri_queue    %5ld ms, ", NowMs() - ts);
        ts = NowMs();

        for (int m = 0; m < k; m++) {
            mask += pri_queue.top() * m;
            pri_queue.pop();
        }
        printf("pop %5ld ms\n", NowMs() - ts);
        check_sum(mask, 0, __FUNCTION__);
    }

    {
        ts = NowMs();long mask = 0;
        d_ary<stype, 4> kary_heap;
        for (int m = 0; m < k2; m++)
            kary_heap.push(a[m + n]);
        printf("\t kary4        %5ld ms, ", NowMs() - ts);
        ts = NowMs();
        for (int m = 0; m < k; m++) {
            mask += kary_heap.top() * m;
            kary_heap.pop();
        }
        printf("pop %5ld ms\n", NowMs() - ts);
        check_sum(mask, 0, __FUNCTION__);
    }

    if (0)
    {
        ts = NowMs();long mask = 0;
        d_ary<stype, 6> kary_heap;
        for (int m = 0; m < k2; m++)
            kary_heap.push(a[m + n]);
        printf("\t karyp6        %5ld ms, ", NowMs() - ts);
        ts = NowMs();
        for (int m = 0; m < k; m++) {
            mask += kary_heap.top() * m;
            kary_heap.pop();
        }
        printf("pop %5ld ms\n", NowMs() - ts);
        check_sum(mask, 0, __FUNCTION__);
    }

    {
        ts = NowMs();long mask = 0;
        d_ary<stype, 8> kary_heap;
        for (int m = 0; m < k2; m++)
            kary_heap.push(a[m + n]);
        printf("\t kary8        %5ld ms, ", NowMs() - ts);
        ts = NowMs();

        for (int m = 0; m < k; m++) {
            mask += kary_heap.top() * m;
            kary_heap.pop();
        }
        printf("pop %5ld ms\n", NowMs() - ts);
        check_sum(mask, 0, __FUNCTION__);
    }

    {
        ts = NowMs();long mask = 0;
        max_heap<stype> my_heap(k);
        for (int m = 0; m < k2; m++)
            my_heap.push(a[m + n]);
        printf("\t my_heap      %5ld ms, ", NowMs() - ts);
        ts = NowMs();

        for (int m = 0; m < k; m++) {
            mask += my_heap.top() * m;
            my_heap.pop();
        }
        printf("pop %5ld ms\n", NowMs() - ts);
        check_sum(mask, 0, __FUNCTION__);
    }


    {
        ts = NowMs();long mask = 0;
        for (int m = 0; m < k2; m++) {
            a[m] = a[m + n];
            std::push_heap(a, a + m);
        }
        printf("\t make_heap    %5ld ms, ", NowMs() - ts);
        ts = NowMs();

        for (int m = k2; m > k; m--) {
            mask += a[0] * (k2 - m);
            std::pop_heap(a, a + m);
        }
        printf("pop %5ld ms\n", NowMs() - ts);
        check_sum(mask, 0, __FUNCTION__);
    }

#if MIN_MAX_HEAP
    {
        ts = NowMs();long mask = 0;
        for (int m = 0; m < k2; m++) {
            a[m] = a[m + n];
            push_minmax_heap(a, a + m, std::greater<stype>());
        }
        printf("\t min_max      %5ld ms, ", NowMs() - ts);
        ts = NowMs();

        for (int m = k2; m > k; m--) {
            mask += a[0] * (k2 - m);
            pop_minmax_heap_min(a, a + m, std::greater<stype>());
        }
        printf("pop %5ld ms\n", NowMs() - ts);
        check_sum(mask, 0, __FUNCTION__);
    }

    if (0)
    {
        ts = NowMs();long mask = 0;
        for (int m = 0; m < k2; m++) {
            a[m] = a[m + n];
            push_dary_heap<5>(a, a + m);
        }
        printf("\t dary_heap5   %5ld ms, ", NowMs() - ts);
        ts = NowMs();

        for (int m = k2; m > k; m--) {
            mask += a[0] * (k2 - m);
            pop_dary_heap<5>(a, a + m);
        }
        printf("pop %5ld ms\n", NowMs() - ts);
        check_sum(mask, 0, __FUNCTION__);
    }

    {
        ts = NowMs();long mask = 0;
        //    make_dary_heap<8>(a, a + k);
        for (int m = 0; m < k2; m++) {
            a[m] = a[m + n];
            push_dary_heap<4>(a, a + m);
        }
        printf("\t dary_heap4   %5ld ms, ", NowMs() - ts);
        ts = NowMs();

        for (int m = k2; m > k; m--) {
            mask += a[0] * (k2 - m);
            pop_dary_heap<4>(a, a + m);
        }
        printf("pop %5ld ms\n", NowMs() - ts);
        check_sum(mask, 0, __FUNCTION__);
    }

    if (0)
    {
        ts = NowMs();long mask = 0;
        for (int m = 0; m < k2; m++) {
            a[m] = a[m + n];
            push_dary_heap<6>(a, a + m);
        }
        printf("\t dary_heap6   %5ld ms, ", NowMs() - ts);
        ts = NowMs();

        for (int m = k2; m > k; m--) {
            mask += a[0] * (k2 - m);
            pop_dary_heap<6>(a, a + m);
        }
        printf("pop %5ld ms\n", NowMs() - ts);
        check_sum(mask, 0, __FUNCTION__);
    }

    {
        ts = NowMs(); long mask = 0;
        for (int m = 0; m < k2; m++) {
            a[m] = a[m + n];
            push_dary_heap<8>(a, a + m);
        }
        printf("\t dary_heap8   %5ld ms, ", NowMs() - ts);
        ts = NowMs();

        for (int m = k2; m > k; m--) {
            mask += a[0] * (k2 - m);
            pop_dary_heap<8>(a, a + m);
        }
        printf("pop %5ld ms\n", NowMs() - ts);
        check_sum(mask, 0, __FUNCTION__);
    }
#endif
}

int main(int argc, char* argv[])
{
    srand(time(nullptr));
    const int max_n = 1024*1024*64;

    printf("\ntopk -k(<=%d) -n(<=%d) -r(1-16) m(1-16) p(pdqs 2 ska_sort 1) s(shuffle) q(push_pop)\n\
            [type = 0-1 rand, 2-3 wavy, 4-5 rand, 6 decrease, 7 incre]\n\n", 1000000, max_n);
    printInfo();

    int n = max_n, k = max_n / 1000, type = 0, shuff = 0;
    int push_pop = 0;
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
            AKI = r;
        } else if (c == 'm' && r > 0) {
            AK2 = r;
        } else if (c == 's') {
            shuff = r;
        } else if (c == 'p') {
#if PDQS || SKA_SORT
            pdqs = r;
#endif
        } else if (c == 'q')
            push_pop = r;
    }

    if (k > n)
        k = n;

    auto arr = static_cast<stype *>(malloc(sizeof(stype) * (n * 2 + k + 8))) + 8;
    auto buff = arr + n;
    assert((size_t)buff % 8 == 0);
    buff[n] = (stype)(-1);

    std::random_device rd;
    std::default_random_engine e(rd());
    std::mt19937_64 rng(rd());

    std::uniform_int_distribution<int64_t> u(0, INT64_MAX / 2);
    std::normal_distribution<> d(1 << (20 + rand() % 10), 1 << 16);
//    std::exponential_distribution<> p(0.1);

    printf("n = %d, k = %d, r1 = %d, r2 = %d, shuff = %d, pdqs = %d, push_pop = %d\n", n, k, AKI, AK2, shuff, pdqs, push_pop);

    for (int j = 0; j <= 8; j ++) {
        int64_t r = 0;
        type = j;
        for (size_t i = 1; i <= (size_t)n; i++) {
            if (type == 0) {
                r = rd() + e();
            } else if (type == 2) {
                r = rng();
            } else if (type == 1) {
                r = rng() + i;
            } else if (type == 3) {
                r = d(rng);
            } else if (type == 4) {
                r = i * i + (512 - rng() % 1024) * (i % 8);
            } else if (type == 5) {
                r = rng() - (int)(i * i) % 1024;
            } else if (type == 6) {
                r = n - i - rng() % 8; // i + 1
            } else if (type == 7) {
                r = i + rng() % 1024 - sin(i/360.0) * i;
            } else if (type == 8) {
                r = rng() % (100);
            }
            buff[i] = (stype)r + 1;
        }

        if (shuff) {
            std::shuffle(buff + 1, buff + n + 1 , rng);
        }

        gsum = gmaxe = 0;
        printf("data type = %d\n", type);
        stl_nth(arr, n, k);

        if (sizeof(arr[0]) <= sizeof(int))
            bucket_sort((uint*)arr, n, k);

#if __cplusplus
        stl_priqueue(arr, n, k);
        kary_heap(arr, n, k);
        stl_makeheap(arr, n, k);
        ktprime_heap(arr, n, k);
#endif

#if MIN_MAX_HEAP
        dary_heap(arr, n, k);
//        min_max_heap(arr, n, k);
#endif
//        pdqs = 2;
        merge_inplace(arr, n, k);
        merge_sort(arr, n, k);
//        pdqs = 1; merge_sort(arr, n, k);
//        pdqs = 0; merge_sort(arr, n, k);

        if (push_pop >= 1) {
            //test flow win
            printf("\n\ttest push_pop\n");
            gsum = 0;
            auto k = n / 256;
            stlpri_pushpop(arr, n, k);
            kary_pushpop(arr, n, k);

#if MIN_MAX_HEAP
            dary_pushpop(arr, n, k);
            //minmax_pushpop(arr, n, k);
#endif

            //ktprime_heap3(arr, n, k);
            
            //stlheap_pushpop(arr, n, k);
            //stlheap_pushpop2(arr, n, k);

            ktprime_pushpop(arr, n, k);
            ktprime_pushpop2(arr, n, k);
        }

        if (push_pop > 1) {
            printf("\n\ttest push & pop\n");
            heap_push_pop(arr, n, n / 64);
        }
        putchar('\n');
    }

    free(arr - 8);
    return 0;
}

// https://www.jdoodle.com/online-compiler-c++
// https://repl.it/repls/ProductivePowerlessSlope
// https://www.tutorialspoint.com/online_cpp_compiler.php
// http://rextester.com/l/cpp_online_compiler_gcc
// https://www.onlinegdb.com/
// https://www.ideone.com/3rCdob

