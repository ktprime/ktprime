
#pragma once

#include <vector>

template<typename T, int D = 4, typename COMP = std::less<T>>
class d_ary
{
public:
    explicit d_ary(int capacity=4) { heap.reserve(capacity); }
    bool empty() const { return heap.size() == 0; }
    size_t size() const { return  heap.size(); }
    void clear() { heap.clear(); }

    T& top()  { return heap.front(); }
    const T& top() const { return heap.front(); }

    T& last() { return heap.back(); }
    const T& last() const { return heap.back(); }

    void push(T&& i)
    {
        heap.push_back(std::move(i));
        fixUp(heap.size() - 1);
    }

    void push(const T& i)
    {
        heap.push_back(i);
        fixUp(heap.size() - 1);
    }

    void pop()
    {
        heap[0] = heap.back();
        heap.pop_back();
        fixDown(0);
    }

    void pop_push(const T& i)
    {
        heap[0] = i;
        fixDown(0);
    }

private:
    std::vector<T> heap;
    COMP compare;

    int parent(int i) { return (i - 1) / D; };
    int child(int i, int j) {
        int o = D * i + j;
        return o < size() ? o : -1;
    }

    int smallestChild(int c)
    {
        int min = D * c + 1;
        if (min + D - 1 < size()) {
            if (D == 4) {
                const int min1 = compare(heap[min + 0], heap[min + 1]) ? min + 1 : min + 0;
                const int min2 = compare(heap[min + 2], heap[min + 3]) ? min + 3 : min + 2;
                return compare(heap[min1], heap[min2]) ? min2 : min1;
            } else if (D == 8) {
                const int min1 = compare(heap[min + 0], heap[min + 1]) ? min + 1 : min + 0;
                const int min2 = compare(heap[min + 2], heap[min + 3]) ? min + 3 : min + 2;
                const int min3 = compare(heap[min + 4], heap[min + 5]) ? min + 5 : min + 4;
                const int min4 = compare(heap[min + 6], heap[min + 7]) ? min + 7 : min + 6;

                const int min5 = compare(heap[min1], heap[min2]) ? min2 : min1;
                const int min6 = compare(heap[min3], heap[min4]) ? min4 : min3;
                return compare(heap[min5], heap[min6]) ? min6 : min5;
            } else if (D == 6) {
                const int min1 = compare(heap[min + 0], heap[min + 1]) ? min + 1 : min + 0;
                const int min2 = compare(heap[min + 2], heap[min + 3]) ? min + 3 : min + 2;
                const int min3 = compare(heap[min + 4], heap[min + 5]) ? min + 5 : min + 4;
                const int min4 = compare(heap[min1 + 0], heap[min2 + 0]) ? min2 : min1;
                return compare(heap[min3], heap[min4]) ? min4 : min3;
            }
        }

        //SIMD
        const auto chr = std::min((int)size(), min + D);
        for (auto chi = min + 1; chi < chr; ++chi) {
            if (compare(heap[min], heap[chi]))
                min = chi;
        }
        return min;
    }

    void fixUp(int c)
    {
        auto v = heap[c];
        auto p = (c - 1) / D;

        while (compare(heap[p], v)) {
            heap[c] = std::move(heap[p]);
            c = p;
            if (p == 0)
                break;
            p = (c - 1) / D;
        }
        heap[c] = std::move(v);
    }

    void fixDown(int p)
    {
        int min = D * p + 1;
        if (min >= size())
            return ;

        int chimin = smallestChild(p);
        auto v = heap[p];
        while (compare(v, heap[chimin])) {
            heap[p] = heap[chimin];
            p = chimin;
            chimin = smallestChild(chimin);
            if (chimin >= size())
                break;
        }
        heap[p] = v;
    }
};


