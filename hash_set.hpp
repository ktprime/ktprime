
// emhash6::HashSet for C++11
// version 1.2.0
// https://github.com/ktprime/ktprime/blob/master/hash_set.hpp
//
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.
// SPDX-License-Identifier: MIT
// Copyright (c) 2019-2019 Huang Yuanbing & bailuzhou AT 163.com
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE


// From
// NUMBER OF PROBES / LOOKUP       Successful            Unsuccessful
// Quadratic collision resolution   1 - ln(1-L) - L/2    1/(1-L) - L - ln(1-L)
// Linear collision resolution     [1+1/(1-L)]/2         [1+1/(1-L)2]/2
//
// -- enlarge_factor --           0.10  0.50  0.60  0.75  0.80  0.90  0.99
// QUADRATIC COLLISION RES.
//    probes/successful lookup    1.05  1.44  1.62  2.01  2.21  2.85  5.11
//    probes/unsuccessful lookup  1.11  2.19  2.82  4.64  5.81  11.4  103.6
// LINEAR COLLISION RES.
//    probes/successful lookup    1.06  1.5   1.75  2.5   3.0   5.5   50.5
//    probes/unsuccessful lookup  1.12  2.5   3.6   8.5   13.0  50.0

#pragma once

#include <cstring>
#include <cstdlib>
#include <type_traits>
#include <cassert>
#include <utility>
#include <cstdint>
#include <functional>
#include <iterator>

#ifdef  GET_KEY
    #undef  NEXT_BUCKET
    #undef  GET_KEY
    #undef NEW_KEY
#endif

// likely/unlikely
#if (__GNUC__ >= 4 || __clang__)
#    define EMHASH_LIKELY(condition) __builtin_expect(condition, 1)
#    define EMHASH_UNLIKELY(condition) __builtin_expect(condition, 0)
#else
#    define EMHASH_LIKELY(condition) condition
#    define EMHASH_UNLIKELY(condition) condition
#endif

#define NEW_KEY(key, bucket)    new(_pairs + bucket) PairT(key, bucket), _num_filled ++; SET_BIT(bucket)
#define COPY_KEY(key, bucket)   _pairs[bucket].first = key
#if EMHASH_CACHE_LINE_SIZE < 32
    #define EMHASH_CACHE_LINE_SIZE 64
#endif

#define GET_KEY(p,n)     p[n].first
#define NEXT_BUCKET(s,n) s[n].second


#define MASK_BIT         32
#define MASK_N(n)        1 << (n % MASK_BIT)
#define SET_BIT(bucket) _bitmask[bucket / MASK_BIT] &= ~(MASK_N(bucket))
#define CLS_BIT(bucket) _bitmask[bucket / MASK_BIT] |= MASK_N(bucket)
#define IS_SET(bucket)  _bitmask[bucket / MASK_BIT] & (MASK_N(bucket))

#if _WIN32 || _MSC_VER > 1400
    #include <intrin.h>
#endif

namespace emhash {

constexpr uint32_t INACTIVE = 0xFFFFFFFF;

inline static uint32_t CTZ64(const uint64_t n)
{
#if _MSC_VER > 1400 || _WIN32
    unsigned long index;
#if __x86_64__ || __amd64__ || _M_X64
    _BitScanForward64(&index, n);
#else
    _BitScanForward(&index, n);
#endif
#elif __GNUC__ || __clang__
    uint32_t index = __builtin_ctzll(n);
#elif 1
    uint64_t index;
#if __GNUC__ || __clang__
    __asm__("bsfq %1, %0\n" : "=r" (index) : "rm" (n) : "cc");
#else
    __asm
    {
        bsfq eax, n
        mov index, eax
    }
#endif
#endif

    return (uint32_t)index;
}

/// A cache-friendly hash table with open addressing, linear probing and power-of-two capacity
template <typename KeyT, typename HashT = std::hash<KeyT>, typename EqT = std::equal_to<KeyT>>
class HashSet
{

private:
    typedef  HashSet<KeyT, HashT> htype;
    typedef  std::pair<KeyT, uint32_t> PairT;

public:
    typedef size_t   size_type;
    typedef KeyT     value_type;
    typedef KeyT&    reference;
    typedef KeyT*     pointer;

    class iterator
    {
    public:
        //typedef std::forward_iterator_tag iterator_category;
        typedef size_t                    difference_type;
        typedef size_t                    distance_type;

        iterator() { }
        iterator(htype* hash_set, uint32_t bucket) : _set(hash_set), _bucket(bucket) { }

        iterator& operator++()
        {
            this->goto_next_element();
            return *this;
        }

        iterator operator++(int)
        {
            auto old_index = _bucket;
            this->goto_next_element();
            return {_set, old_index};
        }

        reference operator*() const
        {
            return _set->GET_KEY(_pairs, _bucket);
        }

        pointer operator->() const
        {
            return &(_set->GET_KEY(_pairs, _bucket));
        }

        bool operator==(const iterator& rhs) const
        {
            return this->_bucket == rhs._bucket;
        }

        bool operator!=(const iterator& rhs) const
        {
            return this->_bucket != rhs._bucket;
        }

    private:
        void goto_next_element()
        {
            auto _bitmask = _set->_bitmask;
            do {
                _bucket++;
            } while (IS_SET(_bucket));
        }

    public:
        htype* _set;
        uint32_t  _bucket;
    };

    class const_iterator
    {
    public:
        //typedef std::forward_iterator_tag iterator_category;
        typedef size_t                    difference_type;
        typedef size_t                    distance_type;

        const_iterator() { }
        const_iterator(const iterator& proto) : _set(proto._set), _bucket(proto._bucket) {  }
        const_iterator(const htype* hash_set, uint32_t bucket) : _set(hash_set), _bucket(bucket) {  }

        const_iterator& operator++()
        {
            this->goto_next_element();
            return *this;
        }

        const_iterator operator++(int)
        {
            auto old_index = _bucket;
            this->goto_next_element();
            return {_set, old_index};
        }

        const reference operator*() const
        {
            return _set->GET_KEY(_pairs, _bucket);
        }

        const pointer operator->() const
        {
            return &(_set->GET_KEY(_pairs, _bucket));
        }

        bool operator==(const const_iterator& rhs) const
        {
            return this->_bucket == rhs._bucket;
        }

        bool operator!=(const const_iterator& rhs) const
        {
            return this->_bucket != rhs._bucket;
        }

    private:
        void goto_next_element()
        {
            auto _bitmask = _set->_bitmask;
            do {
                _bucket++;
            } while (IS_SET(_bucket));
        }

    public:
        const htype* _set;
        uint32_t  _bucket;
    };

    // ------------------------------------------------------------------------

    void init(uint32_t bucket)
    {
        _num_buckets = 0;
        _mask = 0;
        _pairs = nullptr;
        _bitmask = nullptr;
        _num_filled = 0;
        max_load_factor(0.90f);
        reserve(bucket);
    }

    HashSet(uint32_t bucket = 4)
    {
        init(bucket);
    }

    HashSet(const HashSet& other)
    {
        _pairs = (PairT*)malloc((2 + other._num_buckets) * sizeof(PairT) + other._num_buckets / 8 + sizeof(uint64_t));
        clone(other);
    }

    HashSet(HashSet&& other)
    {
        init(1);
        *this = std::move(other);
    }

    HashSet(std::initializer_list<value_type> il, int n = 8)
    {
        init((uint32_t)il.size());
        for (auto begin = il.begin(); begin != il.end(); ++begin)
            insert(*begin);
    }

    HashSet& operator=(const HashSet& other)
    {
        if (this == &other)
            return *this;

        if (!std::is_pod<KeyT>::value)
            clearkv();

        if (_num_buckets != other._num_buckets) {
            free(_pairs);
            _pairs = (PairT*)malloc((2 + other._num_buckets) * sizeof(PairT) + other._num_buckets / 8 + sizeof(uint64_t));
        }

        clone(other);
        return *this;
    }

    HashSet& operator=(HashSet&& other)
    {
        this->swap(other);
        return *this;
    }

    ~HashSet()
    {
        if (!std::is_pod<KeyT>::value)
            clearkv();

        free(_pairs);
    }

    void clone(const HashSet& other)
    {
        _hasher      = other._hasher;
        _num_buckets = other._num_buckets;
        _num_filled  = other._num_filled;
        _mask        = other._mask;
        _loadlf      = other._loadlf;
        _bitmask     = (uint32_t*)(_pairs + 2 + _num_buckets);

        auto opairs = other._pairs;
        if (std::is_pod<KeyT>::value) {
            memcpy(_pairs, opairs, other._num_buckets * sizeof(PairT));
        } else {
            for (uint32_t bucket = 0; bucket < _num_buckets; bucket++) {
                auto next_bucket = NEXT_BUCKET(_pairs, bucket) = NEXT_BUCKET(opairs, bucket);
                if (next_bucket != INACTIVE)
                    new(_pairs + bucket) PairT(opairs[bucket]);
            }
        }
        memcpy(_pairs + _num_buckets, opairs + _num_buckets, 2 * sizeof(PairT) + _num_buckets / 8 + sizeof(uint64_t));
    }
    void swap(HashSet& other)
    {
        std::swap(_hasher, other._hasher);
        std::swap(_pairs, other._pairs);
        std::swap(_num_buckets, other._num_buckets);
        std::swap(_num_filled, other._num_filled);
        std::swap(_mask, other._mask);
        std::swap(_loadlf, other._loadlf);
        std::swap(_bitmask, other._bitmask);
    }

    // -------------------------------------------------------------

    iterator begin()
    {
        uint32_t bucket = 0;
        while (IS_SET(bucket)) {
            ++bucket;
        }
        return {this, bucket};
    }

    const_iterator cbegin() const
    {
        uint32_t bucket = 0;
        while (IS_SET(bucket)) {
            ++bucket;
        }
        return {this, bucket};
    }

    const_iterator begin() const
    {
        return cbegin();
    }

    iterator end()
    {
        return {this, _num_buckets};
    }

    const_iterator cend() const
    {
        return {this, _num_buckets};
    }

    const_iterator end() const
    {
        return {this, _num_buckets};
    }

    size_type size() const
    {
        return _num_filled;
    }

    bool empty() const
    {
        return _num_filled == 0;
    }

    // Returns the number of buckets.
    size_type bucket_count() const
    {
        return _num_buckets;
    }

    /// Returns average number of elements per bucket.
    float load_factor() const
    {
        return static_cast<float>(_num_filled) / (_num_buckets + 1);
    }

    HashT hash_function() const
    {
        return _hasher;
    }

    EqT key_eq() const
    {
        return _eq;
    }

    constexpr float max_load_factor() const
    {
        return (1 << 17) / (float)_loadlf;
    }

    void max_load_factor(float value)
    {
        if (value < 0.995f && value > 0.2f)
            _loadlf = (uint32_t)((1 << 17) / value);
    }

    constexpr size_type max_size() const
    {
        return (1 << 30) / sizeof(PairT);
    }

    constexpr size_type max_bucket_count() const
    {
        return (1 << 30) / sizeof(PairT);
    }

#ifdef EMHASH_STATIS
    //Returns the bucket number where the element with key k is located.
    size_type bucket(const KeyT& key) const
    {
        const auto bucket = hash_bucket(key);
        const auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        if (next_bucket == INACTIVE)
            return 0;
        else if (bucket == next_bucket)
            return bucket + 1;

        const auto& bucket_key = GET_KEY(_pairs, bucket);
        return hash_bucket(bucket_key) + 1;
    }

    //Returns the number of elements in bucket n.
    size_type bucket_size(const size_type bucket) const
    {
        auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        if (next_bucket == INACTIVE)
            return 0;

        const auto& bucket_key = GET_KEY(_pairs, bucket);
        next_bucket = hash_bucket(bucket_key);
        uint32_t ibucket_size = 1;

        //iterator each item in current main bucket
        while (true) {
            const auto nbucket = NEXT_BUCKET(_pairs, next_bucket);
            if (nbucket == next_bucket) {
                break;
            }
            ibucket_size ++;
            next_bucket = nbucket;
        }
        return ibucket_size;
    }

    size_type get_main_bucket(const uint32_t bucket) const
    {
        auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        if (next_bucket == INACTIVE)
            return INACTIVE;

        const auto& bucket_key = GET_KEY(_pairs, bucket);
        const auto main_bucket = hash_bucket(bucket_key);
        return main_bucket;
    }

    int get_cache_info(uint32_t bucket, uint32_t next_bucket) const
    {
        auto pbucket = reinterpret_cast<size_t>(&_pairs[bucket]);
        auto pnext   = reinterpret_cast<size_t>(&_pairs[next_bucket]);
        if (pbucket / 64 == pnext / 64)
            return 0;
        auto diff = pbucket > pnext ? (pbucket - pnext) : pnext - pbucket;
        if (diff < 127 * 64)
            return diff / 64 + 1;
        return 127;
    }

    int get_bucket_info(const uint32_t bucket, uint32_t steps[], const uint32_t slots) const
    {
        auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        if (next_bucket == INACTIVE)
            return -1;

        const auto& bucket_key = GET_KEY(_pairs, bucket);
        const auto main_bucket = hash_bucket(bucket_key);
        if (main_bucket != bucket)
            return 0;
        else if (next_bucket == bucket)
            return 1;

        steps[get_cache_info(bucket, next_bucket) % slots] ++;
        uint32_t ibucket_size = 2;
        //find a new empty and linked it to tail
        while (true) {
            const auto nbucket = NEXT_BUCKET(_pairs, next_bucket);
            if (nbucket == next_bucket)
                break;

            steps[get_cache_info(nbucket, next_bucket) % slots] ++;
            ibucket_size ++;
            next_bucket = nbucket;
        }
        return ibucket_size;
    }

    void dump_statis() const
    {
        uint32_t buckets[129] = {0};
        uint32_t steps[129]   = {0};
        for (uint32_t bucket = 0; bucket < _num_buckets; ++bucket) {
            auto bsize = get_bucket_info(bucket, steps, 128);
            if (bsize > 0)
                buckets[bsize] ++;
        }

        uint32_t sumb = 0, collision = 0, sumc = 0, finds = 0, sumn = 0;
        puts("============== buckets size ration =========");
        for (uint32_t i = 0; i < sizeof(buckets) / sizeof(buckets[0]); i++) {
            const auto bucketsi = buckets[i];
            if (bucketsi == 0)
                continue;
            sumb += bucketsi;
            sumn += bucketsi * i;
            collision += bucketsi * (i - 1);
            finds += bucketsi * i * (i + 1) / 2;
            printf("  %2u  %8u  %.2lf  %.2lf\n", i, bucketsi, bucketsi * 100.0 * i / _num_filled, sumn * 100.0 / _num_filled);
        }

        puts("========== collision miss ration ===========");
        for (uint32_t i = 0; i < sizeof(steps) / sizeof(steps[0]); i++) {
            sumc += steps[i];
            if (steps[i] <= 2)
                continue;
            printf("  %2u  %8u  %.2lf  %.2lf\n", i, steps[i], steps[i] * 100.0 / collision, sumc * 100.0 / collision);
        }

        if (sumb == 0)  return;
        printf("    _num_filled/bucket_size/packed collision/cache_miss/hit_find = %u/%.2lf/%zd/ %.2lf%%/%.2lf%%/%.2lf\n",
                _num_filled, _num_filled * 1.0 / sumb, sizeof(PairT), (collision * 100.0 / _num_filled), (collision - steps[0]) * 100.0 / _num_filled, finds * 1.0 / _num_filled);
        assert(sumn == _num_filled);
        assert(sumc == collision);
    }
#endif

    // ------------------------------------------------------------

    iterator find(const KeyT& key) noexcept
    {
        return {this, find_filled_bucket(key)};
    }

    const_iterator find(const KeyT& key) const noexcept
    {
        return {this, find_filled_bucket(key)};
    }

    bool contains(const KeyT& key) const noexcept
    {
        return find_filled_bucket(key) != _num_buckets;
    }

    size_type count(const KeyT& key) const noexcept
    {
        return find_filled_bucket(key) == _num_buckets ? 0 : 1;
    }

    /// Returns a pair consisting of an iterator to the inserted element
    /// (or to the element that prevented the insertion)
    /// and a bool denoting whether the insertion took place.
    std::pair<iterator, bool> insert(const KeyT& key)
    {
        check_expand_need();
        const auto bucket = find_or_allocate(key);
        if (NEXT_BUCKET(_pairs, bucket) == INACTIVE) {
            NEW_KEY(key, bucket);
            return { {this, bucket}, true };
        } else {
            return { {this, bucket}, false };
        }
    }

    std::pair<iterator, bool> insert(KeyT&& key)
    {
        check_expand_need();
        const auto bucket = find_or_allocate(key);
        if (NEXT_BUCKET(_pairs, bucket) == INACTIVE) {
            NEW_KEY(std::move(key), bucket);
            return { {this, bucket}, true };
        } else {
            return { {this, bucket}, false };
        }
    }

#if 0
    template <typename Iter>
    inline void insert(Iter begin, Iter end)
    {
        reserve(end - begin + _num_filled);
        for (; begin != end; ++begin) {
            insert(*begin);
        }
    }
#endif

    void insert(std::initializer_list<value_type> ilist)
    {
        reserve((uint32_t)ilist.size() + _num_filled);
        for (auto begin = ilist.begin(); begin != ilist.end(); ++begin) {
            insert(*begin);
        }
    }

    template <typename Iter>
    inline void insert(Iter begin, Iter end)
    {
        Iter citbeg = begin;
        Iter citend = begin;
        reserve(end - begin + _num_filled);
        for (; begin != end; ++begin) {
            if (try_insert_mainbucket(*begin) == INACTIVE) {
                std::swap(*begin, *citend++);
            }
        }

        for (; citbeg != citend; ++citbeg) {
            auto& key = *citbeg;
            const auto bucket = find_or_allocate(key);
            if (NEXT_BUCKET(_pairs, bucket) == INACTIVE) {
                NEW_KEY(key, bucket);
            }
        }
    }

    template <typename Iter>
    inline void insert_unique(Iter begin, Iter end)
    {
        reserve(end - begin + _num_filled);
        for (; begin != end; ++begin) {
            insert_unique(*begin);
        }
    }

    /// Same as above, but contains(key) MUST be false
    uint32_t insert_unique(const KeyT& key)
    {
        check_expand_need();
        auto bucket = find_unique_bucket(key);
        NEW_KEY(key, bucket);
        return bucket;
    }

    uint32_t insert_unique(KeyT&& key)
    {
        check_expand_need();
        auto bucket = find_unique_bucket(key);
        NEW_KEY(std::move(key), bucket);
        return bucket;
    }

    //not
    template <class... Args>
    inline std::pair<iterator, bool> emplace(Args&&... args)
    {
        return insert(std::forward<Args>(args)...);
    }

    //no any optimize for position
    template <class... Args>
    iterator emplace_hint(const_iterator position, Args&&... args)
    {
        return insert(std::forward<Args>(args)...).first;
    }
    std::pair<iterator, bool> try_emplace(const value_type& k)
    {
        return insert(k).first;
    }
    template <class... Args>
    inline std::pair<iterator, bool> emplace_unique(Args&&... args)
    {
        return insert_unique(std::forward<Args>(args)...);
    }

    uint32_t try_insert_mainbucket(const KeyT& key)
    {
        auto bucket = hash_bucket(key);
        auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        if (next_bucket == INACTIVE) {
            NEW_KEY(key, bucket);
            return bucket;
        } else if(_eq(key, GET_KEY(_pairs, bucket)))
            return bucket;

        return INACTIVE;
    }

    void insert_or_assign(const KeyT& key)
    {
        check_expand_need();
        const auto bucket = find_or_allocate(key);
        // Check if inserting a new value rather than overwriting an old entry
        if (NEXT_BUCKET(_pairs, bucket) == INACTIVE) {
            NEW_KEY(key, bucket);
        } else {
            NEXT_BUCKET(_pairs, bucket) = bucket;
        }
    }

    void clear_bucket(uint32_t bucket)
    {
        _pairs[bucket].~PairT();
        NEXT_BUCKET(_pairs, bucket) = INACTIVE;
        _num_filled --;
        CLS_BIT(bucket);
    }

    // -------------------------------------------------------
    /// Erase an element from the hash table.
    /// return 0 if element was not found
    size_type erase(const KeyT& key) noexcept
    {
        const auto bucket = erase_key(key);
        if (bucket == INACTIVE)
            return 0;

        clear_bucket(bucket);
        return 1;
    }

    //iterator erase(const_iterator begin_it, const_iterator end_it)
    iterator erase(const_iterator cit) noexcept
    {
        iterator it(this, cit._bucket);
        return erase(it);
    }

    /// Erase an element typedef an iterator.
    /// Returns an iterator to the next element (or end()).
    iterator erase(iterator it) noexcept
    {
        const auto bucket = erase_bucket(it._bucket);
        clear_bucket(bucket);
        //erase from main bucket, return main bucket as next
        return (bucket == it._bucket) ? ++it : it;
    }

    void _erase(const_iterator it)
    {
        const auto bucket = erase_bucket(it._bucket);
        clear_bucket(bucket);
    }

    void clearkv()
    {
        for (uint32_t bucket = 0; _num_filled > 0; ++bucket) {
            if (NEXT_BUCKET(_pairs, bucket) != INACTIVE)
                clear_bucket(bucket);
        }
    }

    /// Remove all elements, keeping full capacity.
    void clear()
    {
        if (!std::is_pod<KeyT>::value || sizeof(PairT) > EMHASH_CACHE_LINE_SIZE || _num_filled < _num_buckets / 4)
            clearkv();
        else {
            memset(_pairs, INACTIVE, sizeof(_pairs[0]) * _num_buckets);
            memset(_bitmask, INACTIVE, _num_buckets / 8);
            _bitmask[_num_buckets / MASK_BIT] &= (1 << _num_buckets % MASK_BIT) - 1;
        }
        _num_filled = 0;
    }

    void shrink_to_fit()
    {
        rehash(_num_filled);
    }

    /// Make room for this many elements
    bool reserve(uint32_t num_elems) noexcept
    {
        const auto required_buckets = (uint32_t)(((uint64_t)num_elems) * _loadlf >> 17);
        if (EMHASH_LIKELY(required_buckets < _mask))
            return false;

        rehash(required_buckets + 2);
        return true;
    }

private:
    void rehash(uint32_t required_buckets) noexcept
    {
        if (required_buckets < _num_filled)
            return ;

        uint32_t num_buckets = _num_filled > 65536 ? (1 << 16) : 4;
        while (num_buckets < required_buckets) { num_buckets *= 2; }

        const auto num_byte = num_buckets / 8;
        auto new_pairs = (PairT*)malloc((2 + num_buckets) * sizeof(PairT) + num_byte + sizeof(uint64_t));
        auto old_num_filled  = _num_filled;
        auto old_pairs = _pairs;

        _num_filled  = 0;
        _num_buckets = num_buckets;
        _mask        = num_buckets - 1;

        _bitmask     = (uint32_t*)(new_pairs + 2 + num_buckets);

        if (sizeof(PairT) <= EMHASH_CACHE_LINE_SIZE / 2)
            memset(new_pairs, INACTIVE, sizeof(_pairs[0]) * num_buckets);
        else
            for (uint32_t bucket = 0; bucket < num_buckets; bucket++)
                NEXT_BUCKET(new_pairs, bucket) = INACTIVE;

        /***************** ----------------------**/
        memset(_bitmask, INACTIVE, num_byte);
        memset(((char*)_bitmask) + num_byte, 0, sizeof(uint64_t));
        //set high bit to zero
        _bitmask[num_buckets / MASK_BIT] &= (1 << num_buckets % MASK_BIT) - 1;
        /***************** ----------------------**/
        memset(new_pairs + _num_buckets, 0, sizeof(PairT) * 2);
        _pairs       = new_pairs;
        for (uint32_t src_bucket = 0; _num_filled < old_num_filled; src_bucket++) {
            if (NEXT_BUCKET(old_pairs, src_bucket) == INACTIVE)
                continue;

            auto& key = GET_KEY(old_pairs, src_bucket);
            const auto bucket = find_unique_bucket(key);
            NEW_KEY(std::move(key), bucket);
            old_pairs[src_bucket].~PairT();
        }

#if EMHASH_REHASH_LOG
        if (_num_filled > 100000) {
            char buff[255] = {0};
            sprintf(buff, "    _num_filled/K/pack= %u/%s/%zd/",
                    _num_filled, typeid(KeyT).name(), sizeof(_pairs[0]));
#if EMHASH_TAF_LOG
            static uint32_t ihashs = 0;
            FDLOG() << "|hash_nums = " << ihashs ++ << "|" <<__FUNCTION__ << "|" << buff << endl;
#else
            puts(buff);
#endif
        }
#endif

        free(old_pairs);
        assert(old_num_filled == _num_filled);
    }

private:
    // Can we fit another element?
    inline bool check_expand_need()
    {
        return reserve(_num_filled);
    }

    uint32_t erase_key(const KeyT& key)
    {
        const auto bucket = hash_bucket(key);
        auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        if (next_bucket == INACTIVE)
            return INACTIVE;

        const auto eqkey = _eq(key, GET_KEY(_pairs, bucket));
        if (next_bucket == bucket) {
            return eqkey ? bucket : INACTIVE;
         } else if (eqkey) {
            const auto nbucket = NEXT_BUCKET(_pairs, next_bucket);
            if (std::is_pod<KeyT>::value)
                GET_KEY(_pairs, bucket) = GET_KEY(_pairs, next_bucket);
            else
                std::swap(GET_KEY(_pairs, bucket), GET_KEY(_pairs, next_bucket));
            NEXT_BUCKET(_pairs, bucket) = (nbucket == next_bucket) ? bucket : nbucket;
            return next_bucket;
        } else if (EMHASH_UNLIKELY(bucket != hash_bucket(GET_KEY(_pairs, bucket))))
            return INACTIVE;

        auto prev_bucket = bucket;
        while (true) {
            const auto nbucket = NEXT_BUCKET(_pairs, next_bucket);
            if (_eq(key, GET_KEY(_pairs, next_bucket))) {
                NEXT_BUCKET(_pairs, prev_bucket) = (nbucket == next_bucket) ? prev_bucket : nbucket;
                return next_bucket;
            }

            if (nbucket == next_bucket)
                break;
            prev_bucket = next_bucket;
            next_bucket = nbucket;
        }

        return INACTIVE;
    }

    uint32_t erase_bucket(const uint32_t bucket)
    {
        const auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        const auto main_bucket = hash_bucket(GET_KEY(_pairs, bucket));
        if (bucket == main_bucket) {
            if (bucket != next_bucket) {
                const auto nbucket = NEXT_BUCKET(_pairs, next_bucket);
                if (!std::is_pod<KeyT>::value)
                    std::swap(GET_KEY(_pairs, bucket), GET_KEY(_pairs, next_bucket));
                else
                    GET_KEY(_pairs, bucket) = GET_KEY(_pairs, next_bucket);
                NEXT_BUCKET(_pairs, bucket) = (nbucket == next_bucket) ? bucket : nbucket;
            }
            return next_bucket;
        }

        const auto prev_bucket = find_prev_bucket(main_bucket, bucket);
        NEXT_BUCKET(_pairs, prev_bucket) = (bucket == next_bucket) ? prev_bucket : next_bucket;
        return bucket;
    }

    // Find the bucket with this key, or return bucket size
    uint32_t find_filled_bucket(const KeyT& key) const
    {
        const auto bucket = hash_bucket(key);
        auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        const auto& bucket_key = GET_KEY(_pairs, bucket);
        if (next_bucket == INACTIVE)
            return _num_buckets;
        else if (_eq(key, bucket_key))
            return bucket;
        else if (next_bucket == bucket)
            return _num_buckets;
//        else if (bucket != (hash_bucket(bucket_key)))
//            return _num_buckets;

        while (true) {
            if (_eq(key, GET_KEY(_pairs, next_bucket)))
                return next_bucket;

            const auto nbucket = NEXT_BUCKET(_pairs, next_bucket);
            if (nbucket == next_bucket)
                break;
            next_bucket = nbucket;
        }

        return _num_buckets;
    }

    //main --> prev --> bucekt -->next --> new
    uint32_t kickout_bucket(const uint32_t main_bucket, const uint32_t bucket)
    {
        const auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        const auto new_bucket  = find_empty_bucket(next_bucket);
        const auto prev_bucket = find_prev_bucket(main_bucket, bucket);
        NEXT_BUCKET(_pairs, prev_bucket) = new_bucket;
        new(_pairs + new_bucket) PairT(std::move(_pairs[bucket]));
        SET_BIT(new_bucket);
        if (next_bucket == bucket)
            NEXT_BUCKET(_pairs, new_bucket) = new_bucket;
         _pairs[bucket].~PairT(); NEXT_BUCKET(_pairs, bucket) = INACTIVE;
        return bucket;
    }

/*
** inserts a new key into a hash table; first, check whether key's main
** bucket/position is free. If not, check whether colliding node/bucket is in its main
** position or not: if it is not, move colliding bucket to an empty place and
** put new key in its main position; otherwise (colliding bucket is in its main
** position), new key goes to an empty position.
*/
    uint32_t find_or_allocate(const KeyT& key)
    {
        const auto bucket = hash_bucket(key);
        const auto& bucket_key = GET_KEY(_pairs, bucket);
        auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        if (next_bucket == INACTIVE || _eq(key, bucket_key))
            return bucket;

        //check current bucket_key is in main bucket or not
        const auto main_bucket = hash_bucket(bucket_key);
        if (main_bucket != bucket)
            return kickout_bucket(main_bucket, bucket);
        else if (next_bucket == bucket)
            return NEXT_BUCKET(_pairs, next_bucket) = find_empty_bucket(next_bucket);

        //find next linked bucket and check key
        while (true) {
            if (_eq(key, GET_KEY(_pairs, next_bucket))) {
#if EMHASH_LRU_SET
                std::swap(GET_KEY(_pairs, bucket), GET_KEY(_pairs, next_bucket));
                return bucket;
#else
                return next_bucket;
#endif
            }

            const auto nbucket = NEXT_BUCKET(_pairs, next_bucket);
            if (nbucket == next_bucket)
                break;
            next_bucket = nbucket;
        }

        //find a new empty and link it to tail
        const auto new_bucket = find_empty_bucket(next_bucket);
        return NEXT_BUCKET(_pairs, next_bucket) = new_bucket;
    }

    // key is not in this map. Find a place to put it.
    uint32_t find_empty_bucket(const uint32_t bucket_from)
    {
        const auto bucket1 = bucket_from + 1;
        if (NEXT_BUCKET(_pairs, bucket1) == INACTIVE)
            return bucket1;

        const auto bucket2 = bucket_from + 2;
        if (NEXT_BUCKET(_pairs, bucket2) == INACTIVE)
            return bucket2;

        const auto boset = bucket_from % 8;
        const uint64_t bmask = *(uint64_t*)((uint8_t*)_bitmask + bucket_from / 8) >> boset;
        if (bmask != 0) {
            return bucket_from + CTZ64(bmask);
        }

        const auto qmask = (64 + _num_buckets - 1) / 64 - 1;
        for (uint32_t last = 2, next2 = (bucket_from + _num_filled) & qmask; ; last ++, next2 = (next2 + last) & qmask) {
            const auto bmask2 = *((uint64_t*)_bitmask + next2);
            if (bmask2 != 0) {
                return next2 * 64 + CTZ64(bmask2);
            }
        }

        return 0;
    }

    uint32_t find_last_bucket(uint32_t main_bucket) const
    {
        auto next_bucket = NEXT_BUCKET(_pairs, main_bucket);
        if (next_bucket == main_bucket)
            return main_bucket;

        while (true) {
            const auto nbucket = NEXT_BUCKET(_pairs, next_bucket);
            if (nbucket == next_bucket)
                return next_bucket;
            next_bucket = nbucket;
        }
    }

    uint32_t find_prev_bucket(uint32_t main_bucket, const uint32_t bucket) const
    {
        auto next_bucket = NEXT_BUCKET(_pairs, main_bucket);
        if (next_bucket == bucket)
            return main_bucket;

        while (true) {
            const auto nbucket = NEXT_BUCKET(_pairs, next_bucket);
            if (nbucket == bucket)
                return next_bucket;
            next_bucket = nbucket;
        }
    }

    uint32_t find_unique_bucket(const KeyT& key)
    {
        const auto bucket = hash_bucket(key);
        auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        if (next_bucket == INACTIVE)
            return bucket;

        //check current bucket_key is in main bucket or not
        const auto main_bucket = hash_bucket(GET_KEY(_pairs, bucket));
        if (main_bucket != bucket)
            return kickout_bucket(main_bucket, bucket);
        else if (next_bucket != bucket)
            next_bucket = find_last_bucket(next_bucket);

        //find a new empty and link it to tail
        return NEXT_BUCKET(_pairs, next_bucket) = find_empty_bucket(next_bucket);
    }

    static inline uint32_t hash32(uint32_t key)
    {
#if 0
        key = ((key >> 16) ^ key) * 0x45d9f3b;
        key = ((key >> 16) ^ key) * 0x45d9f3b; //0x119de1f3
//        key = ((key >> 13) ^ key) * 0xc2b2ae35;
        key = (key >> 16) ^ key;
        return key;
#elif 1
        uint64_t const r = key * UINT64_C(2654435769);
        const uint32_t h = static_cast<uint32_t>(r >> 32);
        const uint32_t l = static_cast<uint32_t>(r);
        return h + l;
#elif 1
        key += ~(key << 15);
        key ^= (key >> 10);
        key += (key << 3);
        key ^= (key >> 6);
        key += ~(key << 11);
        key ^= (key >> 16);
        return key;
#endif
    }

    static inline uint64_t hash64(uint64_t key)
    {
#if __SIZEOF_INT128__ && _MPCLMUL
        //uint64_t const inline clmul_mod(const uint64_t& i,const uint64_t& j)
        __m128i I{}; I[0] ^= key;
        __m128i J{}; J[0] ^= UINT64_C(0xde5fb9d2630458e9);
        __m128i M{}; M[0] ^= 0xb000000000000000ull;
        __m128i X = _mm_clmulepi64_si128(I,J,0);
        __m128i A = _mm_clmulepi64_si128(X,M,0);
        __m128i B = _mm_clmulepi64_si128(A,M,0);
        return A[0]^A[1]^B[1]^X[0]^X[1];
#elif __SIZEOF_INT128__
        constexpr uint64_t k = UINT64_C(11400714819323198485);
        __uint128_t r = key; r *= k;
        return (uint32_t)(r >> 64) + (uint32_t)r;
#elif 1
        uint64_t const r = key * UINT64_C(0xca4bcaa75ec3f625);
        const uint32_t h = static_cast<uint32_t>(r >> 32);
        const uint32_t l = static_cast<uint32_t>(r);
        return h + l;
#elif 1
        //MurmurHash3Mixer
        uint64_t h = key;
        h ^= h >> 33;
        h *= 0xff51afd7ed558ccd;
        h ^= h >> 33;
        h *= 0xc4ceb9fe1a85ec53;
        h ^= h >> 33;
        return static_cast<size_t>(h);
#elif 1
        uint64_t x = key;
        x = (x ^ (x >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
        x = (x ^ (x >> 27)) * UINT64_C(0x94d049bb133111eb);
        x = x ^ (x >> 31);
        return x;
#endif
    }

    template<typename UType, typename std::enable_if<std::is_integral<UType>::value, uint32_t>::type = 0>
    inline uint32_t hash_bucket(const UType key) const
    {
#ifdef EMHASH_FIBONACCI_HASH
        if (sizeof(UType) <= sizeof(uint32_t))
            return hash32(key) & _mask;
        else
            return uint32_t(hash64(key) & _mask);
#elif EMHASH_SAFE_HASH
        if (_hash_inter > 0) {
            if (sizeof(UType) <= sizeof(uint32_t))
                return hash32(key) & _mask;
            else
                return uint32_t(hash64(key) & _mask);
        }
        return _hasher(key) & _mask;
#elif EMHASH_IDENTITY_HASH
        return (key + (key >> 20)) & _mask;
#else
        return _hasher(key) & _mask;
#endif
    }

    template<typename UType, typename std::enable_if<!std::is_integral<UType>::value, uint32_t>::type = 0>
    inline uint32_t hash_bucket(const UType& key) const
    {
#ifdef EMHASH_FIBONACCI_HASH
        return (_hasher(key) * 11400714819323198485ull) & _mask;
#elif EMHASH_STD_STRING
        uint32_t hash = 0;
        if (key.size() < 32) {
            for (const auto c : key) hash = c + hash * 131;
        } else {
            for (int i = 0, j = 1; i < key.size(); i += j++)
                hash = key[i] + hash * 131;
        }
        return hash & _mask;
#else
        return _hasher(key) & _mask;
#endif
    }

private:

    //the first cache line packed
    HashT     _hasher;
    EqT       _eq;
    uint32_t  _loadlf;
    uint32_t  _num_buckets;
    uint32_t  _mask;

    uint32_t  _num_filled;
    PairT*    _pairs;
    uint32_t* _bitmask;
};
} // namespace emhash
#if __cplusplus >= 201103L
//template <class Key, typename Hash = std::hash<Key>> using ktprime_hashset = emhash6::HashSet<Key, Hash>;
#endif
