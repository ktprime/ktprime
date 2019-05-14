// By Huang Yuanbing 2019
// bailuzhou@163.com
// https://github.com/ktprime/ktprime/blob/master/hash_table6.hpp

// LICENSE:
//   This software is dual-licensed to the public domain and under the following
//   license: you are granted a perpetual, irrevocable license to copy, modify,
//   publish, and distribute this file as you see fit.

/*****
    This a very fast and efficent c++ hash map implemention and foucus optimization with small k/v.
some Features as fllow:
1. combine Linear Probing and Quadratic Probing
2. open Addressing with Linked collsion slot
3. specially optimization with cpu cache line size
4. very fast/good rehash algorithm
5. high load factor(0.9) and can save memory than unordered_map.
6. can set lru get/set by set marco
7.
*/

#pragma once

#include <cstring>
#include <cassert>
#if EMILIB_TAF_LOG
    #include "servant/AutoLog.h"
    #include "servant/RollLogHelper.h"
#endif

#ifdef GET_KEY
    #undef  GET_KEY
    #undef  GET_VAL
    #undef  NEXT_BUCKET
    #undef  GET_PVAL
#endif

#ifndef EMILIB_ORDER_INDEX
    #define EMILIB_ORDER_INDEX 2
#endif
#if EMILIB_CACHE_LINE_SIZE < 32
    #define EMILIB_CACHE_LINE_SIZE 64
#endif

#if EMILIB_ORDER_INDEX == 0
    #define GET_KEY(p,n)     p[n].second.first
    #define GET_VAL(p,n)     p[n].second.second
    #define NEXT_BUCKET(s,n) s[n].first
    #define GET_PVAL(s,n)    s[n].second
    #define NEW_KVALUE(key, value, bucket) new(_pairs + bucket) PairT(bucket, std::pair<KeyT, ValueT>(key, value))
#elif EMILIB_ORDER_INDEX == 1
    #define GET_KEY(p,n)     p[n].first.first
    #define GET_VAL(p,n)     p[n].first.second
    #define NEXT_BUCKET(s,n) s[n].second
    #define GET_PVAL(s,n)    s[n].first
    #define NEW_KVALUE(key, value, bucket) new(_pairs + bucket) PairT(std::pair<KeyT, ValueT>(key, value), bucket)
#else
    #define GET_KEY(p,n)     p[n].first
    #define GET_VAL(p,n)     p[n].second
    #define NEXT_BUCKET(s,n) s[n].nextbucket
    #define GET_PVAL(s,n)    s[n]
    #define NEW_KVALUE(key, value, bucket) new(_pairs + bucket) PairT(key, value, bucket)
#endif

namespace emilib1 {
/// like std::equal_to but no need to #include <functional>

template <typename First, typename Second>
struct myPair {
    myPair(const First& key, const Second& value, uint32_t bucket)
        :first(key), second(value)
    {
        nextbucket = bucket;
    }

    myPair(const std::pair<First,Second>& pair)
        :first(pair.first), second(pair.second)
    {
        nextbucket = -1;
    }

    myPair(std::pair<First, Second>&& pair)
        :first(std::move(pair.first)), second(std::move(pair.second))
    {
        nextbucket = -1;
    }

    myPair(const myPair& pairT)
        :first(pairT.first), second(pairT.second)
    {
        nextbucket = pairT.nextbucket;
    }

    myPair(myPair&& pairT)
        :first(std::move(pairT.first)), second(std::move(pairT.second))
    {
        nextbucket = pairT.nextbucket;
    }

    myPair& operator = (myPair&& pairT)
    {
        first = std::move(pairT.first);
        second = std::move(pairT.second);
        nextbucket = pairT.nextbucket;
        return *this;
    }

    void swap(myPair<First, Second>& o)
    {
        std::swap(first, o.first);
        std::swap(second, o.second);
    }

    First  first;
    uint32_t nextbucket;
    Second second;
};

/// A cache-friendly hash table with open addressing, linear probing and power-of-two capacity
template <typename KeyT, typename ValueT, typename HashT = std::hash<KeyT>>
class HashMap
{
    static constexpr unsigned INACTIVE = 0xFFFFFFFF;

private:
    typedef  HashMap<KeyT, ValueT, HashT> MyType;

#if EMILIB_ORDER_INDEX == 0
    typedef std::pair<uint32_t, std::pair<KeyT, ValueT>> PairT;
#elif EMILIB_ORDER_INDEX == 1
    typedef std::pair<std::pair<KeyT, ValueT>, uint32_t> PairT;
#else
    typedef myPair<KeyT, ValueT>                    PairT;
#endif

public:
    typedef  size_t       size_type;
    typedef  PairT        value_type;
    typedef  PairT&       reference;
    typedef  const PairT& const_reference;

    class iterator
    {
    public:
        typedef std::forward_iterator_tag iterator_category;
        typedef size_t                    difference_type;
        typedef size_t                    distance_type;

#if EMILIB_ORDER_INDEX > 1
        typedef PairT                     value_type;
#else
        typedef std::pair<KeyT, ValueT>   value_type;
#endif

        typedef value_type*               pointer;
        typedef value_type&               reference;

        iterator() { }

        iterator(MyType* hash_map, uint32_t bucket) : _map(hash_map), _bucket(bucket)
        {
        }

        iterator& operator++()
        {
            this->goto_next_element();
            return *this;
        }

        iterator operator++(int)
        {
            auto old_index = _bucket;
            this->goto_next_element();
            return iterator(_map, old_index);
        }

        reference operator*() const
        {
            return _map->GET_PVAL(_pairs, _bucket);
        }

        pointer operator->() const
        {
            return &(_map->GET_PVAL(_pairs, _bucket));
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
            do {
                _bucket++;
            } while (_bucket < _map->_num_buckets && _map->NEXT_BUCKET(_pairs, _bucket) == INACTIVE);
        }

    public:
        MyType* _map;
        uint32_t  _bucket;
    };

    class const_iterator
    {
    public:
        typedef std::forward_iterator_tag iterator_category;
        typedef size_t                    difference_type;
        typedef size_t                    distance_type;
#if EMILIB_ORDER_INDEX > 1
        typedef PairT                     value_type;
#else
        typedef std::pair<KeyT, ValueT>   value_type;
#endif

        typedef value_type*               pointer;
        typedef value_type&               reference;

        const_iterator() { }

        const_iterator(iterator proto) : _map(proto._map), _bucket(proto._bucket)
        {
        }

        const_iterator(const MyType* hash_map, uint32_t bucket) : _map(hash_map), _bucket(bucket)
        {
        }

        const_iterator& operator++()
        {
            this->goto_next_element();
            return *this;
        }

        const_iterator operator++(int)
        {
            auto old_index = _bucket;
            this->goto_next_element();
            return const_iterator(_map, old_index);
        }

        reference operator*() const
        {
            return _map->GET_PVAL(_pairs, _bucket);
        }

        pointer operator->() const
        {
            return &(_map->GET_PVAL(_pairs, _bucket));
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
            do {
                _bucket++;
            } while (_bucket < _map->_num_buckets && _map->NEXT_BUCKET(_pairs, _bucket) == INACTIVE);
        }

    public:
        const MyType* _map;
        uint32_t  _bucket;
    };

    // ------------------------------------------------------------------------

    void init()
    {
        _num_buckets = 0;
        _num_filled = 0;
        _mask = 0;
        _pairs = nullptr;
        _max_load_factor = 0.90f;
        _load_threshold = (uint32_t)(4 * _max_load_factor);
    }

    HashMap(uint32_t bucket = 8)
    {
        init();
        reserve(bucket);
    }

    HashMap(const HashMap& other)
    {
        _hasher      = other._hasher;
        _num_buckets = other._num_buckets;
        _num_filled  = other._num_filled;
        _mask        = other._mask;
        _max_load_factor = other._max_load_factor;
        _load_threshold  = other._load_threshold;

        _pairs = (PairT*)malloc(_num_buckets * sizeof(PairT));

        auto old_pairs = other._pairs;
        if (sizeof(PairT) <= 24) {
            memcpy(_pairs, old_pairs, _num_buckets * sizeof(PairT));
        } else {
            for (uint32_t bucket = 0; bucket < _num_buckets; bucket++) {
                auto state = NEXT_BUCKET(_pairs, bucket) = NEXT_BUCKET(old_pairs, bucket);
                if (state != INACTIVE)
                    new(_pairs + bucket) PairT(old_pairs[bucket]);
            }
        }
    }

    HashMap(HashMap&& other)
    {
        init();
        reserve(1);
        *this = std::move(other);
    }

    HashMap(std::initializer_list<std::pair<KeyT, ValueT>> il)
    {
        init();
        reserve(il.size());
        for (auto begin = il.begin(); begin != il.end(); ++begin)
            insert(*begin);
    }

    HashMap& operator=(const HashMap& other)
    {
        if (this == &other)
            return *this;
        HashMap tmp(other);
        this->swap(tmp);
        return *this;
    }

    HashMap& operator=(HashMap&& other)
    {
        this->swap(other);
        return *this;
    }

    ~HashMap()
    {
        clear();
        if (_pairs)
            free(_pairs);
    }

    void swap(HashMap& other)
    {
        std::swap(_hasher, other._hasher);
        std::swap(_pairs, other._pairs);
        std::swap(_num_buckets, other._num_buckets);
        std::swap(_num_filled, other._num_filled);
        std::swap(_mask, other._mask);
        std::swap(_max_load_factor, other._max_load_factor);
        std::swap(_load_threshold, other._load_threshold);
    }

    // -------------------------------------------------------------

    iterator begin()
    {
        uint32_t bucket = 0;
        while (bucket < _num_buckets && NEXT_BUCKET(_pairs, bucket) == INACTIVE) {
            ++bucket;
        }
        return iterator(this, bucket);
    }

    const_iterator cbegin() const
    {
        uint32_t bucket = 0;
        while (bucket < _num_buckets && NEXT_BUCKET(_pairs, bucket) == INACTIVE) {
            ++bucket;
        }
        return const_iterator(this, bucket);
    }

    const_iterator begin() const
    {
        return cbegin();
    }

    iterator end()
    {
        return iterator(this, _num_buckets);
    }

    const_iterator cend() const
    {
        return const_iterator(this, _num_buckets);
    }

    const_iterator end() const
    {
        return cend();
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
        return static_cast<float>(_num_filled) / static_cast<float>(_num_buckets);
    }

    HashT hash_function() const
    {
        return _hasher;
    }

    constexpr float max_load_factor() const
    {
        return _max_load_factor;
    }

    void max_load_factor(float value)
    {
        if (value < 0.95 && value > 0.2)
            _max_load_factor = value;
    }

    //... depend on system
    constexpr size_type max_size() const
    {
        return (1 << 30) / sizeof(PairT);
    }

    constexpr size_type max_bucket_count() const
    {
        return (1 << 30) / sizeof(PairT);
    }

    //Returns the bucket number where the element with key k is located.
    //
    size_type bucket(const KeyT& key) const
    {
        const auto bucket = hash_key(key);
        const auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        if (next_bucket == INACTIVE)
            return 0;
        if (bucket == next_bucket)
            return bucket + 1;

        const auto& bucket_key = GET_KEY(_pairs, bucket);
        return hash_key(bucket_key) + 1;
    }

    //Returns the number of elements in bucket n.
    //bucket < _num_buckets, count items in bucket,
    //if bucket is not the main bucket, search from main bucket
    size_type bucket_size(const size_type bucket) const
    {
        auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        if (next_bucket == INACTIVE)
            return 0;

        const auto& bucket_key = GET_KEY(_pairs, bucket);
        next_bucket = hash_key(bucket_key);
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

#ifdef EMILIB_STATIS
    //bucket < _num_buckets
    size_type get_main_bucket(const uint32_t bucket) const
    {
        auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        if (next_bucket == INACTIVE)
            return INACTIVE;

        const auto& bucket_key = GET_KEY(_pairs, bucket);
        const auto main_bucket = hash_key(bucket_key);
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

    int get_bucket_info(const uint32_t bucket, int steps[], const int slots) const
    {
        auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        if (next_bucket == INACTIVE)
            return INACTIVE;

        const auto& bucket_key = GET_KEY(_pairs, bucket);
        const auto main_bucket = hash_key(bucket_key);
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
        puts("===============  buckets ration ========= ");
        for (uint32_t i = 0; i < sizeof(buckets) / sizeof(buckets[0]); i++) {
            const auto bucketsi = buckets[i];
            if (bucketsi == 0)
                continue;
            sumb += bucketsi;
            sumn += bucketsi * i;
            collision += bucketsi * (i - 1);
            finds += bucketsi * i * (i + 1) / 2;
            printf("  %2d  %8d  %.2lf  %.2lf\n", i, bucketsi, bucketsi * 100.0 * i / _num_filled, sumn * 100.0 / _num_filled);
        }

        puts("========== collision cache miss ========= ");
        for (uint32_t i = 0; i < sizeof(steps) / sizeof(steps[0]); i++) {
            sumc += steps[i];
            if (steps[i] <= 2)
                continue;
            //            printf("  %2d  %8d  %.2lf  %.2lf\n", i, steps[i], steps[i] * 100.0 / collision, sumc * 100.0 / collision);
        }

        if (sumb == 0)  return;
        printf("    _num_filled/bucket_size/packed collision/cache_miss/hit_find = %u/%.2lf/%zd/ %.2lf%%/%.2lf%%/%.2lf\n",
                _num_filled, _num_filled * 1.0 / sumb, sizeof(PairT), (collision * 100.0 / _num_filled), (collision - steps[0]) * 100.0 / _num_filled, finds * 1.0 / _num_filled);
        assert(sumn == _num_filled);
        assert(sumc == collision);
    }
#endif

    /****
      std::pair<iterator, iterator> equal_range(const KeyT & key)
      {
      iterator found = find(key);
      if (found == end())
      return {found, found};
      else
      return {found, std::next(found)};
      }*/

    // ------------------------------------------------------------

    iterator find(const KeyT& key)
    {
        auto bucket = find_filled_bucket(key);
        if (bucket == INACTIVE) {
            return end();
        }
        return iterator(this, bucket);
    }

    const_iterator find(const KeyT& key) const
    {
        auto bucket = find_filled_bucket(key);
        if (bucket == INACTIVE) {
            return end();
        }
        return const_iterator(this, bucket);
    }

    bool contains(const KeyT& key) const
    {
        return find_filled_bucket(key) != INACTIVE;
    }

    size_t count(const KeyT& key) const
    {
        return find_filled_bucket(key) != INACTIVE ? 1 : 0;
    }

    /// Returns the matching ValueT or nullptr if k isn't found.
    bool try_get(const KeyT& key, ValueT& val)
    {
        auto bucket = find_filled_bucket(key);
        if (bucket != INACTIVE) {
            val = GET_VAL(_pairs, bucket);
            return true;
        }
        else {
            return false;
        }
    }

    /// Returns the matching ValueT or nullptr if k isn't found.
    ValueT* try_get(const KeyT& key)
    {
        auto bucket = find_filled_bucket(key);
        if (bucket != INACTIVE) {
            return &GET_VAL(_pairs, bucket);
        }
        else {
            return nullptr;
        }
    }

    /// Const version of the above
    const ValueT* try_get(const KeyT& key) const
    {
        auto bucket = find_filled_bucket(key);
        if (bucket != INACTIVE) {
            return &GET_VAL(_pairs, bucket);
        }
        else {
            return nullptr;
        }
    }

    /// Convenience function.
    const ValueT get_or_return_default(const KeyT& key) const
    {
        const ValueT* ret = try_get(key);
        if (ret) {
            return *ret;
        }
        else {
            return ValueT();
        }
    }

    // -----------------------------------------------------

    /// Returns a pair consisting of an iterator to the inserted element
    /// (or to the element that prevented the insertion)
    /// and a bool denoting whether the insertion took place.
    std::pair<iterator, bool> insert(const KeyT& key, const ValueT& value)
    {
        auto bucket = find_or_allocate(key);
        if (NEXT_BUCKET(_pairs, bucket) != INACTIVE) {
            return { iterator(this, bucket), false };
        }
        else {
            if (check_expand_need())
                bucket = insert_main_bucket(key, true);

            NEW_KVALUE(key, value, bucket);
            _num_filled++;
            return { iterator(this, bucket), true };
        }
    }

    std::pair<iterator, bool> insert(KeyT&& key, ValueT&& value)
    {
        auto bucket = find_or_allocate(key);
        if (NEXT_BUCKET(_pairs, bucket) != INACTIVE) {
            return { iterator(this, bucket), false };
        }
        else {
            if (check_expand_need())
                bucket = insert_main_bucket(key, true);

            NEW_KVALUE(std::move(key), std::move(value), bucket);
            _num_filled++;
            return { iterator(this, bucket), true };
        }
    }

    inline std::pair<iterator, bool> insert(const std::pair<KeyT, ValueT>& p)
    {
        return insert(p.first, p.second);
    }

    inline std::pair<iterator, bool> insert(std::pair<KeyT, ValueT>&& p)
    {
        return insert(std::move(p.first), std::move(p.second));
    }

    inline void insert(const_iterator begin, const_iterator end)
    {
        for (; begin != end; ++begin) {
            insert(begin->first, begin->second);
        }
    }

    inline void insert_unique(const_iterator begin, const_iterator end)
    {
        for (; begin != end; ++begin) {
            insert_unique(begin->first, begin->second);
        }
    }

    /// Same as above, but contains(key) MUST be false
    uint32_t insert_unique(const KeyT& key, const ValueT& value)
    {
        check_expand_need();
        auto bucket = insert_main_bucket(key, true);
        NEW_KVALUE(key, value, bucket);
        _num_filled++;
        return bucket;
    }

    uint32_t insert_unique(KeyT&& key, ValueT&& value)
    {
        check_expand_need();
        auto bucket = insert_main_bucket(key, true);
        NEW_KVALUE(std::move(key), std::move(value), bucket);
        _num_filled++;
        return bucket;
    }

    inline uint32_t insert_unique(std::pair<KeyT, ValueT>&& p)
    {
        return insert_unique(std::move(p.first), std::move(p.second));
    }

    inline uint32_t insert_unique(std::pair<KeyT, ValueT>& p)
    {
        return insert_unique(p.first, p.second);
    }

    //not
    template <class... Args>
        inline std::pair<iterator, bool> emplace(Args&&... args)
        {
            return insert(std::forward<Args>(args)...);
        }

    template <class... Args>
        inline std::pair<iterator, bool> emplace_unique(Args&&... args)
        {
            return insert_unique(std::forward<Args>(args)...);
        }

    uint32_t try_insert_mainbucket(const KeyT& key, const ValueT& value)
    {
        const auto bucket = hash_key(key);
        auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        if (next_bucket != INACTIVE)
            return INACTIVE;

        check_expand_need();
        NEW_KVALUE(key, value, bucket);
        _num_filled++;
        return bucket;
    }

    void insert_or_assign(const KeyT& key, ValueT&& value)
    {
        check_expand_need();

        auto bucket = find_or_allocate(key);
        // Check if inserting a new value rather than overwriting an old entry
        if (NEXT_BUCKET(_pairs, bucket) != INACTIVE) {
            GET_VAL(_pairs, bucket) = value;
        }
        else {
            NEW_KVALUE(key, value, bucket);
            _num_filled++;
        }
    }

    /// Return the old value or ValueT() if it didn't exist.
    ValueT set_get(const KeyT& key, const ValueT& new_value)
    {
        check_expand_need();

        auto bucket = find_or_allocate(key);

        // Check if inserting a new value rather than overwriting an old entry
        if (NEXT_BUCKET(_pairs, bucket) != INACTIVE) {
            ValueT old_value = GET_VAL(_pairs, bucket);
            GET_VAL(_pairs, bucket) = new_value;
            return old_value;
        }
        else {
            NEW_KVALUE(key, new_value, bucket);
            _num_filled++;
            return ValueT();
        }
    }

    /// Like std::map<KeyT,ValueT>::operator[].
    ValueT& operator[](const KeyT& key)
    {
        auto bucket = find_or_allocate(key);
        /* Check if inserting a new value rather than overwriting an old entry */
        if (NEXT_BUCKET(_pairs, bucket) == INACTIVE) {
            if (check_expand_need())
                bucket = insert_main_bucket(key, true);

            NEW_KVALUE(key, ValueT(), bucket);
            _num_filled++;
        }

        return GET_VAL(_pairs, bucket);
    }

    // -------------------------------------------------------

    /// Erase an element from the hash table.
    /// return false if element was not found
    bool erase(const KeyT& key)
    {
        auto bucket = erase_from_bucket(key);
        if (bucket == INACTIVE) {
            return false;
        }

        NEXT_BUCKET(_pairs, bucket) = INACTIVE;
        _pairs[bucket].~PairT();
        _num_filled -= 1;

#ifdef EMILIB_AUTO_SHRINK
        if (_num_buckets > 254 && _num_buckets > 4 * _num_filled)
            rehash(_num_filled / max_load_factor()  + 2);
#endif
        return true;
    }

    /// Erase an element typedef an iterator.
    /// Returns an iterator to the next element (or end()).
    iterator erase(iterator it)
    {
        auto bucket = erase_from_bucket(it->first);
        if (bucket == INACTIVE) {
            return ++it;
        }

        NEXT_BUCKET(_pairs, bucket) = INACTIVE;
        _pairs[bucket].~PairT();
        _num_filled -= 1;
        if (bucket == it._bucket)
            ++it;

#ifdef EMILIB_AUTO_SHRINK
        if (_num_buckets > 254 && _num_buckets > 4 * _num_filled) {
            rehash(_num_filled * max_load_factor() + 2);
            it = begin();
        }
#endif
        return it;
    }

    /// Remove all elements, keeping full capacity.
    void clear()
    {
        for (uint32_t bucket = 0; _num_filled > 0; ++bucket) {
            if (NEXT_BUCKET(_pairs, bucket) != INACTIVE) {
                NEXT_BUCKET(_pairs, bucket) = INACTIVE;
                _pairs[bucket].~PairT();
                _num_filled -= 1;
            }
        }
        _num_filled = 0;
    }

    /// Make room for this many elements
    bool reserve(uint32_t required_buckets)
    {
        if (required_buckets < _load_threshold)
            return false;
        else if (required_buckets < _num_filled)
            return false;

        rehash(required_buckets);
        return true;
    }

    /// Make room for this many elements
    void rehash(uint32_t required_buckets)
    {
        uint32_t num_buckets = 8;
        uint32_t shift = 3;
        while (num_buckets < required_buckets) { num_buckets *= 2;  shift ++;}
        if (num_buckets <= _num_buckets) {
            num_buckets = 2 * _num_buckets; shift ++;
        }

        assert(num_buckets * _max_load_factor + 2 >= _num_filled);
        auto new_pairs = (PairT*)malloc(num_buckets * sizeof(PairT));
/*
        if (!new_pairs) {
            throw std::bad_alloc();
        }
**/
        auto old_num_filled  = _num_filled;
        auto old_num_buckets = _num_buckets;
        auto old_pairs = _pairs;
        auto reset = 0;

        _num_filled  = 0;
        _num_buckets = num_buckets;
        _mask        = num_buckets - 1;
        _pairs       = new_pairs;

        for (uint32_t bucket = 0; bucket < num_buckets; bucket++)
            NEXT_BUCKET(_pairs, bucket) = INACTIVE;

        uint32_t collision = 0, mbucket = 0;
        //set all main bucket first
        for (uint32_t src_bucket = 0; src_bucket < old_num_buckets; src_bucket++) {
            if (NEXT_BUCKET(old_pairs, src_bucket) == INACTIVE)
                continue;

            const auto main_bucket = hash_key(GET_KEY(old_pairs, src_bucket));
            auto& next_bucket = NEXT_BUCKET(_pairs, main_bucket);
            auto& src_pair = old_pairs[src_bucket];
            if (next_bucket == INACTIVE) {
                new(_pairs + main_bucket) PairT(std::move(src_pair));
                next_bucket = main_bucket;
                mbucket ++;
            }
            else {
                //move collision bucket to head for better cache performance
                new(old_pairs + collision) PairT(std::move(src_pair));
                NEXT_BUCKET(old_pairs, collision++) = main_bucket;
            }
            src_pair.~PairT();
            _num_filled += 1;
            if (_num_filled >= old_num_filled)
                break;
        }

#if 1
        //reset all collisions bucket, not linke new bucket after main bucket beause of cache miss
        for (uint32_t src_bucket = 0; src_bucket < collision; src_bucket++) {
            auto& old_pair = old_pairs[src_bucket];
            const auto main_bucket = NEXT_BUCKET(old_pairs, src_bucket);
            auto& next_bucket = NEXT_BUCKET(_pairs, main_bucket);
            if (main_bucket == next_bucket)
            {
                const auto new_bucket = find_empty_bucket(main_bucket);
                new(_pairs + new_bucket) PairT(std::move(old_pair)); old_pair.~PairT();
                NEXT_BUCKET(_pairs, new_bucket) = next_bucket = new_bucket;
            }
            else
            {
                const auto last_bucket = find_last_bucket(next_bucket);//how to fast find the last bucket ?
                const auto new_bucket  = find_empty_bucket(last_bucket);
                new(_pairs + new_bucket) PairT(std::move(old_pair)); old_pair.~PairT();
                NEXT_BUCKET(_pairs, new_bucket) = NEXT_BUCKET(_pairs, last_bucket) = new_bucket;
                reset++;
            }
        }
#elif 0
        //reset all collisions bucket
        for (uint32_t bucket = 0; bucket < collision; bucket++) {
            auto new_bucket = insert_main_bucket(GET_KEY(old_pairs, bucket), false);
            auto& src_pair = old_pairs[bucket];
            new(_pairs + new_bucket) PairT(std::move(src_pair)); src_pair.~PairT();
            NEXT_BUCKET(_pairs, new_bucket) = new_bucket;
        }
#endif

#ifdef EMILIB_REHASH_LOG
        if (_num_filled > 0) {
            char buff[255] = {0};
            sprintf(buff, "    _num_filled/aver_size/K.V/pack/collision = %u/%2.lf/%s.%s/%zd/%.2lf%%",
                    _num_filled, double (_num_filled) / mbucket, typeid(KeyT).name(), typeid(ValueT).name(), sizeof(_pairs[0]), (collision * 100.0 / _num_filled));
#if EMILIB_TAF_LOG
            static uint32_t ihashs = 0;
            FDLOG() << "EMILIB_ORDER_INDEX = " << EMILIB_ORDER_INDEX << "|hash_nums = " << ihashs ++ << "|" <<__FUNCTION__ << "|" << buff << endl;
#else
            puts(buff);
#endif
        }
#endif

        _load_threshold = uint32_t(_num_buckets * max_load_factor());
        free(old_pairs);
        assert(old_num_filled == _num_filled);
    }

private:
    // Can we fit another element?
    inline bool check_expand_need()
    {
        return reserve(_num_filled);
    }

    uint32_t erase_from_bucket(const KeyT& key)
    {
        const auto bucket = hash_key(key);
        auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        if (next_bucket == INACTIVE)
            return INACTIVE;
        else if (next_bucket == bucket) {
            if (key == GET_KEY(_pairs, bucket))
                return bucket;
            return INACTIVE;
        }
        else if (key == GET_KEY(_pairs, bucket)) {
            const auto nbucket = NEXT_BUCKET(_pairs, next_bucket);
            GET_PVAL(_pairs, bucket).swap(GET_PVAL(_pairs, next_bucket));
            NEXT_BUCKET(_pairs, bucket) = (nbucket == next_bucket) ? bucket : nbucket;
            return next_bucket;
        }

        auto prev_bucket = bucket;
        while (true) {
            const auto nbucket = NEXT_BUCKET(_pairs, next_bucket);
            if (key == GET_KEY(_pairs, next_bucket)) {
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

    // Find the bucket with this key, or return INACTIVE
    uint32_t find_filled_bucket(const KeyT& key) const
    {
        const auto bucket = hash_key(key);
        auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        if (next_bucket == INACTIVE)
            return INACTIVE;
        else if (key == GET_KEY(_pairs, bucket))
            return bucket;
        else if (next_bucket == bucket)
            return INACTIVE;

        //find next linked bucket
#if EMILIB_LRU_FIND
        auto prev_bucket = bucket;
#endif
        while (true) {
            if (key == GET_KEY(_pairs, next_bucket)) {
#if EMILIB_LRU_FIND
                GET_PVAL(_pairs, next_bucket).swap(GET_PVAL(_pairs, prev_bucket));
                return prev_bucket;
#else
                return next_bucket;
#endif
            }
            const auto nbucket = NEXT_BUCKET(_pairs, next_bucket);
            if (nbucket == next_bucket)
                break;
#if EMILIB_LRU_FIND
            prev_bucket = next_bucket;
#endif
            next_bucket = nbucket;
        }

        return INACTIVE;
    }

    uint32_t reset_main_bucket(const uint32_t main_bucket, const uint32_t bucket)
    {
        const auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        const auto new_bucket  = find_empty_bucket(next_bucket);
        const auto prev_bucket = find_prev_bucket(main_bucket, bucket);
        NEXT_BUCKET(_pairs, prev_bucket) = new_bucket;
        new(_pairs + new_bucket) PairT(std::move(_pairs[bucket])); _pairs[bucket].~PairT();
        NEXT_BUCKET(_pairs, new_bucket) = (next_bucket == bucket) ? new_bucket : next_bucket;
        NEXT_BUCKET(_pairs, bucket) = INACTIVE;
        return new_bucket;
    }

    // Find the bucket with this key, or return a good empty bucket to place the key in.
    // In the latter case, the bucket is expected to be filled.
    // If the bucket opt by other Key who's main bucket is not in this postion, kick it out
    // and move it to a new empty postion.
    uint32_t find_or_allocate(const KeyT& key)
    {
        const auto bucket = hash_key(key);
        auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        const auto& bucket_key = GET_KEY(_pairs, bucket);
        if (next_bucket == INACTIVE || bucket_key == key)
            return bucket;
        else if (next_bucket == bucket && bucket == hash_key(bucket_key))
            return NEXT_BUCKET(_pairs, next_bucket) = find_empty_bucket(next_bucket);

        //find next linked bucket and check key
        while (true) {
            if (key == GET_KEY(_pairs, next_bucket)) {
#if EMILIB_LRU_SET
                GET_PVAL(_pairs, next_bucket).swap(GET_PVAL(_pairs, bucket));
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

        //check current bucket_key is in main bucket or not
        const auto main_bucket = hash_key(bucket_key);
        if (main_bucket != bucket) {
            reset_main_bucket(main_bucket, bucket);
            return bucket;
        }

        //find a new empty and link it to tail
        const auto new_bucket = find_empty_bucket(next_bucket);
        return NEXT_BUCKET(_pairs, next_bucket) = new_bucket;
    }

    // key is not in this map. Find a empty place to put it.
    // combine linear probing and quadratic probing
    uint32_t find_empty_bucket(uint32_t bucket_from) const
    {
        const auto bucket = (bucket_from + 1) & _mask;
        if (NEXT_BUCKET(_pairs, bucket) == INACTIVE)
            return bucket;

        auto slot = 2;
        constexpr bool small_kv = sizeof(PairT) * 4 < EMILIB_CACHE_LINE_SIZE * 3;
        if (small_kv) {
            //at most only one cache line miss for linear probing
            const auto bucket_address = (uint32_t)(reinterpret_cast<size_t>(&NEXT_BUCKET(_pairs, bucket_from)) % EMILIB_CACHE_LINE_SIZE);
            const auto line_probe_length = (uint32_t)((EMILIB_CACHE_LINE_SIZE * 2 - bucket_address) / sizeof(PairT));

            for (; slot <= line_probe_length; ++slot) {
                const auto bucket = (bucket_from + slot) & _mask;
                if (NEXT_BUCKET(_pairs, bucket) == INACTIVE)
                    return bucket;
            }
        }

        //use quadratic probing to accerate find speed for high load factor and clustering
        //if (line_probe_length > 6 && _num_filled * 16 > 11 * _mask) bucket_from += _num_buckets / 2;
        bucket_from += (slot * slot - slot) / 2 + 1;

        for ( ; ; ++slot) {
            const auto bucket1 = bucket_from & _mask;
            if (NEXT_BUCKET(_pairs, bucket1) == INACTIVE)
                return bucket1;
#if 0
            //check adjacent slot is empty and the slot is in the same cache line which can reduce cpu cache miss.
            const auto cache_offset = reinterpret_cast<size_t>(&NEXT_BUCKET(_pairs,bucket1)) % EMILIB_CACHE_LINE_SIZE;
            if (cache_offset + sizeof(PairT) < EMILIB_CACHE_LINE_SIZE) {
                const auto bucket2 = (bucket_from + 1) & _mask;
                if (NEXT_BUCKET(_pairs, bucket2) == INACTIVE)
                    return bucket2;
            }
#endif
            if (slot > 6)
                bucket_from += _num_buckets / 3, slot = 1;

            bucket_from += slot;
        }
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

    uint32_t insert_main_bucket(const KeyT& key, bool check_main)
    {
        const auto bucket = hash_key(key);
        auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        if (next_bucket == INACTIVE)
            return bucket;

        const auto& bucket_key = GET_KEY(_pairs, bucket);
        const auto main_bucket = hash_key(bucket_key);
        if (main_bucket == bucket) {
            if (next_bucket == bucket)
                return NEXT_BUCKET(_pairs, bucket) = find_empty_bucket(next_bucket);
        } else if (check_main) {
            //check current bucket_key is linked in main bucket
            reset_main_bucket(main_bucket, bucket);
            return bucket;
        }

        //find a new empty and linked it to tail
        uint32_t last_bucket = next_bucket;
        while (true) {
            const auto nbucket = NEXT_BUCKET(_pairs, next_bucket);
            if (nbucket == next_bucket) {
                last_bucket = nbucket;
                break;
            }
            next_bucket = nbucket;
        }

        return NEXT_BUCKET(_pairs, last_bucket) = find_empty_bucket(last_bucket);
    }

    //Thomas Wang's function
    //https://gist.github.com/badboy/6267743
    static inline uint32_t hash32(uint32_t key)
    {
#if 0
        return (key ^ 2166136261U)  * 16777619UL;
#else
        key += ~(key << 15);
        key ^= (key >> 10);
        key += (key << 3);
        key ^= (key >> 6);
        key += ~(key << 11);
        key ^= (key >> 16);
        return key;
#endif
    }

    //64 bit to 32 bit Hash Functions
    static inline uint32_t hash64(uint64_t key)
    {
#if 0
        return (key ^ 14695981039346656037ULL) * 1099511628211ULL;
#else
        key = (~key) + (key << 18); // key = (key << 18) - key - 1;
        key = key ^ (key >> 31);
        key = key * 21; // key = (key + (key << 2)) + (key << 4);
        key = key ^ (key >> 11);
        key = key + (key << 6);
        key = key ^ (key >> 22);
        return (uint32_t) key;
#endif
    }

    template<typename UType, typename std::enable_if<std::is_integral<UType>::value, long>::type = 0>
        //    template<class UType, class = typename std::enable_if<std::is_integral<UType>::value, int>::type>
    inline uint32_t hash_key(const UType key) const
    {
#ifdef EMILIB_FIBONACCI_HASH
        return (key * 2654435761ull) & _mask;
#elif EMILIB_HASH32
        if (sizeof(key) <= sizeof(uint32_t))
            return hash32(key) & _mask;
        else
            return hash64(key) & _mask;
#elif EMILIB_FAST_HASH
        return (uint32_t)key & _mask;
#else
        return _hasher(key) & _mask;
#endif
    }

    template<typename UType, typename std::enable_if<!std::is_integral<UType>::value, long>::type = 0>
    inline uint32_t hash_key(const UType& key) const
    {
#ifdef EMILIB_FIBONACCI_HASH
        return (_hasher(key) * 11400714819323198485ull) & _mask;
#else
        return _hasher(key) & _mask;
#endif
    }

private:
    HashT   _hasher;
    PairT*  _pairs;

    uint32_t  _num_buckets;
    uint32_t  _num_filled;
    uint32_t  _mask;  // _num_buckets minus one

    float     _max_load_factor;
    uint32_t  _load_threshold;
};

} // namespace emilib

#if __cplusplus > 199711
//template <class Key, class Val> using emihash = emilib1::HashMap<Key, Val, std::hash<Key>>;
#endif
