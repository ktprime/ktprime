// By Huang Yuanbing 2019
// bailuzhou@163.com
// https://github.com/ktprime/ktprime/blob/master/hash_table5.hpp

// LICENSE:
//   This software is dual-licensed to the public domain and under the following
//   license: you are granted a perpetual, irrevocable license to copy, modify,
//   publish, and distribute this file as you see fit.


#pragma once

#include <cstring>
#include <cassert>

#if EMILIB_TAF_LOG
    #include "servant/AutoLog.h"
    #include "servant/RollLogHelper.h"
#endif

#ifdef  GET_KEY
    #undef  GET_KEY
    #undef  GET_VAL
    #undef  NEXT_BUCKET
    #undef  GET_PVAL
    #undef  hash_bucket
#endif

// likely/unlikely
#if (__GNUC__ >= 4 || __clang__)
#    define EMILIB_LIKELY(condition) __builtin_expect(condition, 1)
#    define EMILIB_UNLIKELY(condition) __builtin_expect(condition, 0)
#else
#    define EMILIB_LIKELY(condition) condition
#    define EMILIB_UNLIKELY(condition) condition
#endif

#define hash_bucket(key)  ((uint32_t)_hasher(key) & _mask)

#ifndef EMILIB_BUCKET_INDEX
    #define EMILIB_BUCKET_INDEX 1
#endif
#if EMILIB_CACHE_LINE_SIZE < 32
    #define EMILIB_CACHE_LINE_SIZE 64
#endif

#if EMILIB_BUCKET_INDEX == 0
    #define GET_KEY(p,n)     p[n].second.first
    #define GET_VAL(p,n)     p[n].second.second
    #define NEXT_BUCKET(s,n) s[n].first
    #define GET_PVAL(s,n)    s[n].second
    #define NEW_KVALUE(key, value, bucket) new(_pairs + bucket) PairT(bucket, std::pair<KeyT, ValueT>(key, value))
#elif EMILIB_BUCKET_INDEX == 2
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

constexpr static uint32_t INACTIVE = 0xFFFFFFFF;

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
        nextbucket = INACTIVE;
    }

    myPair(std::pair<First, Second>&& pair)
        :first(std::move(pair.first)), second(std::move(pair.second))
    {
        nextbucket = INACTIVE;
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

    Second second;//int
    uint32_t nextbucket;
    First first; //long
};// __attribute__ ((packed));

/// A cache-friendly hash table with open addressing, linear probing and power-of-two capacity
template <typename KeyT, typename ValueT, typename HashT = std::hash<KeyT>>
class HashMap
{

private:
    typedef  HashMap<KeyT, ValueT, HashT> MyType;

#if EMILIB_BUCKET_INDEX == 0
    typedef std::pair<uint32_t, std::pair<KeyT, ValueT>> PairT;
#elif EMILIB_BUCKET_INDEX == 2
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

#if EMILIB_BUCKET_INDEX == 1
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
            return {_map, old_index};
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
#if EMILIB_BUCKET_INDEX == 1
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
            return {_map, old_index};
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
        _mask = 0;
        _pairs = nullptr;
        _num_filled = 0;
        max_load_factor(0.9f);
    }

    HashMap(uint32_t bucket = 4)
    {
        init();
        reserve(bucket);
    }

    HashMap(const HashMap& other)
    {
        _pairs = (PairT*)malloc(other._num_buckets * sizeof(PairT));
        clone(other);
    }

    void clone(const HashMap& other)
    {
        _hasher      = other._hasher;
        _num_buckets = other._num_buckets;
        _num_filled  = other._num_filled;
        _mask        = other._mask;
        _loadlf      = other._loadlf;

        if (std::is_pod<KeyT>::value && std::is_pod<ValueT>::value) {
            memcpy(_pairs, other._pairs, other._num_buckets * sizeof(PairT));
        }
        else {
            auto old_pairs = other._pairs;
            for (uint32_t bucket = 0; bucket < _num_buckets; bucket++) {
                auto next_bucket = NEXT_BUCKET(_pairs, bucket) = NEXT_BUCKET(old_pairs, bucket);
                if (next_bucket != INACTIVE)
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
        reserve((uint32_t)il.size());
        for (auto begin = il.begin(); begin != il.end(); ++begin)
            insert(*begin);
    }

    HashMap& operator=(const HashMap& other)
    {
        if (this == &other)
            return *this;

        if (!std::is_pod<KeyT>::value || !std::is_pod<ValueT>::value)
            clearkv();

        if (_num_buckets < other._num_buckets) {
            free(_pairs);
            _pairs = (PairT*)malloc(other._num_buckets * sizeof(PairT));
        }

        clone(other);
        return *this;
    }

    HashMap& operator=(HashMap&& other)
    {
        this->swap(other);
        return *this;
    }

    ~HashMap()
    {
        if (!std::is_pod<KeyT>::value || !std::is_pod<ValueT>::value)
            clearkv();

        free(_pairs);
    }

    void swap(HashMap& other)
    {
        std::swap(_hasher, other._hasher);
        std::swap(_pairs, other._pairs);
        std::swap(_num_buckets, other._num_buckets);
        std::swap(_num_filled, other._num_filled);
        std::swap(_mask, other._mask);
        std::swap(_loadlf, other._loadlf);
    }

    // -------------------------------------------------------------

    iterator begin() const
    {
        uint32_t bucket = 0;
        while (bucket < _num_buckets && NEXT_BUCKET(_pairs, bucket) == INACTIVE) {
            ++bucket;
        }
        return {this, bucket};
    }

    const_iterator cbegin() const
    {
        uint32_t bucket = 0;
        while (bucket < _num_buckets && NEXT_BUCKET(_pairs, bucket) == INACTIVE) {
            ++bucket;
        }
        return {this, bucket};
    }

    const_iterator begin()
    {
        return cbegin();
    }

    iterator end() const
    {
        return {this, _num_buckets};
    }

    const_iterator cend() const
    {
        return {this, _num_buckets};
    }

    const_iterator end()
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
        return static_cast<float>(_num_filled) / static_cast<float>(_num_buckets);
    }

    HashT hash_function() const
    {
        return _hasher;
    }

    constexpr float max_load_factor() const
    {
        return (1 << 20) / _loadlf;
    }

    void max_load_factor(float value)
    {
        if (value < 0.95 && value > 0.2)
            _loadlf = (uint32_t)((1 << 20) / value);
    }

    constexpr size_type max_size() const
    {
        return (1 << 30) / sizeof(PairT);
    }

    constexpr size_type max_bucket_count() const
    {
        return (1 << 30) / sizeof(PairT);
    }

    //Returns the bucket number where the element with key k is located.
    size_type bucket(const KeyT& key) const
    {
        const auto bucket = hash_bucket(key);
        const auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        if (next_bucket == INACTIVE)
            return 0;
        if (bucket == next_bucket)
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

#ifdef EMILIB_STATIS
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

    iterator find(const KeyT& key) const
    {
        return {this, find_filled_bucket(key)};
    }

    const_iterator find(const KeyT& key)
    {
        return {this, find_filled_bucket(key)};
    }

    bool contains(const KeyT& key) const
    {
        return find_filled_bucket(key) != _num_buckets;
    }

    size_type count(const KeyT& key) const
    {
        return find_filled_bucket(key) == _num_buckets ? 0 : 1;
    }

    /// Returns the matching ValueT or nullptr if k isn't found.
    bool try_get(const KeyT& key, ValueT& val) const
    {
        const auto bucket = find_filled_bucket(key);
        const auto find = bucket != _num_buckets;
        if (find) {
            val = GET_VAL(_pairs, bucket);
        }
        return find;
    }

    /// Returns the matching ValueT or nullptr if k isn't found.
    ValueT* try_get(const KeyT& key)
    {
        const auto bucket = find_filled_bucket(key);
        return bucket == _num_buckets ? nullptr : &GET_VAL(_pairs, bucket);
    }

    /// Const version of the above
    const ValueT* try_get(const KeyT& key) const
    {
        const auto bucket = find_filled_bucket(key);
        return bucket == _num_buckets ? nullptr : &GET_VAL(_pairs, bucket);
    }

    /// Convenience function.
    const ValueT get_or_return_default(const KeyT& key) const
    {
        const auto bucket = find_filled_bucket(key);
        return bucket == _num_buckets ? ValueT() : GET_VAL(_pairs, bucket);
    }

    // -----------------------------------------------------

    /// Returns a pair consisting of an iterator to the inserted element
    /// (or to the element that prevented the insertion)
    /// and a bool denoting whether the insertion took place.
    std::pair<iterator, bool> insert(const KeyT& key, const ValueT& value)
    {
        auto bucket = find_or_allocate(key);
        const auto find = NEXT_BUCKET(_pairs, bucket) == INACTIVE;
        if (find) {
            if (EMILIB_UNLIKELY(check_expand_need()))
                bucket = find_unique_bucket(key);

            NEW_KVALUE(key, value, bucket); _num_filled++;
        }
        return { {this, bucket}, find };
    }

    std::pair<iterator, bool> insert(KeyT&& key, ValueT&& value)
    {
        auto bucket = find_or_allocate(key);
        const auto find = NEXT_BUCKET(_pairs, bucket) == INACTIVE;
        if (find) {
            if (check_expand_need())
                bucket = find_unique_bucket(key);

            NEW_KVALUE(std::move(key), std::move(value), bucket); _num_filled++;
        }
        return { {this, bucket}, find };
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
        auto bucket = find_unique_bucket(key);
        NEW_KVALUE(key, value, bucket); _num_filled++;
        return bucket;
    }

    uint32_t insert_unique(KeyT&& key, ValueT&& value)
    {
        check_expand_need();
        auto bucket = find_unique_bucket(key);
        NEW_KVALUE(std::move(key), std::move(value), bucket); _num_filled++;
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
        const auto bucket = hash_bucket(key);
        auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        if (next_bucket != INACTIVE)
            return INACTIVE;

        check_expand_need();
        NEW_KVALUE(key, value, bucket); _num_filled++;
        return bucket;
    }

    void insert_or_assign(const KeyT& key, ValueT&& value)
    {
        check_expand_need();

        auto bucket = find_or_allocate(key);
        // Check if inserting a new value rather than overwriting an old entry
        if (NEXT_BUCKET(_pairs, bucket) == INACTIVE) {
            NEW_KVALUE(key, value, bucket); _num_filled++;
        }
        else {
            GET_VAL(_pairs, bucket) = value;
        }
    }

    /// Return the old value or ValueT() if it didn't exist.
    ValueT set_get(const KeyT& key, const ValueT& new_value)
    {
        check_expand_need();

        auto bucket = find_or_allocate(key);

        // Check if inserting a new value rather than overwriting an old entry
        if (NEXT_BUCKET(_pairs, bucket) == INACTIVE) {
            NEW_KVALUE(key, new_value, bucket); _num_filled++;
            return ValueT();
        }
        else {
            ValueT old_value = GET_VAL(_pairs, bucket);
            GET_VAL(_pairs, bucket) = new_value;
            return old_value;
        }
    }

    /// Like std::map<KeyT,ValueT>::operator[].
    ValueT& operator[](const KeyT& key)
    {
        auto bucket = find_or_allocate(key);
        /* Check if inserting a new value rather than overwriting an old entry */
        if (NEXT_BUCKET(_pairs, bucket) == INACTIVE) {
            if (EMILIB_UNLIKELY(check_expand_need()))
                bucket = find_unique_bucket(key);

            NEW_KVALUE(key, ValueT(), bucket); _num_filled++;
        }

        return GET_VAL(_pairs, bucket);
    }

    // -------------------------------------------------------
    //reset remove key which hash_value is not in current bucket
    bool reset_bucket_key(uint32_t mask_bucket)
    {
        if (std::is_integral<KeyT>::value)
        {
            GET_KEY(_pairs, mask_bucket) += 1;
            while (mask_bucket == hash_bucket(GET_KEY(_pairs, mask_bucket)))
                GET_KEY(_pairs, mask_bucket) += rand() + 1;
        }
        return true;
    }

    /// Erase an element from the hash table.
    /// return false if element was not found
    size_type erase(const KeyT& key)
    {
//         if (_num_filled == 0)
//            return 0;
        const auto bucket = erase_by_key(key);
        if (bucket == INACTIVE)
            return 0;

        NEXT_BUCKET(_pairs, bucket) = INACTIVE; _pairs[bucket].~PairT(); _num_filled -= 1;
        return 1;
    }

    /// Erase an element typedef an iterator.
    /// Returns an iterator to the next element (or end()).
    iterator erase(const_iterator it)
    {
#if 0
        // we assume that it always points to a valid entry, and not end().
        assert(this == it._map);
        if (it._bucket >= _num_buckets)
            return end();
        else if (INACTIVE == NEXT_BUCKET(_pairs, it._bucket)) {
            return ++it;
        }
#endif
        //assert(it->first == GET_KEY(_pairs, it._bucket));
        const auto bucket = erase_from_bucket(it._bucket);
        NEXT_BUCKET(_pairs, bucket) = INACTIVE; _pairs[bucket].~PairT(); _num_filled -= 1;
        //erase from main bucket, return main bucket as next
        if (bucket != it._bucket)
            return {this, it._bucket};

        iterator itnext = {this, it._bucket};
        return ++itnext;
    }

    void clearkv()
    {
        for (uint32_t bucket = 0; _num_filled > 0; ++bucket) {
            if (NEXT_BUCKET(_pairs, bucket) != INACTIVE) {
                NEXT_BUCKET(_pairs, bucket) = INACTIVE; _pairs[bucket].~PairT(); _num_filled -= 1;
            }
        }
    }

    /// Remove all elements, keeping full capacity.
    void clear()
    {
        if (_num_filled > _num_buckets / 4 && std::is_pod<KeyT>::value && std::is_pod<ValueT>::value)
            memset(_pairs, INACTIVE, sizeof(_pairs[0]) * _num_buckets);
        else
            clearkv();
        _num_filled = 0;
    }

    /// Make room for this many elements
    bool reserve(uint32_t num_elems)
    {
        //auto required_buckets = (uint32_t)(((uint64_t)num_elems * _loadlf) >> 20) + 2;
        const auto required_buckets = num_elems * 10 / 8 + 2;
        if (EMILIB_LIKELY(required_buckets <= _num_buckets))
            return false;

        rehash(required_buckets);
        return true;
    }

    /// Make room for this many elements
    void rehash(uint32_t required_buckets)
    {
        uint32_t num_buckets = 4;
        if (required_buckets >= 1024)
            num_buckets = 1024 * 2;
        while (num_buckets < required_buckets) { num_buckets *= 2; }

        //assert(num_buckets > _num_filled);
        auto new_pairs = (PairT*)malloc(num_buckets * sizeof(PairT));
        auto old_num_filled  = _num_filled;
        auto old_num_buckets = _num_buckets;
        auto old_pairs = _pairs;

        _num_filled  = 0;
        _num_buckets = num_buckets;
        _mask        = num_buckets - 1;
        _pairs       = new_pairs;

        if (sizeof(PairT) <= sizeof(int64_t) * 4)
            memset(_pairs, INACTIVE, sizeof(_pairs[0]) * num_buckets);
        else
            for (uint32_t bucket = 0; bucket < num_buckets; bucket++)
                NEXT_BUCKET(_pairs, bucket) = INACTIVE;

        uint32_t collision = 0;
        //set all main bucket first
        for (uint32_t src_bucket = 0; src_bucket < old_num_buckets; src_bucket++) {
            if (NEXT_BUCKET(old_pairs, src_bucket) == INACTIVE)
                continue;

            const auto main_bucket = hash_bucket(GET_KEY(old_pairs, src_bucket));
            auto& next_bucket = NEXT_BUCKET(_pairs, main_bucket);
            auto& old_pair = old_pairs[src_bucket];
            if (next_bucket == INACTIVE) {
                new(_pairs + main_bucket) PairT(std::move(old_pair)); old_pair.~PairT();
                next_bucket = main_bucket;
            }
            else {
                //move collision bucket to head for better cache performance
                NEXT_BUCKET(old_pairs, collision++) = src_bucket;
            }
            _num_filled += 1;
            if (EMILIB_UNLIKELY(_num_filled >= old_num_filled))
                break;
        }

        //reset all collisions bucket
        for (uint32_t colls = 0; colls < collision; colls++) {
            const auto src_bucket = NEXT_BUCKET(old_pairs, colls);
            const auto main_bucket = hash_bucket(GET_KEY(old_pairs, src_bucket));
            auto& old_pair = old_pairs[src_bucket];
            {
                auto next_bucket = NEXT_BUCKET(_pairs, main_bucket);
                //check current bucket_key is in main bucket or not
                if (next_bucket != main_bucket)
                    next_bucket = find_last_bucket(next_bucket);
                //find a new empty and link it to tail
                auto new_bucket = NEXT_BUCKET(_pairs, next_bucket) = find_empty_bucket(next_bucket);
                new(_pairs + new_bucket) PairT(std::move(old_pair)); old_pair.~PairT();
                NEXT_BUCKET(_pairs, new_bucket) = new_bucket;
            }
        }

#if EMILIB_REHASH_LOG
        if (_num_filled > 100000) {
            auto mbucket = _num_filled - collision;
            char buff[255] = {0};
            sprintf(buff, "    _num_filled/aver_size/K.V/pack/collision = %u/%2.lf/%s.%s/%zd/%.2lf%%",
                    _num_filled, double (_num_filled) / mbucket, typeid(KeyT).name(), typeid(ValueT).name(), sizeof(_pairs[0]), (collision * 100.0 / _num_filled));
#if EMILIB_TAF_LOG
            static uint32_t ihashs = 0;
            FDLOG() << "EMILIB_BUCKET_INDEX = " << EMILIB_BUCKET_INDEX << "|hash_nums = " << ihashs ++ << "|" <<__FUNCTION__ << "|" << buff << endl;
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

    uint32_t erase_by_key(const KeyT& key)
    {
        const auto bucket = hash_bucket(key);
        auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        if (next_bucket == INACTIVE)
            return INACTIVE;

        const auto eqkey = key == GET_KEY(_pairs, bucket);
        if (next_bucket == bucket)
            return eqkey ? bucket : INACTIVE;
        else if (eqkey) {
            const auto nbucket = NEXT_BUCKET(_pairs, next_bucket);
            GET_PVAL(_pairs, bucket).swap(GET_PVAL(_pairs, next_bucket));
            NEXT_BUCKET(_pairs, bucket) = (nbucket == next_bucket) ? bucket : nbucket;
            return next_bucket;
        }
        else if (EMILIB_UNLIKELY(bucket != hash_bucket(GET_KEY(_pairs, bucket))))
            return INACTIVE;

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

    uint32_t erase_from_bucket(const uint32_t bucket)
    {
        const auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        const auto main_bucket = hash_bucket(GET_KEY(_pairs, bucket));
        if (bucket == main_bucket) {
            //more than one bucket
            if (bucket != next_bucket) {
                const auto nbucket = NEXT_BUCKET(_pairs, next_bucket);
                GET_PVAL(_pairs, bucket).swap(GET_PVAL(_pairs, next_bucket));
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
#if 0
        const auto& bucket_key = GET_KEY(_pairs, bucket);
        if (key == bucket_key && next_bucket != INACTIVE)
            return bucket;
        else if (next_bucket == INACTIVE || next_bucket == bucket)
            return _num_buckets;
        else if (EMILIB_UNLIKELY(bucket != hash_bucket(bucket_key)))
            return _num_buckets;
#else
        if (next_bucket == INACTIVE)
            return _num_buckets;
        else if (key == GET_KEY(_pairs, bucket))
            return bucket;
        else if (next_bucket == bucket)
            return _num_buckets;
//        else if (hash_bucket(GET_KEY(_pairs, bucket)) != bucket)
//            return _num_buckets;
#endif

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

        return _num_buckets;
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
        if (next_bucket == INACTIVE || key == bucket_key)
            return bucket;

        //check current bucket_key is in main bucket or not
        const auto main_bucket = hash_bucket(bucket_key);
        if (main_bucket != bucket) {
            reset_main_bucket(main_bucket, bucket);
            return bucket;
        }
        else if (next_bucket == bucket)
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

        //find a new empty and link it to tail
        const auto new_bucket = find_empty_bucket(next_bucket);
        return NEXT_BUCKET(_pairs, next_bucket) = new_bucket;
    }

    // key is not in this map. Find a place to put it.
    uint32_t find_empty_bucket(uint32_t bucket_from) const
    {
        const auto bucket = (++bucket_from) & _mask;
        if (NEXT_BUCKET(_pairs, bucket) == INACTIVE)
            return bucket;
#if 1
        const auto bucket_address = (uint32_t)(reinterpret_cast<size_t>(&NEXT_BUCKET(_pairs, bucket_from)) % EMILIB_CACHE_LINE_SIZE);
        const auto max_probe_length = (uint32_t)((EMILIB_CACHE_LINE_SIZE * 2 - bucket_address + sizeof(int)) / sizeof(PairT));
#else
        constexpr auto max_probe_length = (EMILIB_CACHE_LINE_SIZE * 2 / sizeof(PairT)) + 2;//cpu cache line 64 byte,2-3 cache line miss
#endif
        for (uint32_t slot = 1; ; ++slot) {
            const auto bucket = (bucket_from + slot) & _mask;
            if (NEXT_BUCKET(_pairs, bucket) == INACTIVE)
                return bucket;
            else if (slot >= max_probe_length) {
                const auto bucket1 = (bucket + slot * slot) & _mask; //switch to square search
                if (NEXT_BUCKET(_pairs, bucket1) == INACTIVE)
                    return bucket1;
#if 0
                const auto cache_offset = (uint32_t)reinterpret_cast<size_t>(&NEXT_BUCKET(_pairs,bucket1)) % EMILIB_CACHE_LINE_SIZE;
                if (cache_offset + sizeof(PairT) < EMILIB_CACHE_LINE_SIZE) {
                    const auto bucket2 = (bucket_from + 1) & _mask;
                    if (NEXT_BUCKET(_pairs, bucket2) == INACTIVE)
                        return bucket2;
                }
#else
                const auto bucket2 = (bucket1 + 1) & _mask;
                if (NEXT_BUCKET(_pairs, bucket2) == INACTIVE)
                    return bucket2;
#endif
                else if (slot > 6 /*|| max_probe_length > 5*/)
                    bucket_from += _num_filled;
            }
        }
    }

    uint32_t find_last_bucket(uint32_t main_bucket) const
    {
        auto next_bucket = NEXT_BUCKET(_pairs, main_bucket);
        if (EMILIB_LIKELY(next_bucket == main_bucket))
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
        if (EMILIB_LIKELY(next_bucket == bucket))
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
        if (EMILIB_UNLIKELY(main_bucket != bucket)) {
            reset_main_bucket(main_bucket, bucket);
            return bucket;
        }
        else if (next_bucket != bucket)
            next_bucket = find_last_bucket(next_bucket);

        //find a new empty and link it to tail
        return NEXT_BUCKET(_pairs, next_bucket) = find_empty_bucket(next_bucket);
    }

private:

    //the first cache line packed
    HashT     _hasher;
    uint32_t  _loadlf;
    uint32_t  _num_buckets;
    uint32_t  _mask;
    //uint32_t  _pack[12];

    uint32_t  _num_filled;
    PairT*    _pairs;
};

} // namespace emilib
#if __cplusplus > 199711
//template <class Key, class Val> using emihash = emilib1::HashMap<Key, Val, std::hash<Key>>;
#endif

