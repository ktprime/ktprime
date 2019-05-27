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
    #define GET_KEY(p,n)     p[n]._mypair.first
    #define GET_VAL(p,n)     p[n]._mypair.second
    #define NEXT_BUCKET(s,n) s[n]._mypair._ibucket
    #define GET_PVAL(s,n)    s[n]._mypair
    #define NEW_KVALUE(key, value, bucket) new(_pairs + bucket) PairT(key, value, bucket)
#endif

namespace emilib1 {
template <typename First, typename Second>
struct pair {
    typedef First  first_type;
    typedef Second second_type;

    // pair constructors are explicit so we don't accidentally call this ctor when we don't have to.
    pair(const First& firstArg, const Second& secondArg)
        : first( firstArg )
        , second( secondArg ) { _ibucket = 0xFFFFFFFF; }

    pair(First&& firstArg, Second&& secondArg)
        : first( std::move(firstArg) )
        , second( std::move(secondArg) ) { _ibucket = 0xFFFFFFFF; }

    void swap(pair<First, Second>& o) {
        std::swap(first, o.first);
        std::swap(second, o.second);
    }

    Second second;//int
    uint32_t _ibucket;
    First first; //long
};// __attribute__ ((packed));

/// A cache-friendly hash table with open addressing, linear probing and power-of-two capacity
template <typename KeyT, typename ValueT, typename HashT = std::hash<KeyT>>
class HashMap
{
    constexpr static uint32_t INACTIVE = 0xFFFFFFFF;

private:
    typedef  HashMap<KeyT, ValueT, HashT> MyType;

#if EMILIB_BUCKET_INDEX == 0
    typedef std::pair<uint32_t, std::pair<KeyT, ValueT>> PairT;
#elif EMILIB_BUCKET_INDEX == 2
    typedef std::pair<std::pair<KeyT, ValueT>, uint32_t> PairT;
#else

    struct PairT
    {
        explicit PairT(const KeyT& key, const ValueT& value, uint32_t bucket)
            :_mypair(key, value)
        {
            _mypair._ibucket = bucket;
        }

        explicit PairT(PairT& pairT)
            :_mypair(pairT._mypair.first, pairT._mypair.second)
        {
            _mypair._ibucket = pairT._mypair._ibucket;
        }

        explicit PairT(PairT&& pairT) noexcept
            :_mypair(pairT._mypair.first, pairT._mypair.second)
        {
            _mypair._ibucket = pairT._mypair._ibucket;
        }

        pair<KeyT, ValueT> _mypair;
    };
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
        typedef pair<KeyT, ValueT>        value_type;
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
            } while (_bucket < _map->_num_buckets && (int)_map->NEXT_BUCKET(_pairs, _bucket) < 0);
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
        typedef pair<KeyT, ValueT>        value_type;
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
            } while (_bucket < _map->_num_buckets && (int)_map->NEXT_BUCKET(_pairs, _bucket) < 0);
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
        _moves = 0;
        _pairs = nullptr;
        _pempty = nullptr;
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
        _pempty = nullptr;
        clone(other);
    }

    void clone(const HashMap& other)
    {
        _hasher      = other._hasher;
        _num_buckets = other._num_buckets;
        _num_filled  = other._num_filled;
        _mask        = other._mask;
        _moves       = other._moves;
        _loadlf      = other._loadlf;

        if (std::is_integral<KeyT>::value && std::is_trivially_copyable<ValueT>::value) {
            memcpy(_pairs, other._pairs, other._num_buckets * sizeof(PairT));
        }
        else {
            auto old_pairs = other._pairs;
            for (uint32_t bucket = 0; bucket < _num_buckets; bucket++) {
                auto next_bucket = NEXT_BUCKET(_pairs, bucket) = NEXT_BUCKET(old_pairs, bucket);
                if ((int)next_bucket >= 0)
                    new(_pairs + bucket) PairT(old_pairs[bucket]);
            }
        }

        if (other._pempty) {
            _pempty = (uint32_t*)malloc(other._pempty[0] * sizeof(uint32_t));
            memcpy(_pempty, other._pempty, other._pempty[0] * sizeof(uint32_t));
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

        clear();
        if (_num_buckets < other._num_buckets) {
            free(_pairs);
            _pairs = (PairT*)malloc(other._num_buckets * sizeof(PairT));
        }

        if (_pempty) {
            free(_pempty); _pempty = nullptr;
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
        if (!std::is_pod<ValueT>::value || !std::is_pod<KeyT>::value)
        {
            for (uint32_t bucket = 0; _num_filled > 0; ++bucket) {
                auto& next_bucket = NEXT_BUCKET(_pairs, bucket);
                if ((int)next_bucket >= 0) {
                    _pairs[bucket].~PairT(); _num_filled -= 1;
                }
            }
        }

        free(_pairs);
        if (_pempty)
            free(_pempty);
    }

    void swap(HashMap& other)
    {
        std::swap(_hasher, other._hasher);
        std::swap(_pairs, other._pairs);
        std::swap(_num_buckets, other._num_buckets);
        std::swap(_num_filled, other._num_filled);
        std::swap(_mask, other._mask);
        std::swap(_loadlf, other._loadlf);
        std::swap(_moves, other._moves);
        std::swap(_pempty, other._pempty);
    }

    // -------------------------------------------------------------

    iterator begin() const
    {
        uint32_t bucket = 0;
        while (bucket < _num_buckets && (int)NEXT_BUCKET(_pairs, bucket) < 0) {
            ++bucket;
        }
        return {this, bucket};
    }

    const_iterator cbegin() const
    {
        uint32_t bucket = 0;
        while (bucket < _num_buckets && (int)NEXT_BUCKET(_pairs, bucket) < 0) {
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
        return  {this, _num_buckets};
    }

    const_iterator cend() const
    {
        return {this, _num_buckets};
    }

    const_iterator end()
    {
        return  {this, _num_buckets};
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
        return  (1 << 20) / _loadlf;
    }

    void max_load_factor(float value)
    {
        if (value < 0.99 && value > 0.2)
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
        if ((int)next_bucket < 0)
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
        if ((int)next_bucket < 0)
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
        if ((int)next_bucket < 0)
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
        if ((int)next_bucket < 0)
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
        auto bucket = find_filled_bucket(key);
        if (bucket == INACTIVE) {
            bucket = _num_buckets;
        }
        return {this, bucket};
    }

    const_iterator find(const KeyT& key)
    {
        auto bucket = find_filled_bucket(key);
        if (bucket == INACTIVE) {
            bucket = _num_buckets;
            //return const_iterator(this, _num_buckets);
        }
        return {this, bucket};
    }

    bool contains(const KeyT& key) const
    {
        return find_filled_bucket(key) != INACTIVE;
    }

    size_type count(const KeyT& key) const
    {
        return find_filled_bucket(key) == INACTIVE ? 0 : 1;
    }

    /// Returns the matching ValueT or nullptr if k isn't found.
    bool try_get(const KeyT& key, ValueT& val) const
    {
        const auto bucket = find_filled_bucket(key);
        const auto find = bucket != INACTIVE;
        if (find) {
            val = GET_VAL(_pairs, bucket);
        }
        return find;
    }

    /// Returns the matching ValueT or nullptr if k isn't found.
    ValueT* try_get(const KeyT& key)
    {
        const auto bucket = find_filled_bucket(key);
        if ((int)bucket < 0) {
            return nullptr;
        }
        else {
            return &GET_VAL(_pairs, bucket);
        }
    }

    /// Const version of the above
    const ValueT* try_get(const KeyT& key) const
    {
        const auto bucket = find_filled_bucket(key);
        if ((int)bucket < 0) {
            return nullptr;
        }
        else {
            return &GET_VAL(_pairs, bucket);
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
        const auto find = (int)NEXT_BUCKET(_pairs, bucket) < 0;
        if (find) {
            if (EMILIB_UNLIKELY(check_expand_need()))
                bucket = find_unique_bucket(key);

            NEW_KVALUE(key, value, bucket); _num_filled++;
        }
        return {{this, bucket}, find };
    }

    std::pair<iterator, bool> insert(KeyT&& key, ValueT&& value)
    {
        auto bucket = find_or_allocate(key);
        const auto find = (int)NEXT_BUCKET(_pairs, bucket) < 0;
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
        if ((int)next_bucket < 0)
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
        if ((int)NEXT_BUCKET(_pairs, bucket) < 0) {
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
        if ((int)NEXT_BUCKET(_pairs, bucket) < 0) {
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
        if ((int)NEXT_BUCKET(_pairs, bucket) < 0) {
            if (EMILIB_UNLIKELY(reserve(_num_filled)))
                bucket = find_unique_bucket(key);

            NEW_KVALUE(key, ValueT(), bucket); _num_filled++;
        }

        return GET_VAL(_pairs, bucket);
    }

    // -------------------------------------------------------
    //reset remove key which hash_value is not in current bucket
    bool reset_bucket_key(const uint32_t mask_bucket)
    {
        if (std::is_integral<KeyT>::value)
        {
            auto& key = GET_KEY(_pairs, mask_bucket);
            while (mask_bucket == hash_bucket(key)) {
                memset(&key, rand(), sizeof(key));
                key = GET_KEY(_pairs, mask_bucket);
            }
        }

        return true;
    }

    /// Erase an element from the hash table.
    /// return false if element was not found
    size_type erase(const KeyT& key)
    {
        //         if (_num_filled == 0)
        //            return 0;
        const auto bucket = erase_from_key(key);
        if (bucket == INACTIVE)
            return 0;

        NEXT_BUCKET(_pairs, bucket) = INACTIVE; _pairs[bucket].~PairT(); _num_filled -= 1;
//      push_pempty(bucket);
        return 1;
    }

    /// Erase an element typedef an iterator.
    /// Returns an iterator to the next element (or end()).
    //iterator erase(const iterator& it)
    iterator erase(const_iterator it)
    {
#if 0
        if (it._bucket >= _num_buckets)
            return end();
        else if (INACTIVE == NEXT_BUCKET(_pairs, it._bucket)) {
            assert(false);
            return ++it;
        }
#endif
        //assert(it->first == GET_KEY(_pairs, it._bucket));
        const auto bucket = erase_from_bucket(it._bucket);
        NEXT_BUCKET(_pairs, bucket) = INACTIVE; _pairs[bucket].~PairT(); _num_filled -= 1;
//      push_pempty(bucket);
        //erase from main bucket, return main bucket as next
        if (bucket != it._bucket)
            return {this, it._bucket};

        reset_bucket_key(bucket);
        iterator itnext = {this, it._bucket};
        return ++itnext;
    }

    /// Remove all elements, keeping full capacity.
    void clear()
    {
        _moves = 0;
        if (_pempty) {
            _pempty[1] = 0;
        }

        if (_num_filled > _num_buckets / 4 && std::is_integral<KeyT>::value && std::is_trivially_copyable<ValueT>::value) {
            _num_filled = 0;
            memset(_pairs, INACTIVE, sizeof(_pairs[0]) * _num_buckets);
            reset_bucket_key(hash_bucket(GET_KEY(_pairs, 0)));
            return;
        }

        for (uint32_t bucket = 0; _num_filled > 0; ++bucket) {
            auto& next_bucket = NEXT_BUCKET(_pairs, bucket);
            if ((int)next_bucket >= 0) {
                if (std::is_integral<KeyT>::value && bucket == hash_bucket(GET_KEY(_pairs, bucket)))
                    reset_bucket_key(bucket);
                next_bucket = INACTIVE; _pairs[bucket].~PairT(); _num_filled -= 1;
            }
        }
    }

   /***
     * 100         98                98
     * free_bucket free_size 1 2 .....n
   */
    void set_pempty(uint32_t free_buckets)
    {
        _pempty = (uint32_t*)malloc(free_buckets * sizeof(uint32_t) + 8);
        _pempty[0] = free_buckets;

        auto& empty_size = _pempty[1];
        empty_size = 0;

        for (uint32_t bucket = 0; bucket < _num_buckets && empty_size + 2 < free_buckets; ++bucket) {
            if ((int)NEXT_BUCKET(_pairs, bucket) < 0)
                _pempty[++empty_size + 1] = bucket;
        }
    }

    //TODO
    void push_pempty(uint32_t empty_bucket)
    {
#if EMILIB_HIGH_LOAD
        if (_pempty)
        {
            auto &empty_size = _pempty[1];
            if (_pempty[0] > empty_size + 2)
                _pempty[++empty_size + 1] = empty_bucket;
        }
#endif
    }

    uint32_t pop_pempty()
    {
        auto& empty_size = _pempty[1];
        for (; empty_size > 0; ) {
            const auto bucket = _pempty[--empty_size + 2];
            if ((int)NEXT_BUCKET(_pairs, bucket) < 0)
                return bucket;
        }

        empty_size = 0;
        for (uint32_t bucket = 0; bucket < _num_buckets && empty_size + 2 < _pempty[0]; ++bucket) {
            if ((int)NEXT_BUCKET(_pairs, bucket) < 0)
                _pempty[++empty_size + 1] = bucket;
        }

        return _pempty[--empty_size + 2];
    }

    /// Make room for this many elements
    bool reserve(uint32_t num_elems)
    {
        //auto required_buckets = (uint32_t)(((uint64_t)num_elems * _loadlf) >> 20) + 2;
        const auto required_buckets = num_elems * 10 / 8 + 2;
        if (EMILIB_LIKELY(required_buckets <= _num_buckets))
            return false;

#if EMILIB_HIGH_LOAD > 10000
        if (_num_filled > EMILIB_HIGH_LOAD) {
            const auto left = _num_buckets - _num_filled;
            if (!_pempty) {
                set_pempty(std::min(left + 2, _num_buckets * 2 / 10));
                return false;
            }
            else if (left > 1000) {
                return false;
            }
            free(_pempty); _pempty = nullptr;
        }
#endif

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
        auto new_pairs = (PairT*)malloc((num_buckets + 1) * sizeof(PairT));
        auto old_num_filled  = _num_filled;
        auto old_num_buckets = _num_buckets;
        auto old_pairs = _pairs;

        _num_filled  = 0;
        _num_buckets = num_buckets;
        _mask        = num_buckets - 1;
        _pairs       = new_pairs;

#if EMILIB_SAFE_HASH
        auto mbucket = 0;
        //adjust hash function if bad hash function, alloc more memory
        if (_moves == 0 && old_num_filled > 10000)
        {
            for (uint32_t src_bucket = 0; src_bucket < old_num_buckets; src_bucket++) {
                if (NEXT_BUCKET(old_pairs, src_bucket) == src_bucket)
                    mbucket ++;
            }
            if (mbucket * 4 < old_num_filled) {
                _moves = 1;
                while (mbucket < (old_num_filled >> _moves))
                    _moves ++;
            }
        }
#endif

        if (sizeof(PairT) <= 64) {
            memset(_pairs, INACTIVE, sizeof(_pairs[0]) * num_buckets);
            if (std::is_integral<KeyT>::value)
                reset_bucket_key(hash_bucket(GET_KEY(_pairs, 0)));
        }
        else {
            for (uint32_t bucket = 0; bucket < num_buckets; bucket++) {
                reset_bucket_key(bucket);
                NEXT_BUCKET(_pairs, bucket) = INACTIVE;
            }
        }

        uint32_t collision = 0;
        //set all main bucket first
        for (uint32_t src_bucket = 0; src_bucket < old_num_buckets; src_bucket++) {
            if ((int)NEXT_BUCKET(old_pairs, src_bucket) < 0)
                continue;

            const auto main_bucket = hash_bucket(GET_KEY(old_pairs, src_bucket));
            auto& next_bucket = NEXT_BUCKET(_pairs, main_bucket);
            auto& old_pair = old_pairs[src_bucket];
            if ((int)next_bucket < 0) {
                new(_pairs + main_bucket) PairT(std::move(old_pair)); old_pair.~PairT();
                next_bucket = main_bucket;
            }
            else {
                //move collision bucket to head for better cache performance
#if EMILIB_REHASH_COPY
                new(old_pairs + collision) PairT(std::move(old_pair)); old_pair.~PairT();
                NEXT_BUCKET(old_pairs, collision++) = main_bucket;
#else
                NEXT_BUCKET(old_pairs, collision++) = src_bucket;
#endif
            }
            _num_filled ++;
        }

        //reset all collisions bucket, not linke new bucket after main bucket beause of cache miss
        for (uint32_t colls = 0; colls < collision; colls++) {
#if EMILIB_REHASH_COPY
            const auto src_bucket = colls;
            const auto main_bucket = NEXT_BUCKET(old_pairs, src_bucket);
#else
            const auto src_bucket = NEXT_BUCKET(old_pairs, colls);
            const auto main_bucket = hash_bucket(GET_KEY(old_pairs, src_bucket));
#endif
            //
            auto& old_pair = old_pairs[src_bucket];
            {
                auto next_bucket = NEXT_BUCKET(_pairs, main_bucket);
                //check current bucket_key is in main bucket or not
#if 1
                if (next_bucket != main_bucket)
                    next_bucket = find_last_bucket(next_bucket);
                //find a new empty and link it to tail
                auto new_bucket = NEXT_BUCKET(_pairs, next_bucket) = find_empty_bucket(next_bucket);
                new(_pairs + new_bucket) PairT(std::move(old_pair)); old_pair.~PairT();
                NEXT_BUCKET(_pairs, new_bucket) = new_bucket;
#else
                auto new_bucket = find_empty_bucket(next_bucket);
                NEXT_BUCKET(_pairs, main_bucket) = new_bucket;
                new(_pairs + new_bucket) PairT(std::move(old_pair)); old_pair.~PairT();
                NEXT_BUCKET(_pairs, new_bucket)  = next_bucket == main_bucket ? new_bucket : next_bucket;
#endif
            }
        }

#if EMILIB_REHASH_LOG
        if (_num_filled > 10000) {
            auto mbucket = _num_filled - collision;
            char buff[255] = {0};
            sprintf(buff, "    _num_filled/_moves/aver_size/K.V/pack/collision = %u/%u/%.2lf/%s.%s/%zd/%.2lf%%",
                    _num_filled, _moves, (double)_num_filled / mbucket, typeid(KeyT).name(), typeid(ValueT).name(), sizeof(_pairs[0]), (collision * 100.0 / _num_filled));
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

    uint32_t erase_from_key(const KeyT& key)
    {
        const auto bucket = hash_bucket(key);
        auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        if ((int)next_bucket < 0)
            return INACTIVE;

        const auto bqKey = key == GET_KEY(_pairs, bucket);
        if (next_bucket == bucket) {
            if (!bqKey)
                return INACTIVE;
            reset_bucket_key(bucket);
            return bucket;
        }
        else if (bqKey) {
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

    // Find the bucket with this key, or return INACTIVE
    uint32_t find_filled_bucket(const KeyT& key) const
    {
        const auto bucket = hash_bucket(key);
        auto next_bucket = NEXT_BUCKET(_pairs, bucket);

        if (std::is_integral<KeyT>::value) {
            //no need check exist for performance
            if (key == GET_KEY(_pairs, bucket))
                return bucket;
            else if ((int)next_bucket < 0 || next_bucket == bucket)
                return INACTIVE;
        }
        else {
            if ((int)next_bucket < 0)
                return INACTIVE;
            else if (key == GET_KEY(_pairs, bucket))
                return bucket;
            else if (next_bucket == bucket)
                return INACTIVE;
        }

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
        if ((int)next_bucket < 0 || key == bucket_key)
            return bucket;

        //check current bucket_key is in main bucket or not
        const auto main_bucket = hash_bucket(bucket_key);
        if (EMILIB_UNLIKELY(main_bucket != bucket)) {
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
    uint32_t find_empty_bucket(uint32_t bucket_from)
    {
        const auto bucket = (++bucket_from) & _mask;
        if ((int)NEXT_BUCKET(_pairs, bucket) < 0)
            return bucket;

#if 1
        const auto bucket_address = (uint32_t)(reinterpret_cast<size_t>(&NEXT_BUCKET(_pairs, bucket_from)) % EMILIB_CACHE_LINE_SIZE);
        const auto max_probe_length = (uint32_t)((EMILIB_CACHE_LINE_SIZE * 2 - bucket_address + sizeof(int)) / sizeof(PairT));
#else
        constexpr auto max_probe_length = (EMILIB_CACHE_LINE_SIZE * 2 / sizeof(PairT)) + 2;//cpu cache line 64 byte,2-3 cache line miss
#endif
        for (uint32_t slot = 1; ; ++slot) {
            const auto bucket = (bucket_from + slot) & _mask;
            if ((int)NEXT_BUCKET(_pairs, bucket) < 0)
                return bucket;
            else if (slot >= max_probe_length) {
#if EMILIB_HIGH_LOAD
                if (_pempty)
                    return pop_pempty();
#endif
                const auto bucket1 = (bucket + slot * slot) & _mask; //switch to square search
                if ((int)NEXT_BUCKET(_pairs, bucket1) < 0)
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
                if (EMILIB_UNLIKELY((int)NEXT_BUCKET(_pairs, bucket2) < 0))
                    return bucket2;
#endif
                else if (slot > 6 /*|| max_probe_length > 5*/)
                    bucket_from += _num_buckets / 2;
            }
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

    uint32_t find_unique_bucket(const KeyT& key)
    {
        const auto bucket = hash_bucket(key);
        auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        if ((int)next_bucket < 0) {
            return bucket;
        }

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

    static inline uint32_t unhash(uint32_t key)
    {
#if 1
        uint64_t const r = key * UINT64_C(0xca4bcaa75ec3f625);
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
//       return (key >> 33 ^ key ^ key << 11);
#if 0
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

    template<typename UType, typename std::enable_if<std::is_integral<UType>::value, int>::type = 0>
    inline uint32_t hash_bucket(const UType key) const
    {
#if EMILIB_IDENTITY_HASH
        return ((uint32_t)key >> 3) & _mask;
#elif EMILIB_SAFE_HASH
        //if (_moves > 0) {
            if (sizeof(UType) <= sizeof(uint32_t))
                return unhash(key) & _mask;
            else
                return hash64(key) & _mask;
       // }
#else
        return _hasher(key) & _mask;
#endif
    }

    template<typename UType, typename std::enable_if<!std::is_integral<UType>::value, int>::type = 0>
    inline uint32_t hash_bucket(const UType& key) const
    {
#ifdef EMILIB_FIBONACCI_HASH
        return (_hasher(key) * 11400714819323198485ull) & _mask;
#else
        return _hasher(key) & _mask;
#endif
    }

    private:

    HashT     _hasher;
    PairT*    _pairs;
    uint32_t  _num_buckets;
    uint32_t  _mask;
    uint32_t  _loadlf;

    uint32_t  _num_filled;
    uint32_t  _moves;   //hash mask

    uint32_t*  _pempty;
};

} // namespace emilib

#if __cplusplus > 199711
//template <class Key, class Val> using emihash = emilib1::HashMap<Key, Val, std::hash<Key>>;
#endif

#undef EMILIB_BUCKET_INDEX
