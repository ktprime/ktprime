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
    #undef  hash_key
    #undef  NEXT_BUCKET
    #undef  GET_KEY
#endif

// likely/unlikely
#if (__GNUC__ >= 4 || __clang__)
#    define EMILIB_LIKELY(condition) __builtin_expect(condition, 1)
#    define EMILIB_UNLIKELY(condition) __builtin_expect(condition, 0)
#else
#    define EMILIB_LIKELY(condition) condition
#    define EMILIB_UNLIKELY(condition) condition
#endif

#define NEW_KVALUE(key, bucket) new(_pairs + bucket) PairT(key, bucket)
#define hash_key(key)  (uint32_t)_hasher(key) & _mask

#if EMILIB_CACHE_LINE_SIZE < 32
    #define EMILIB_CACHE_LINE_SIZE 64
#endif

#define GET_KEY(p,n)     p[n].first
#define NEXT_BUCKET(s,n) s[n].second

namespace emilib5 {
/// A cache-friendly hash table with open addressing, linear probing and power-of-two capacity
template <typename KeyT, typename HashT = std::hash<KeyT>>
class HashSet
{
    constexpr static uint32_t INACTIVE = 0xFFFFFFFF;

private:
    typedef  HashSet<KeyT, HashT> MyType;
    typedef  std::pair<KeyT, int> PairT;

public:
    typedef size_t   size_type;
    typedef KeyT     value_type;
    typedef KeyT&    reference;
    typedef const KeyT& const_reference;

    class iterator
    {
    public:
        typedef std::forward_iterator_tag iterator_category;
        typedef size_t                    difference_type;
        typedef size_t                    distance_type;
        typedef KeyT                      value_type;
        typedef value_type*               pointer;
        typedef value_type&               reference;

        iterator() { }

        iterator(MyType* hash_set, uint32_t bucket) : _set(hash_set), _bucket(bucket)
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
            return iterator(_set, old_index);
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
            do {
                _bucket++;
            } while (_bucket < _set->_num_buckets && _set->NEXT_BUCKET(_pairs, _bucket) == INACTIVE);
        }

    public:
        MyType* _set;
        uint32_t  _bucket;
    };

    class const_iterator
    {
    public:
        typedef std::forward_iterator_tag iterator_category;
        typedef size_t                    difference_type;
        typedef size_t                    distance_type;
        typedef const KeyT                value_type;
        typedef value_type*               pointer;
        typedef value_type&               reference;

        const_iterator() { }

        const_iterator(iterator proto) : _set(proto._set), _bucket(proto._bucket)
        {
        }

        const_iterator(const MyType* hash_set, uint32_t bucket) : _set(hash_set), _bucket(bucket)
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
            return const_iterator(_set, old_index);
        }

        reference operator*() const
        {
            return _set->GET_KEY(_pairs, _bucket);
        }

        pointer operator->() const
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
            do {
                _bucket++;
            } while (_bucket < _set->_num_buckets && _set->NEXT_BUCKET(_pairs, _bucket) == INACTIVE);
        }

    public:
        const MyType* _set;
        uint32_t  _bucket;
    };

    // ------------------------------------------------------------------------

    void init()
    {
        _num_buckets = 0;
        _num_filled = 0;
        _mask = 0;
        _pairs = nullptr;
        max_load_factor(0.9f);
    }

    HashSet(uint32_t bucket = 4)
    {
        init();
        reserve(bucket);
    }

    HashSet(const HashSet& other)
    {
        _pairs = (PairT*)malloc(other._num_buckets * sizeof(PairT));
        clone(other);
    }

    void clone(const HashSet& other)
    {
        _hasher      = other._hasher;
        _num_buckets = other._num_buckets;
        _num_filled  = other._num_filled;
        _mask        = other._mask;
        _loadlf      = other._loadlf;

        if (std::is_integral<KeyT>::value) {
            memcpy(_pairs, other._pairs, other._num_buckets * sizeof(PairT));
        }
        else {
            auto old_pairs = other._pairs;
            for (uint32_t bucket = 0; bucket < _num_buckets; bucket++) {
                auto state = NEXT_BUCKET(_pairs, bucket) = NEXT_BUCKET(old_pairs, bucket);
                if (state != INACTIVE)
                    new(_pairs + bucket) PairT(old_pairs[bucket]);
            }
        }
    }

    HashSet(HashSet&& other)
    {
        init();
        reserve(1);
        *this = std::move(other);
    }

    HashSet(std::initializer_list<KeyT> il)
    {
        init();
        reserve((uint32_t)il.size());
        for (auto begin = il.begin(); begin != il.end(); ++begin)
            insert(*begin);
    }

    HashSet& operator=(const HashSet& other)
    {
        if (this == &other)
            return *this;

        clear();
        if (_num_buckets < other._num_buckets) {
            free(_pairs);
            _pairs = (PairT*)malloc(other._num_buckets * sizeof(PairT));
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
        clear();
        free(_pairs);
    }

    void swap(HashSet& other)
    {
        std::swap(_hasher, other._hasher);
        std::swap(_pairs, other._pairs);
        std::swap(_num_buckets, other._num_buckets);
        std::swap(_num_filled, other._num_filled);
        std::swap(_mask, other._mask);
        std::swap(_loadlf, other._loadlf);
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
        return  (1 << 20) / _loadlf;
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

    int get_bucket_info(const uint32_t bucket, uint32_t steps[], const uint32_t slots) const
    {
        auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        if (next_bucket == INACTIVE)
            return -1;

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
            //printf("  %2d  %8d  %.2lf  %.2lf\n", i, steps[i], steps[i] * 100.0 / collision, sumc * 100.0 / collision);
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
        return find_filled_bucket(key) == INACTIVE ? 0 : 1;
    }

    /// Returns a pair consisting of an iterator to the inserted element
    /// (or to the element that prevented the insertion)
    /// and a bool denoting whether the insertion took place.
    std::pair<iterator, bool> insert(const KeyT& key)
    {
        auto bucket = find_or_allocate(key);
        if (NEXT_BUCKET(_pairs, bucket) != INACTIVE) {
            return { iterator(this, bucket), false };
        }
        else {
            if (EMILIB_UNLIKELY(check_expand_need()))
                bucket = insert_main_bucket(key);

            NEW_KVALUE(key, bucket);
            _num_filled++;
            return { iterator(this, bucket), true };
        }
    }

    std::pair<iterator, bool> insert(KeyT&& key)
    {
        auto bucket = find_or_allocate(key);
        if (NEXT_BUCKET(_pairs, bucket) != INACTIVE) {
            return { iterator(this, bucket), false };
        }
        else {
            if (check_expand_need())
                bucket = insert_main_bucket(key);

            NEW_KVALUE(std::move(key), bucket);
            _num_filled++;
            return { iterator(this, bucket), true };
        }
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
    uint32_t insert_unique(const KeyT& key)
    {
        check_expand_need();
        auto bucket = insert_main_bucket(key);
        NEW_KVALUE(key, bucket);
        _num_filled++;
        return bucket;
    }

    uint32_t insert_unique(KeyT&& key)
    {
        check_expand_need();
        auto bucket = insert_main_bucket(key);
        NEW_KVALUE(std::move(key), bucket);
        _num_filled++;
        return bucket;
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


    void insert_or_assign(const KeyT& key)
    {
        check_expand_need();

        auto bucket = find_or_allocate(key);
        // Check if inserting a new value rather than overwriting an old entry
        if (NEXT_BUCKET(_pairs, bucket) != INACTIVE) {
            NEXT_BUCKET(_pairs, bucket) = bucket;
        }
        else {
            NEW_KVALUE(key, bucket);
            _num_filled++;
        }
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

        return it;
    }

    /// Remove all elements, keeping full capacity.
    void clear()
    {
        if (std::is_integral<KeyT>::value) {
            _num_filled = 0;
            memset(_pairs, INACTIVE, sizeof(_pairs[0]) * _num_buckets);
            return;
        }

        for (uint32_t bucket = 0; _num_filled > 0; ++bucket) {
            if (NEXT_BUCKET(_pairs, bucket) != INACTIVE) {
                NEXT_BUCKET(_pairs, bucket) = INACTIVE;
                _pairs[bucket].~PairT();
                _num_filled -= 1;
            }
        }
    }

    /// Make room for this many elements
    bool reserve(uint32_t num_elems)
    {
        //auto required_buckets = (uint32_t)(((size_t)num_elems * _loadlf) >> 20) + 2;
        auto required_buckets = num_elems * 10 / 8 + 2;
        if (EMILIB_LIKELY(required_buckets <= _num_buckets))
            return false;

        rehash(required_buckets);
        return true;
    }

    /// Make room for this many elements
    void rehash(uint32_t required_buckets)
    {
        uint32_t num_buckets = 4;
        while (num_buckets < required_buckets) { num_buckets *= 2; }

        //assert(num_buckets > _num_filled);
        auto new_pairs = (PairT*)malloc(num_buckets * sizeof(PairT));
        auto old_num_filled  = _num_filled;
        auto old_num_buckets = _num_buckets;
        auto old_pairs = _pairs;
        auto reset = 0;

        _num_filled  = 0;
        _num_buckets = num_buckets;
        _mask        = num_buckets - 1;
        _pairs       = new_pairs;

        if (sizeof(PairT) <= 64)
            memset(_pairs, INACTIVE, sizeof(_pairs[0]) * num_buckets);
        else
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
#if 1
                new(old_pairs + collision) PairT(std::move(src_pair));
                NEXT_BUCKET(old_pairs, collision++) = main_bucket;
#else
                NEXT_BUCKET(old_pairs, collision++) = src_bucket;
#endif
            }
            src_pair.~PairT();
            _num_filled += 1;
            if (EMILIB_UNLIKELY(_num_filled >= old_num_filled))
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
            auto src_bucket = NEXT_BUCKET(old_pairs, bucket);
            auto& src_pair = old_pairs[src_bucket];
            auto new_bucket = insert_main_bucket(GET_KEY(old_pairs, src_bucket));
            new(_pairs + new_bucket) PairT(std::move(src_pair)); src_pair.~PairT();
            NEXT_BUCKET(_pairs, new_bucket) = new_bucket;
        }
#endif

#if EMILIB_REHASH_LOG
        if (_num_filled > 100000) {
            char buff[255] = {0};
            sprintf(buff, "    _num_filled/aver_size/K.V/pack/collision = %u/%2.lf/%s/%zd/%.2lf%%",
                    _num_filled, double (_num_filled) / mbucket, typeid(KeyT).name(), sizeof(_pairs[0]), (collision * 100.0 / _num_filled));
#if EMILIB_TAF_LOG
            static uint32_t ihashs = 0;
            FDLOG() << "|hash_nums = " << ihashs ++ << "|" <<__FUNCTION__ << "|" << buff << endl;
#else
            puts(buff);
#endif
        }
#endif
        if (old_pairs)
        free(old_pairs);
//      assert(old_num_filled == _num_filled);
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

        const auto bqKey = key == GET_KEY(_pairs, bucket);
        if (next_bucket == bucket)
            return bqKey ? bucket : INACTIVE;
        else if (bqKey) {
            const auto nbucket = NEXT_BUCKET(_pairs, next_bucket);
            std::swap(GET_KEY(_pairs, bucket), GET_KEY(_pairs, next_bucket));
            NEXT_BUCKET(_pairs, bucket) = (nbucket == next_bucket) ? bucket : nbucket;
            return next_bucket;
        }
        else if (bucket != (hash_key(GET_KEY(_pairs, bucket))))
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

    // Find the bucket with this key, or return INACTIVE
    uint32_t find_filled_bucket(const KeyT& key) const
    {
        const auto bucket = hash_key(key);
        auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        const auto& bucket_key = GET_KEY(_pairs, bucket);
#if 0
        if (key == bucket_key && next_bucket != INACTIVE)
            return bucket;
       else if (next_bucket == bucket || next_bucket == INACTIVE)
            return INACTIVE;
        else if (EMILIB_UNLIKELY(bucket != (hash_key(bucket_key))))
            return INACTIVE;
#else
        if (next_bucket == INACTIVE)
            return INACTIVE;
        else if (key == bucket_key)
            return bucket;
        else if (next_bucket == bucket)
            return INACTIVE;
        else if (EMILIB_UNLIKELY(bucket != (hash_key(bucket_key))))
            return INACTIVE;
#endif

        //find next linked bucket
#if EMILIB_LRU_FIND
        auto prev_bucket = bucket;
#endif
        while (true) {
            if (key == GET_KEY(_pairs, next_bucket)) {
#if EMILIB_LRU_FIND
                GET_KEY(_pairs, next_bucket).swap(GET_KEY(_pairs, prev_bucket));
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
    // If the bucket opt by other key who's main bucket is not in this postion, kick it out
    // and move it to a new empty postion.
    uint32_t find_or_allocate(const KeyT& key)
    {
        const auto bucket = hash_key(key);

        const auto& bucket_key = GET_KEY(_pairs, bucket);
        auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        if (next_bucket == INACTIVE || key == bucket_key)
            return bucket;

        //check current bucket_key is in main bucket or not
        const auto main_bucket = hash_key(bucket_key);
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
    uint32_t find_empty_bucket(uint32_t bucket_from) const
    {
        const auto bucket = (++bucket_from) & _mask;
        if (NEXT_BUCKET(_pairs, bucket) == INACTIVE)
            return bucket;
#if 1
        const auto bucket_address = (uint32_t)(reinterpret_cast<size_t>(&NEXT_BUCKET(_pairs, bucket_from)) % EMILIB_CACHE_LINE_SIZE);
        const auto max_probe_length = (uint32_t)((EMILIB_CACHE_LINE_SIZE * 2 - bucket_address + sizeof(int)) / sizeof(PairT));
#else
        constexpr auto max_probe_length = (EMILIB_CACHE_LINE_SIZE * 1 / sizeof(PairT)) + 1;//cpu cache line 64 byte,2-3 cache line miss
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
                else if (slot > 6 || max_probe_length > 5)
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

    uint32_t insert_main_bucket(const KeyT& key)
    {
        const auto bucket = hash_key(key);
        const auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        if (next_bucket == INACTIVE)
            return bucket;

        const auto main_bucket = hash_key(GET_KEY(_pairs, bucket));
        const auto last_bucket = find_last_bucket(next_bucket);
        const auto new_bucket  = find_empty_bucket(last_bucket);
        if (main_bucket == bucket) {
            return NEXT_BUCKET(_pairs, last_bucket) = new_bucket;
        }
        else {
            const auto prev_bucket = find_prev_bucket(main_bucket, bucket);
            NEXT_BUCKET(_pairs, prev_bucket) = new_bucket;
            new(_pairs + new_bucket) PairT(std::move(_pairs[bucket])); _pairs[bucket].~PairT();
            NEXT_BUCKET(_pairs, new_bucket) = (next_bucket == bucket) ? new_bucket : next_bucket;
            NEXT_BUCKET(_pairs, bucket) = INACTIVE;
            return bucket;
        }
    }

private:

    HashT     _hasher;
    PairT*    _pairs;
    uint32_t  _num_buckets;
    uint32_t  _mask;
    uint32_t  _loadlf;

    uint32_t  _num_filled;
};

} // namespace emilib
#if __cplusplus > 199711
//template <class Key, class Val> using emihash = emilib1::HashSet<Key, Val, std::hash<Key>>;
#endif

#undef hash_key
