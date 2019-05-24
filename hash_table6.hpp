
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
#include <climits>

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

#define EMILIB_BUCKET_INDEX 1
#if EMILIB_CACHE_LINE_SIZE < 32
    #define EMILIB_CACHE_LINE_SIZE 64
#endif

// likely/unlikely
#if (__GNUC__ >= 4 || __clang__)
#    define EMILIB_LIKELY(condition) __builtin_expect(condition, 1)
#    define EMILIB_UNLIKELY(condition) __builtin_expect(condition, 0)
#else
#    define EMILIB_LIKELY(condition) condition
#    define EMILIB_UNLIKELY(condition) condition
#endif

constexpr unsigned int EMILIB_BUCKET_NONE     = ~0u;
constexpr unsigned int EMILIB_HASH_BIT        = 4;
constexpr unsigned int EMILIB_HASH_MASK       = (1u << EMILIB_HASH_BIT) - 1; //0x000000FF;
constexpr unsigned int EMILIB_BUCKET_MASK     = EMILIB_BUCKET_NONE << EMILIB_HASH_BIT; //~EMILIB_HASH_MASK     //0xFFFFFF00

//static_assert(EMILIB_HASH_MASK == 0xff000000, "not queal");

#if 1
    #define GET_KEY(p,b)       p[b].first
    #define GET_VAL(p,b)       p[b].second

    #define ADDR_BUCKET(s,b)   s[b]._bucket
    #define MAPA_BUCKET(s,b)   s[b]._bucket >> EMILIB_HASH_BIT
    #define NEXT_BUCKET(s,b)   (s[b]._bucket >> EMILIB_HASH_BIT)
    #define CLS_BUCKET(s,b)    (s[b]._bucket |= EMILIB_BUCKET_MASK)
    #define RST_BUCKET(s,b)    s[b]._bucket = EMILIB_BUCKET_NONE;

   //o in [0, 7]

    #define GET_PVAL(p,n)     p[n]
    #define NEW_KVALUE(key, value, bucket) new(_pairs + bucket) PairT(key, value, bucket)
#endif

inline static uint32_t CTZ(const uint32_t n)
{
#if _MSC_VER > 1400
    unsigned long index;
    _BitScanForward(&index, n);
#elif __GNUC__
    uint32_t index = __builtin_ctz(n);
#elif ASM_X86
    stype index;
    #if __GNUC__ || __TINYC__
    __asm__ ("bsfl %1, %0\n" : "=r" (index) : "rm" (n) : "cc");
    #else
    __asm
    {
        bsf eax, n
        mov index, eax
    }
    #endif
#endif

    return (uint32_t)index;
}

namespace emilib3 {
/// like std::equal_to but no need to #include <functional>

template <typename First, typename Second>
struct myPair {
    myPair(const First& key, const Second& value, uint32_t bucket)
        :first(key), second(value)
    {
        _bucket = (_bucket & EMILIB_HASH_MASK) | (bucket << EMILIB_HASH_BIT);
    }

    myPair(const First&& key, const Second&& value, uint32_t bucket)
        :first(std::move(key)), second(std::move(value))
    {
        _bucket = (_bucket & EMILIB_HASH_MASK) | (bucket << EMILIB_HASH_BIT);
    }

    myPair(const std::pair<First,Second>& pair)
        :first(pair.first), second(pair.second)
    {
        _bucket = EMILIB_BUCKET_NONE;
    }

    myPair(std::pair<First, Second>&& pair)
        :first(std::move(pair.first)), second(std::move(pair.second))
    {
        _bucket = EMILIB_BUCKET_NONE;
    }

    myPair(const myPair& pairT)
        :first(pairT.first), second(pairT.second)
    {
        _bucket = (_bucket & EMILIB_HASH_MASK) | (pairT._bucket & EMILIB_BUCKET_MASK);
    }

    myPair(myPair&& pairT)
        :first(std::move(pairT.first)), second(std::move(pairT.second))
    {
        _bucket = (_bucket & EMILIB_HASH_MASK) | (pairT._bucket & EMILIB_BUCKET_MASK);
    }

    myPair& operator = (myPair&& pairT)
    {
        first = std::move(pairT.first);
        second = std::move(pairT.second);
        _bucket = pairT._bucket;
        return *this;
    }

    void swap(myPair<First, Second>& o)
    {
        std::swap(first, o.first);
        std::swap(second, o.second);
    }

    First  first;
    uint32_t _bucket;
    Second second;
};

/// A cache-friendly hash table with open addressing, linear probing and power-of-two capacity
template <typename KeyT, typename ValueT, typename HashT = std::hash<KeyT>>
class HashMap
{
    static constexpr unsigned INACTIVE = EMILIB_BUCKET_MASK >> EMILIB_HASH_BIT;

    inline void SET_BUCKET(uint32_t bucket, uint32_t next)
    {
        _pairs[bucket]._bucket = (_pairs[bucket]._bucket & EMILIB_HASH_MASK) | (next << EMILIB_HASH_BIT);
    }

    inline void SET_MHASH(uint32_t bucket, uint32_t hash_key)
    {
//        _pairs[bucket]._bucket &= (hash_key >> (32 - EMILIB_HASH_BIT)) | EMILIB_BUCKET_MASK;
    }

    bool TST_HASH(uint32_t bucket, uint32_t hash_key) const
    {
/**        hash_key >>= (32 - EMILIB_HASH_BIT);
        const auto hash_all = (_pairs[bucket]._bucket | hash_key) & EMILIB_HASH_MASK;
        if (hash_all == EMILIB_HASH_MASK)
            return false;
        return hash_all != hash_key;
*/
        return false;
    }

    inline void CLEAR_MHASH(uint32_t bucket)
    {
//        _pairs[bucket]._bucket |= EMILIB_HASH_MASK;
    }

    inline void SET_BIT(uint32_t bucket)
    {
        const auto main = bucket % EMILIB_HASH_BIT;
        const auto mask_bucket = bucket - main;
        _pairs[mask_bucket]._bucket &= ~(1 << main);
    }

    inline void CLR_BIT(uint32_t bucket)
    {
        const auto main = bucket % EMILIB_HASH_BIT;
        const auto mask_bucket = bucket - main;
        _pairs[mask_bucket]._bucket |= 1 << main;
    }

    inline uint32_t GET_BIT(uint32_t bucket) const
    {
        //assert(bucket % EMILIB_HASH_BIT == 0 && bucket < _mask);
        return _pairs[bucket]._bucket & EMILIB_HASH_MASK;
    }

    inline uint32_t TST_BIT(uint32_t bucket) const
    {
        auto bmask = _pairs[bucket]._bucket & EMILIB_HASH_MASK;
        if (bmask != 0)
            return bucket + CTZ(bmask);

        bucket = (bucket + EMILIB_HASH_BIT) & _mask;
        bmask = _pairs[bucket]._bucket & EMILIB_HASH_MASK;

        if (bmask > 0) {
            const auto empty_bucket = bucket + CTZ(bmask);
//            assert (NEXT_BUCKET(_pairs, empty_bucket) == INACTIVE);
            return empty_bucket;
        }

        return INACTIVE;
    }


//    inline void RST_BIT(uint32_t bucket) { _pairs[bucket]._bucket |= EMILIB_HASH_MASK; }

private:
    typedef  HashMap<KeyT, ValueT, HashT> MyType;

#if EMILIB_BUCKET_INDEX == 1
    typedef myPair<KeyT, ValueT>                    PairT;
#else
    typedef std::pair<std::pair<KeyT, ValueT>, uint32_t> PairT;
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
            } while (_bucket < _map->_num_buckets && (_map->MAPA_BUCKET(_pairs, _bucket)) == INACTIVE);
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
            } while (_bucket < _map->_num_buckets && (_map->MAPA_BUCKET(_pairs, _bucket)) == INACTIVE);
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
        max_load_factor(0.9f);
    }

    HashMap(uint32_t bucket = 8)
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

        if (std::is_integral<KeyT>::value && std::is_trivially_copyable<ValueT>::value) {
            memcpy(_pairs, other._pairs, other._num_buckets * sizeof(PairT));
        }
        else {
            auto old_pairs = other._pairs;
            for (uint32_t bucket = 0; bucket < _num_buckets; bucket++) {
                const auto next_bucket = NEXT_BUCKET(old_pairs, bucket);
                if (next_bucket != INACTIVE) {
                    new(_pairs + bucket) PairT(old_pairs[bucket]);
                }
                _pairs[bucket]._bucket = old_pairs[bucket]._bucket;
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

        clear();
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
        clear();
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
        return  {this, _num_buckets};
    }

    const_iterator cend() const
    {
        return {this, _num_buckets};
    }

    const_iterator end()
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
    //bucket < _num_buckets, count items in bucket,
    //if bucket is not the main bucket, search from main bucket
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
            printf("  %2d  %8d  %.2lf  %.2lf\n", i, bucketsi, bucketsi * 100.0 * i / _num_filled, sumn * 100.0 / _num_filled);
        }

        puts("========== collision miss ration ===========");
        for (uint32_t i = 0; i < sizeof(steps) / sizeof(steps[0]); i++) {
            sumc += steps[i];
            if (steps[i] <= 2)
                continue;
            printf("  %2d  %8d  %.2lf  %.2lf\n", i, steps[i], steps[i] * 100.0 / collision, sumc * 100.0 / collision);
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
        const auto bucket = find_filled_bucket(key);
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
        const auto find = NEXT_BUCKET(_pairs, bucket) == INACTIVE;
        if (find) {
            if (EMILIB_UNLIKELY(check_expand_need()))
                bucket = find_unique_bucket(key);

            NEW_KVALUE(key, value, bucket); _num_filled++;
            SET_BIT(bucket);
        }
        return { iterator(this, bucket), find };
    }

    std::pair<iterator, bool> insert(KeyT&& key, ValueT&& value)
    {
        auto bucket = find_or_allocate(key);
        const auto find = NEXT_BUCKET(_pairs, bucket) == INACTIVE;
        if (find) {
            if (check_expand_need())
                bucket = find_unique_bucket(key);

            NEW_KVALUE(std::move(key), std::move(value), bucket); _num_filled++;
            SET_BIT(bucket);
        }
        return { iterator(this, bucket), find };
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

#if 0
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
        if (NEXT_BUCKET(_pairs, bucket) != INACTIVE) {
            GET_VAL(_pairs, bucket) = value;
        }
        else {
            NEW_KVALUE(key, value, bucket); _num_filled++;
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
            NEW_KVALUE(key, new_value, bucket); _num_filled++;
            return ValueT();
        }
    }
#endif

    /// Like std::map<KeyT,ValueT>::operator[].
    ValueT& operator[](const KeyT& key)
    {
        auto bucket = find_or_allocate(key);
        /* Check if inserting a new value rather than overwriting an old entry */
        if (NEXT_BUCKET(_pairs, bucket) == INACTIVE) {
            if (EMILIB_UNLIKELY(check_expand_need()))
                bucket = find_unique_bucket(key);

            NEW_KVALUE(key, ValueT(), bucket); _num_filled++;
            SET_BIT(bucket);
        }

        return GET_VAL(_pairs, bucket);
    }

    // -------------------------------------------------------

    /// Erase an element from the hash table.
    /// return false if element was not found
    size_t erase(const KeyT& key)
    {
        const auto bucket = erase_from_key(key);
        if (bucket == INACTIVE)
            return 0;

        CLS_BUCKET(_pairs, bucket); _pairs[bucket].~PairT(); _num_filled -= 1;
        CLR_BIT(bucket);
//        assert(_pairs[bucket]._bucket == EMILIB_BUCKET_NONE);

#ifdef EMILIB_AUTO_SHRINK
        //        if (_num_buckets > 254 && _num_buckets > 4 * _num_filled)
        //            rehash(_num_filled / max_load_factor()  + 2);
#endif
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
            return ++it;
        }
#endif
        //assert(it->first == GET_KEY(_pairs, it._bucket));
        const auto bucket = erase_from_bucket(it._bucket);
        //assert(bucket != INACTIVE);
        CLS_BUCKET(_pairs, bucket); _pairs[bucket].~PairT(); _num_filled -= 1;
        CLR_BIT(bucket);
        //assert(_pairs[bucket]._bucket == EMILIB_BUCKET_NONE);

        //erase from main bucket, return main bucket as next
        if (bucket != it._bucket)
            return {this, it._bucket};

        iterator itnext = {this, it._bucket};
        return ++itnext;
    }

    /// Remove all elements, keeping full capacity.
    void clear()
    {
        if (_num_filled > _num_buckets / 4 && std::is_integral<KeyT>::value && std::is_trivially_copyable<ValueT>::value) {
            _num_filled = 0;
            memset(_pairs, EMILIB_BUCKET_NONE, sizeof(_pairs[0]) * _num_buckets);
            return;
        }

        for (uint32_t bucket = 0; bucket < _num_buckets ; ++bucket) {
            if (NEXT_BUCKET(_pairs, bucket) != INACTIVE) {
                _pairs[bucket].~PairT(); _num_filled -= 1;
            }
            RST_BUCKET(_pairs, bucket);
        }
    }

    /// Make room for this many elements
    bool reserve(uint32_t num_elems)
    {
        //auto required_buckets = (uint32_t)(((size_t)num_elems * _loadlf) >> 20) + 2;
        const auto required_buckets = num_elems * 10 / 8 + 2;
        if (EMILIB_LIKELY(required_buckets <= _num_buckets))
            return false;

        rehash(required_buckets);
        return true;
    }

    /// Make room for this many elements
    void rehash(uint32_t required_buckets)
    {
        uint32_t num_buckets = 8;
        if (required_buckets >= 1024)
            num_buckets = 1024 * 2;
        while (num_buckets < required_buckets) { num_buckets *= 2; }

        assert(num_buckets <= 1 + (EMILIB_BUCKET_MASK >> EMILIB_HASH_BIT));

        auto new_pairs = (PairT*)malloc(num_buckets * sizeof(PairT));
        /*
           if (!new_pairs) {
           throw std::bad_alloc();
           }
         **/
        auto old_num_filled  = _num_filled;
        auto old_num_buckets = _num_buckets;
        auto old_pairs = _pairs;

        _num_filled  = 0;
        _num_buckets = num_buckets;
        _mask        = num_buckets - 1;
        _pairs       = new_pairs;

        if (sizeof(PairT) <= sizeof(int64_t) * 4)
            memset(_pairs, EMILIB_BUCKET_NONE, sizeof(_pairs[0]) * num_buckets);
        else
            for (uint32_t bucket = 0; bucket < num_buckets; bucket++)
                RST_BUCKET(_pairs, bucket);

        uint32_t collision = 0;
        //set all main bucket first
        for (uint32_t src_bucket = 0; src_bucket < old_num_buckets; src_bucket++) {
            if (NEXT_BUCKET(old_pairs, src_bucket) == INACTIVE)
                continue;

            const auto hashkey = (uint32_t)_hasher(GET_KEY(old_pairs, src_bucket));
            const auto main_bucket = hashkey & _mask;

            auto next_bucket = NEXT_BUCKET(_pairs, main_bucket);
            auto& old_pair = old_pairs[src_bucket];
            if (next_bucket == INACTIVE) {
                new(_pairs + main_bucket) PairT(std::move(old_pair)); old_pair.~PairT();
                SET_BUCKET(main_bucket, main_bucket);
                SET_BIT(main_bucket);
            }
            else {
                ADDR_BUCKET(old_pairs, collision++) = src_bucket;
            }

            SET_MHASH(main_bucket, hashkey);
            _num_filled += 1;
            if (_num_filled >= old_num_filled)
                break;
        }

        //reset all collisions bucket
        for (uint32_t colls = 0; colls < collision; colls++) {
            const auto src_bucket = ADDR_BUCKET(old_pairs, colls);

            const auto hashkey = (uint32_t)_hasher(GET_KEY(old_pairs, src_bucket));
            const auto main_bucket = hashkey & _mask;
            //
            auto& old_pair = old_pairs[src_bucket];
            {
                auto next_bucket = NEXT_BUCKET(_pairs, main_bucket);
                //assert(next_bucket != INACTIVE);
                //check current bucket_key is in main bucket or not
#if 1
                if (next_bucket != main_bucket)
                    next_bucket = find_last_bucket(next_bucket);

                //find a new empty and link it to tail
                auto new_bucket = find_empty_bucket(next_bucket);
                //assert(_pairs[new_bucket]._bucket == EMILIB_BUCKET_NONE);

                new(_pairs + new_bucket) PairT(std::move(old_pair)); old_pair.~PairT();
                SET_BUCKET(new_bucket, new_bucket);
                SET_BUCKET(next_bucket, new_bucket);
                SET_BIT(new_bucket);
#else
                auto new_bucket = find_empty_bucket(next_bucket);
                SET_BUCKET(main_bucket, new_bucket);
                new(_pairs + new_bucket) PairT(std::move(old_pair)); old_pair.~PairT();
                SET_BUCKET(new_bucket, next_bucket == main_bucket ? new_bucket : next_bucket);
#endif
            }
        }

#ifdef EMILIB_REHASH_LOG
        if (_num_filled > 0) {
            char buff[255] = {0};
            sprintf(buff, "    _num_filled/K.V/pack/collision = %u/%s.%s/%zd/%.2lf%%", _num_filled, typeid(KeyT).name(), typeid(ValueT).name(), sizeof(_pairs[0]), (collision * 100.0 / _num_filled));
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
        const auto hashkey = (uint32_t)_hasher(key);
        const auto bucket = hashkey & _mask;

        auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        if (next_bucket == INACTIVE)
            return INACTIVE;

        const auto bqKey = key == GET_KEY(_pairs, bucket);
        if (bqKey) {
            CLEAR_MHASH(bucket);
            if (next_bucket == bucket)
                return bucket;

            const auto nbucket = NEXT_BUCKET(_pairs, next_bucket);
            GET_PVAL(_pairs, bucket).swap(GET_PVAL(_pairs, next_bucket));
            SET_BUCKET(bucket, (nbucket == next_bucket) ? bucket : nbucket);
            return next_bucket;
        }
        else if (next_bucket == bucket || TST_HASH(bucket, hashkey))
            return INACTIVE;

        auto prev_bucket = bucket;
        while (true) {
            const auto nbucket = NEXT_BUCKET(_pairs, next_bucket);
            if (key == GET_KEY(_pairs, next_bucket)) {
                SET_BUCKET(prev_bucket, (nbucket == next_bucket) ? prev_bucket : nbucket);
                CLEAR_MHASH(bucket);
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
        //assert(next_bucket != INACTIVE && next_bucket < _num_buckets);
        const auto main_bucket = hash_bucket(GET_KEY(_pairs, bucket));
        CLEAR_MHASH(main_bucket);

        if (bucket == main_bucket) {
            //more than one bucket
            if (bucket != next_bucket) {
                const auto nbucket = NEXT_BUCKET(_pairs, next_bucket);
                GET_PVAL(_pairs, bucket).swap(GET_PVAL(_pairs, next_bucket));
                SET_BUCKET(bucket, (nbucket == next_bucket) ? bucket : nbucket);
            }

            return next_bucket;
        }

        const auto prev_bucket = find_prev_bucket(main_bucket, bucket);
        SET_BUCKET(prev_bucket, (bucket == next_bucket) ? prev_bucket : next_bucket);
        return bucket;
    }

    // Find the bucket with this key, or return INACTIVE
    uint32_t find_filled_bucket(const KeyT& key) const
    {
        const auto hashkey = (uint32_t)_hasher(key);
        const auto bucket  = hashkey & _mask;

        auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        if (next_bucket == INACTIVE)
            return INACTIVE;
        else if (key == GET_KEY(_pairs, bucket))
            return bucket;
        else if (next_bucket == bucket || TST_HASH(bucket, hashkey))
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
        //assert(_pairs[new_bucket]._bucket == EMILIB_BUCKET_NONE);

        const auto prev_bucket = find_prev_bucket(main_bucket, bucket);
        SET_BUCKET(prev_bucket, new_bucket);

        new(_pairs + new_bucket) PairT(std::move(_pairs[bucket])); _pairs[bucket].~PairT();
        SET_BUCKET(new_bucket, (next_bucket == bucket) ? new_bucket : next_bucket);

        CLS_BUCKET(_pairs, bucket);
        SET_BIT(new_bucket);
        return new_bucket;
    }

    // Find the bucket with this key, or return a good empty bucket to place the key in.
    // In the latter case, the bucket is expected to be filled.
    // If the bucket opt by other Key who's main bucket is not in this postion, kick it out
    // and move it to a new empty postion.
    uint32_t find_or_allocate(const KeyT& key)
    {
        const auto hashkey = (uint32_t)_hasher(key);
        const auto bucket = hashkey & _mask;

        const auto& bucket_key = GET_KEY(_pairs, bucket);
        auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        //no need set main bucket
        if (next_bucket == INACTIVE || key == bucket_key) {
            SET_MHASH(bucket, hashkey);
            return bucket;
        }

        //check current bucket_key is in main bucket or not
        const auto main_bucket = hash_bucket(bucket_key);
        if (main_bucket != bucket) {
            reset_main_bucket(main_bucket, bucket);
            SET_MHASH(bucket, hashkey);
            return bucket;
        }
        else if (next_bucket == bucket) {
            const auto new_bucket = find_empty_bucket(next_bucket);
            SET_BUCKET(next_bucket, new_bucket);
            SET_MHASH(bucket, hashkey);
            return new_bucket;
        }
        else if (TST_HASH(bucket, hashkey)) {
            const auto last_bucket = find_last_bucket(next_bucket);
            const auto new_bucket = find_empty_bucket(last_bucket);
            SET_BUCKET(last_bucket, new_bucket);
            SET_MHASH(bucket, hashkey);
            return new_bucket;
        }

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

        SET_MHASH(bucket, hashkey);
        const auto new_bucket = find_empty_bucket(next_bucket);
        SET_BUCKET(next_bucket, new_bucket);
        return new_bucket;
    }

    // key is not in this map. Find a empty place to put it.
    // combine linear probing and quadratic probing
    uint32_t find_empty_bucket(uint32_t bucket_from)
    {
        const auto bucket = ++ bucket_from & _mask;
/**
        if (NEXT_BUCKET(_pairs, bucket) == INACTIVE) {
            return bucket;
        }
**/
        const auto bofset = bucket % EMILIB_HASH_BIT;
        auto mask_bucket = bucket - bofset;
        auto bmask = GET_BIT(mask_bucket) & ~((1 << bofset) - 1);
        if (bmask != 0) {
            return mask_bucket + CTZ(bmask);
        }

        mask_bucket = (mask_bucket + EMILIB_HASH_BIT) & _mask;
        bmask = GET_BIT(mask_bucket);
        if (bmask > 0) {
            return mask_bucket + CTZ(bmask);
        }

        bucket_from = (mask_bucket + EMILIB_HASH_BIT) & _mask;
         for (uint32_t slot = 1; ; ++slot) {
            const auto empty_bucket = TST_BIT(bucket_from);
            if (empty_bucket != INACTIVE)
                return empty_bucket;

            bucket_from += EMILIB_HASH_BIT;
            if (slot > 4)
                bucket_from += _num_buckets / 2;
            bucket_from = (bucket_from + slot * EMILIB_HASH_BIT) & _mask;
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
        const auto hashkey = (uint32_t)_hasher(key);
        const auto bucket = hashkey & _mask;

        auto next_bucket = NEXT_BUCKET(_pairs, bucket);
        SET_MHASH(bucket, hashkey);
        if (next_bucket == INACTIVE) {
            SET_BIT(bucket);
            return bucket;
        }

        //check current bucket_key is in main bucket or not
        const auto main_bucket = hash_bucket(GET_KEY(_pairs, bucket));
        if (EMILIB_UNLIKELY(main_bucket != bucket)) {
            reset_main_bucket(main_bucket, bucket);
            SET_MHASH(bucket, hashkey);
            SET_BIT(bucket);
            return bucket;
        }
        else if (next_bucket != bucket)
            next_bucket = find_last_bucket(next_bucket);

        const auto new_bucket = find_empty_bucket(next_bucket);
        SET_BUCKET(next_bucket, new_bucket);
        SET_BIT(new_bucket);
        return new_bucket;
    }

    //Thomas Wang's function
    //https://gist.github.com/badboy/6267743
    static inline uint32_t hash32(uint32_t key)
    {
        uint64_t const r = key * UINT64_C(0xca4bcaa75ec3f625);
        uint32_t h = static_cast<uint32_t>(r >> 32);
        uint32_t l = static_cast<uint32_t>(r);
        return h + l;
#if 0
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
        key = (~key) + (key << 18); // key = (key << 18) - key - 1;
        key = key ^ (key >> 31);
        key = key * 21; // key = (key + (key << 2)) + (key << 4);
        key = key ^ (key >> 11);
        key = key + (key << 6);
        key = key ^ (key >> 22);
        return (uint32_t) key;
    }

    template<typename UType, typename std::enable_if<std::is_integral<UType>::value, long>::type = 0>
        //    template<class UType, class = typename std::enable_if<std::is_integral<UType>::value, int>::type>
    inline uint32_t hash_bucket(const UType key) const
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
    inline uint32_t hash_bucket(const UType& key) const
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

    uint32_t  _loadlf;
};

} // namespace emilib

#if __cplusplus > 199711
//template <class Key, class Val> using emihash = emilib1::HashMap<Key, Val, std::hash<Key>>;
#endif
