// By Emil Ernerfeldt 2014-2017
// LICENSE:
//   This software is dual-licensed to the public domain and under the following
//   license: you are granted a perpetual, irrevocable license to copy, modify,
//   publish, and distribute this file as you see fit.
//   http://www.ilikebigbits.com/2016_08_28_hash_table.html

//some others
//https://tessil.github.io/2016/08/29/benchmark-hopscotch-map.html

#pragma once

#include <cstdlib>
#include <iterator>
#include <utility>

//#include <loguru.hpp>
#if 1
    #define HASH_FUNC     _hasher
#else
    #define HASH_FUNC(v)   (v & 0xFFFFFFFF) //for integer test.
#endif

#define GET_KEY(p,n) p[n].second.first
#define GET_VAL(p,n) p[n].second.second
#define GET_STA(s,n) s[n].first


namespace emilib2 {
enum State
{
    INACTIVE = -1, // Never been touched
    FILLED = 0   // Is set with key/value
};

/// like std::equal_to but no need to #include <functional>
template<typename T>
struct HashMapEqualTo
{
    constexpr bool operator()(const T& lhs, const T& rhs) const
    {
        return lhs == rhs;
    }
};

/// A cache-friendly hash table with open addressing, linear probing and power-of-two capacity
template <typename KeyT, typename ValueT, typename HashT = std::hash<KeyT>, typename EqT = HashMapEqualTo<KeyT>>
class HashMap
{
private:
    typedef  HashMap<KeyT, ValueT, HashT, EqT> MyType;

    typedef  std::pair<int, std::pair<KeyT, ValueT>> PairT;

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
        typedef std::pair<KeyT, ValueT>   value_type;
        typedef value_type*               pointer;
        typedef value_type&               reference;

        iterator() { }

        iterator(MyType* hash_map, size_t bucket) : _map(hash_map), _bucket(bucket)
        {
        }

        iterator& operator++()
        {
            this->goto_next_element();
            return *this;
        }

        iterator operator++(int)
        {
            size_t old_index = _bucket;
            this->goto_next_element();
            return iterator(_map, old_index);
        }

        reference operator*() const
        {
            return _map->_pairs[_bucket].second;
        }

        pointer operator->() const
        {
            return &(_map->_pairs[_bucket].second);
        }

        bool operator==(const iterator& rhs) const
        {
            //DCHECK_EQ_F(_map, rhs._map);
            return this->_bucket == rhs._bucket;
        }

        bool operator!=(const iterator& rhs) const
        {
            //DCHECK_EQ_F(_map, rhs._map);
            return this->_bucket != rhs._bucket;
        }

    private:
        void goto_next_element()
        {
            //DCHECK_LT_F(_bucket, _map->_num_buckets);
            do {
                _bucket++;
            } while (_bucket < _map->_num_buckets && _map->GET_STA(_pairs,_bucket) == State::INACTIVE);
        }

    //private:
    //    friend class MyType;
    public:
        MyType* _map;
        size_t  _bucket;
    };

    class const_iterator
    {
    public:
        typedef std::forward_iterator_tag iterator_category;
        typedef size_t                    difference_type;
        typedef size_t                    distance_type;
        typedef std::pair<KeyT, ValueT>   value_type;
        typedef value_type*               pointer;
        typedef value_type&               reference;

        const_iterator() { }

        const_iterator(iterator proto) : _map(proto._map), _bucket(proto._bucket)
        {
        }

        const_iterator(const MyType* hash_map, size_t bucket) : _map(hash_map), _bucket(bucket)
        {
        }

        const_iterator& operator++()
        {
            this->goto_next_element();
            return *this;
        }

        const_iterator operator++(int)
        {
            size_t old_index = _bucket;
            this->goto_next_element();
            return const_iterator(_map, old_index);
        }

        reference operator*() const
        {
            return _map->_pairs[_bucket];
        }

        pointer operator->() const
        {
            return _map->_pairs + _bucket;
        }

        bool operator==(const const_iterator& rhs) const
        {
            //DCHECK_EQ_F(_map, rhs._map);
            return this->_bucket == rhs._bucket;
        }

        bool operator!=(const const_iterator& rhs) const
        {
            //DCHECK_EQ_F(_map, rhs._map);
            return this->_bucket != rhs._bucket;
        }

    private:
        void goto_next_element()
        {
            //DCHECK_LT_F(_bucket, _map->_num_buckets);
            do {
                _bucket++;
            } while (_bucket < _map->_num_buckets && _map->GET_STA(_pairs,_bucket) == State::INACTIVE);
        }

    //private:
    //    friend class MyType;
    public:
        const MyType* _map;
        size_t        _bucket;
    };

    // ------------------------------------------------------------------------

    HashMap()
    {
       _num_buckets      = 0;
       _num_filled       = 0;
       _mask             = 0;  // _num_buckets minus one
       _pairs            = nullptr;
       reserve(8);
    }

    HashMap(const HashMap& other)
    {
        reserve(other.size());
        insert(other.cbegin(), other.cend());
    }

    HashMap(HashMap&& other)
    {
        *this = std::move(other);
    }

    HashMap& operator=(const HashMap& other)
    {
        clear();
        reserve(other.size());
        insert(other.cbegin(), other.cend());
        return *this;
    }

    HashMap& operator=(HashMap&& other)
    {
        this->swap(other);
        return *this;
    }

    ~HashMap()
    {
        for (size_t bucket=0; bucket<_num_buckets; ++bucket) {
            if (GET_STA(_pairs,bucket) != State::INACTIVE) {
                _pairs[bucket].~PairT();
            }
        }
        if (_pairs)
            free(_pairs);
    }

    void swap(HashMap& other)
    {
        std::swap(_hasher,           other._hasher);
        std::swap(_eq,               other._eq);
        std::swap(_pairs,            other._pairs);
        std::swap(_num_buckets,      other._num_buckets);
        std::swap(_num_filled,       other._num_filled);
        std::swap(_mask,             other._mask);
    }

    // -------------------------------------------------------------

    iterator begin()
    {
        size_t bucket = 0;
        while (bucket<_num_buckets && GET_STA(_pairs,bucket) == State::INACTIVE) {
            ++bucket;
        }
        return iterator(this, bucket);
    }

    const_iterator cbegin() const
    {
        size_t bucket = 0;
        while (bucket<_num_buckets && GET_STA(_pairs,bucket) == State::INACTIVE) {
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

    size_t size() const
    {
        return _num_filled;
    }

    bool empty() const
    {
        return _num_filled==0;
    }

    // Returns the number of buckets.
    size_t bucket_count() const
    {
        return _num_buckets;
    }

    /// Returns average number of elements per bucket.
    float load_factor() const
    {
        return static_cast<float>(_num_filled) / static_cast<float>(_num_buckets);
    }

    // ------------------------------------------------------------

    template<typename KeyLike>
    iterator find(const KeyLike& key)
    {
        auto bucket = this->find_filled_bucket(key);
        if (bucket == (size_t)-1) {
            return this->end();
        }
        return iterator(this, bucket);
    }

    template<typename KeyLike>
    const_iterator find(const KeyLike& key) const
    {
        auto bucket = this->find_filled_bucket(key);
        if (bucket == (size_t)-1)
        {
            return this->end();
        }
        return const_iterator(this, bucket);
    }

    template<typename KeyLike>
    bool contains(const KeyLike& k) const
    {
        return find_filled_bucket(k) != (size_t)-1;
    }

    template<typename KeyLike>
    size_t count(const KeyLike& k) const
    {
        return find_filled_bucket(k) != (size_t)-1 ? 1 : 0;
    }

    /// Returns the matching ValueT or nullptr if k isn't found.
    template<typename KeyLike>
    ValueT* try_get(const KeyLike& k)
    {
        auto bucket = find_filled_bucket(k);
        if (bucket != (size_t)-1) {
            return &GET_VAL(_pairs,bucket);
        } else {
            return nullptr;
        }
    }

    /// Const version of the above
    template<typename KeyLike>
    const ValueT* try_get(const KeyLike& k) const
    {
        auto bucket = find_filled_bucket(k);
        if (bucket != (size_t)-1) {
            return &GET_VAL(_pairs,bucket);
        } else {
            return nullptr;
        }
    }

    /// Convenience function.
    template<typename KeyLike>
    const ValueT get_or_return_default(const KeyLike& k) const
    {
        const ValueT* ret = try_get(k);
        if (ret) {
            return *ret;
        } else {
            return ValueT();
        }
    }

    // -----------------------------------------------------

    /// Returns a pair consisting of an iterator to the inserted element
    /// (or to the element that prevented the insertion)
    /// and a bool denoting whether the insertion took place.
    std::pair<iterator, bool> insert(const KeyT& key, const ValueT& value)
    {
        check_expand_need();

        auto bucket = find_or_allocate(key);

        if (GET_STA(_pairs,bucket) != State::INACTIVE) {
            return { iterator(this, bucket), false };
        } else {
            new(_pairs + bucket) PairT(State::FILLED, std::pair<KeyT, ValueT>(key, value));
            _num_filled++;
            return { iterator(this, bucket), true };
        }
    }

    std::pair<iterator, bool> insert(const std::pair<KeyT, ValueT>& p)
    {
        return insert(p.first, p.second);
    }

    void insert(const_iterator begin, const_iterator end)
    {
        // TODO: reserve space exactly once.
        for (; begin != end; ++begin) {
            insert(begin->first, begin->second);
        }
    }

    /// Same as above, but contains(key) MUST be false
    void insert_unique(KeyT&& key, ValueT&& value)
    {
        //DCHECK_F(!contains(key));
        check_expand_need();
        auto bucket = find_main_bucket(key);
        new(_pairs + bucket) PairT(State::FILLED, std::pair<KeyT, ValueT>(std::move(key), std::move(value)));
        _num_filled++;
    }

    void insert_unique(std::pair<KeyT, ValueT>&& p)
    {
        insert_unique(std::move(p.first), std::move(p.second));
    }

    void insert_or_assign(const KeyT& key, ValueT&& value)
    {
        check_expand_need();

        auto bucket = find_or_allocate(key);

        // Check if inserting a new value rather than overwriting an old entry
        if (GET_STA(_pairs,bucket) != State::INACTIVE) {
            GET_VAL(_pairs,bucket) = value;
        } else {
            new(_pairs + bucket) PairT(State::FILLED, std::pair<KeyT, ValueT>(key, value));
            _num_filled++;
        }
    }

    /// Return the old value or ValueT() if it didn't exist.
    ValueT set_get(const KeyT& key, const ValueT& new_value)
    {
        check_expand_need();

        auto bucket = find_or_allocate(key);

        // Check if inserting a new value rather than overwriting an old entry
        if (GET_STA(_pairs,bucket) != State::INACTIVE) {
            ValueT old_value = GET_VAL(_pairs,bucket);
            GET_VAL(_pairs,bucket) = new_value.second;
            return old_value;
        } else {
            new(_pairs + bucket) PairT(State::FILLED, std::pair<KeyT, ValueT>(key, new_value));
            _num_filled++;
            return ValueT();
        }
    }

    /// Like std::map<KeyT,ValueT>::operator[].
    ValueT& operator[](const KeyT& key)
    {
        check_expand_need();

        auto bucket = find_or_allocate(key);

        /* Check if inserting a new value rather than overwriting an old entry */
        if (GET_STA(_pairs,bucket) == State::INACTIVE) {
            new(_pairs + bucket) PairT(State::FILLED, std::pair<KeyT, ValueT>(key, ValueT()));
            _num_filled++;
        }

        return GET_VAL(_pairs,bucket);
    }

    // -------------------------------------------------------

    /// Erase an element from the hash table.
    /// return false if element was not found
    bool erase(const KeyT& key)
    {
        auto bucket = erase_from_bucket (key);
        if (bucket != (size_t)-1) {
            GET_STA(_pairs,bucket) = State::INACTIVE;
            _pairs[bucket].~PairT();
            _num_filled -= 1;
            return true;
        } else {
            return false;
        }
    }

    /// Erase an element typedef an iterator.
    /// Returns an iterator to the next element (or end()).
    iterator erase(iterator it)
    {
        //DCHECK_EQ_F(it._map, this);
        //DCHECK_LT_F(it._bucket, _num_buckets);
        auto bucket = it._bucket;
        if (GET_STA(_pairs,bucket) > State::FILLED) {
             bucket = erase_from_bucket(it->first);
        }

        GET_STA(_pairs,bucket) = State::INACTIVE;
        _pairs[bucket].~PairT();
        _num_filled -= 1;
        if (bucket == it._bucket)
           ++it;
        return it;
    }

    /// Remove all elements, keeping full capacity.
    void clear()
    {
        for (size_t bucket=0; bucket<_num_buckets; ++bucket) {
            if (GET_STA(_pairs,bucket) != State::INACTIVE) {
                GET_STA(_pairs,bucket) = State::INACTIVE;
                _pairs[bucket].~PairT();
            }
        }
        _num_filled = 0;
    }

    /// Make room for this many elements
    void reserve(size_t num_elems)
    {
        size_t required_buckets = num_elems + 1;
        if (required_buckets <= _num_buckets) {
            return;
        }
        size_t num_buckets = 4;
        while (num_buckets < required_buckets) { num_buckets *= 2; }

        auto new_pairs  = (PairT*)malloc(num_buckets * sizeof(PairT));

        if (!new_pairs) {
            free(new_pairs);
            throw std::bad_alloc();
        }

        //auto old_num_filled  = _num_filled;
        auto old_num_buckets = _num_buckets;
        auto old_pairs       = _pairs;

        _num_filled  = 0;
        _num_buckets = num_buckets;
        _mask        = _num_buckets - 1;
        _pairs       = new_pairs;

        for (int i = 0; i < num_buckets; i++)
            GET_STA(new_pairs, i) = State::INACTIVE;

        for (size_t src_bucket=0; src_bucket<old_num_buckets; src_bucket++) {
            if (GET_STA(old_pairs, src_bucket) == State::INACTIVE) {
                continue;
            }

            auto& src_pair = old_pairs[src_bucket];
            auto dst_bucket = find_main_bucket(src_pair.second.first);
             new(_pairs + dst_bucket) PairT(std::move(src_pair));
             GET_STA(_pairs,dst_bucket) = State::FILLED;
             _num_filled += 1;
             src_pair.~PairT();
        }

        //DCHECK_EQ_F(old_num_filled, _num_filled);

        free(old_pairs);
    }

private:
    // Can we fit another element?
    void check_expand_need()
    {
        reserve(_num_filled + 1);
    }

    size_t erase_from_bucket(const KeyT& key) const
    {
        //if (empty()) { return (size_t)-1; } // Optimization
        const auto bucket = HASH_FUNC(key) & _mask;
        auto& max_probe_length = GET_STA(_pairs,bucket);
        if (max_probe_length == State::INACTIVE)
            return (size_t)-1;
        else if (_eq(GET_KEY(_pairs,bucket), key)) {
            if (max_probe_length == State::FILLED)
                return bucket;

            //find next bucket and swap position, can not erase the main bucket if max_probe_length > 1
            for (int offset = max_probe_length; offset > 0; offset--) {
                max_probe_length = offset - 1;
                const auto nbucket = (bucket + offset) & _mask;
                if (GET_STA(_pairs,nbucket) != State::INACTIVE && bucket == (HASH_FUNC(GET_KEY(_pairs,nbucket)) & _mask)) {
                    std::swap(_pairs[nbucket].second, _pairs[bucket].second);
                    return nbucket;
                }
            }
            return bucket;
        }
        else if (max_probe_length == State::FILLED)
            return (size_t)-1;

        for (auto offset = max_probe_length; offset > 0; offset--) {
            const auto nbucket = (bucket + offset) & _mask;
            if (GET_STA(_pairs,nbucket) != State::INACTIVE && _eq(GET_KEY(_pairs,nbucket), key)) {
                if (offset == max_probe_length)
                    max_probe_length = offset - 1;
                return nbucket;
            }
        }

        return (size_t)-1;
    }

    // Find the bucket with this key, or return (size_t)-1
    size_t find_filled_bucket(const KeyT& key) const
    {
        //if (empty()) { return (size_t)-1; } // Optimization
        const auto bucket = HASH_FUNC(key) & _mask;
        const auto& max_probe_length = GET_STA(_pairs,bucket);
#if 0
        if (max_probe_length == State::FILLED) {
            if (_eq(GET_KEY(_pairs,bucket), key))
                return bucket;
            else
                return (size_t)-1;
        }
        else if (max_probe_length == State::INACTIVE)
            return (size_t)-1;
#elif 1
        if (max_probe_length == State::INACTIVE)
            return (size_t)-1;
        else if (_eq(GET_KEY(_pairs,bucket), key))
            return bucket;
        else if (max_probe_length == State::FILLED)
            return (size_t)-1;
#endif

        for (int offset = max_probe_length; offset > 0; offset--) {
            const auto nbucket = (bucket + offset) & _mask;
            if (GET_STA(_pairs,nbucket) != State::INACTIVE && _eq(GET_KEY(_pairs,nbucket), key)) {
//                _pairs[nbucket].second.swap(_pairs[bucket].second);
                return bucket;
            }
        }

        return (size_t)-1;
    }

    inline void reset_main_bucket(size_t main_bucket, size_t bucket)
    {
        auto& max_probe_length = GET_STA(_pairs, main_bucket);
        if (max_probe_length != State::INACTIVE) {
            const auto new_bucket = find_empty_bucket(main_bucket + 1);
            new(_pairs + new_bucket) PairT(State::FILLED, std::pair<KeyT, ValueT>(GET_KEY(_pairs, bucket), GET_VAL(_pairs, bucket)));
            if (new_bucket > main_bucket + max_probe_length)
                max_probe_length = int(new_bucket - main_bucket);
            else if (new_bucket < main_bucket)
                max_probe_length = int(_num_buckets + new_bucket - main_bucket);
        }
        else {
            new(_pairs + main_bucket) PairT(State::FILLED, std::pair<KeyT, ValueT>(GET_KEY(_pairs, bucket), GET_VAL(_pairs, bucket)));
        }
    }

    // Find the bucket with this key, or return a good empty bucket to place the key in.
    // In the latter case, the bucket is expected to be filled.
    size_t find_or_allocate(const KeyT& key)
    {
        const auto bucket = HASH_FUNC(key) & _mask;
        auto& max_probe_length = GET_STA(_pairs,bucket);
        const auto& bucket_key = GET_KEY(_pairs,bucket);
        if (max_probe_length == State::INACTIVE)
            return bucket;
        else if (_eq(bucket_key, key))
            return bucket;

        //find main postion
        const auto main_bucket = HASH_FUNC(bucket_key) & _mask;
        if (main_bucket != bucket) {
            reset_main_bucket(main_bucket, bucket);
            GET_STA(_pairs,bucket) = State::INACTIVE;
            return bucket;
        }

        //find exits
        size_t hole = (size_t)-1;
        for (int offset = max_probe_length; offset > 0; offset--) {
            const auto nbucket = (bucket + offset) & _mask;
            if (GET_STA(_pairs,nbucket) != State::INACTIVE) {
                if (_eq(GET_KEY(_pairs,nbucket), key))
                    return nbucket;
            } else if (hole == (size_t)-1)
                hole = nbucket;
        }

        if (hole != (size_t)-1)
            return hole;

        //find a new empty
        for (int offset = max_probe_length + 1; ; ++offset) {
            const auto nbucket = (bucket + offset) & _mask;
            if (GET_STA(_pairs,nbucket) == State::INACTIVE) {
                max_probe_length = offset;
                return nbucket;
            }
        }
    }

    size_t find_main_bucket(const KeyT& key)
    {
        const auto bucket = HASH_FUNC(key) & _mask;
        auto& max_probe_length = GET_STA(_pairs,bucket);
        if (max_probe_length == State::INACTIVE)
            return bucket;

        //find main postion
        const auto main_bucket = HASH_FUNC(GET_KEY(_pairs,bucket)) & _mask;
        if (main_bucket != bucket) {
            reset_main_bucket(main_bucket, bucket);
            return bucket;
        }

        for (int offset = 1; ; ++offset) {
            const auto nbucket = (bucket + offset) & _mask;
            if (GET_STA(_pairs,nbucket) == State::INACTIVE) {
                if (offset > max_probe_length)
                    max_probe_length = offset;
                return nbucket;
            }
        }
    }

    // key is not in this map. Find a place to put it.
    inline size_t find_empty_bucket(size_t bucket_from)
    {
        for (int offset = 0; ; ++offset) {
            const auto bucket = (bucket_from + offset) & _mask;
            if (GET_STA(_pairs,bucket) == State::INACTIVE)
                return bucket;
        }
    }

private:

    HashT   _hasher;
    EqT     _eq;
    PairT*  _pairs            ;
    size_t  _num_buckets      ;
    size_t  _num_filled       ;
    size_t  _mask             ;  // _num_buckets minus one
};

} // namespace emilib
