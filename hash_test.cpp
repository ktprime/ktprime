
#include <random>
#include <map>
#include <ctime>
#include <cassert>
#include <iostream>
#include <string>
#include <algorithm>
#include <array>
#include "sfc64.h"

//#define TP                1
#define HOOD_HASH           1
//#define EMILIB_LRU_SET      1
//#define EMILIB_IDENTITY_HASH 1
//#define EMILIB_REHASH_LOG   1
//#define EMILIB_SAFE_HASH      1
//#define EMILIB_AVX_MEMCPY 1
//#define EMILIB_STATIS    1
//#define EMILIB_HIGH_LOAD 1001000

#ifndef TKey
    #define TKey  1
#endif
#ifndef TVal
    #define TVal  1
#endif

//#define EMILIB_BUCKET_INDEX 0
//#define EMILIB_BUCKET_INDEX  2

#include "hash_table55.hpp"
#include "hash_table57.hpp"
#include "hash_set6.hpp" //

#include "hash_table52.hpp"
#include "hash_table53.hpp"
//#include "hash_table54.hpp"
#include "hash_table64.hpp"

////some others
//https://github.com/ilyapopov/car-race
//https://hpjansson.org/blag/2018/07/24/a-hash-table-re-hash/
//https://attractivechaos.wordpress.com/2018/01/13/revisiting-hash-table-performance/
//https://www.reddit.com/r/cpp/comments/auwbmg/hashmap_benchmarks_what_should_i_add/
//https://www.youtube.com/watch?v=M2fKMP47slQ
//https://yq.aliyun.com/articles/563053

//https://attractivechaos.wordpress.com/2008/08/28/comparison-of-hash-table-libraries/
//https://en.wikipedia.org/wiki/Hash_table
//https://tessil.github.io/2016/08/29/benchmark-hopscotch-map.html
//https://probablydance.com/2017/02/26/i-wrote-the-fastest-hashtable/
//https://martin.ankerl.com/2016/09/21/very-fast-hashmap-in-c-part-3/
//https://andre.arko.net/2017/08/24/robin-hood-hashing/
//http://www.ilikebigbits.com/2016_08_28_hash_table.html
//http://www.idryman.org/blog/2017/05/03/writing-a-damn-fast-hash-table-with-tiny-memory-footprints/

//https://1ykos.github.io/ordered_patch_map/#Figure%201
#if __cplusplus >= 199711L && HOOD_HASH
    #include "./tsl/robin_hood.h"       //https://github.com/martin/robin-hood-hashing/blob/master/src/include/robin_hood.h
#endif

#if __cplusplus >= 201103L || _MSC_VER >= 1600
    #define _CPP11_HASH    1
    #include "./phmap/phmap.h"           //https://github.com/greg7mdp/parallel-hashmap/tree/master/parallel_hashmap  
    #include "./tsl/robin_map.h"        //https://github.com/tessil/robin-map
    #include "./tsl/hopscotch_map.h"    //https://github.com/tessil/hopscotch-map
    #include "./tsl/robin_hood.h"       //https://github.com/martin/robin-hood-hashing/blob/master/src/include/robin_hood.h
    #include "./tsl/flat_hash_map.hpp"  //https://github.com/skarupke/flat_hash_map/blob/master/flat_hash_map.hpp
#if __cplusplus >= 201403L
    #define _CPP14_HASH   1
    #include "./tsl/bytell_hash_map.hpp"//https://github.com/skarupke/flat_hash_map/blob/master/bytell_hash_map.hpp
#endif
    #include "./tsl/robin_set.h"        //https://github.com/tessil/robin-map
#endif

#include <unordered_map>
#include <unordered_set>

#ifdef _WIN32
    # define CONSOLE "CON"
    # include <windows.h>
#else
    # define CONSOLE "/dev/tty"
    #include <unistd.h>
    #include <sys/resource.h>
#endif

#if _WIN32
typedef __int64 int64;
#else
typedef long int64;
#endif
struct RankItem;

#if TKey == 0
    typedef unsigned int         keyType;
    #define TO_KEY(i)   (keyType)i
    #define sKeyType    "int"
    #define KEY_INT
#elif TKey == 1
    typedef int64       keyType;
    #define TO_KEY(i)   (keyType)i
    #define sKeyType    "int64"
    #define KEY_INT
#else
    typedef std::string keyType;
    #define TO_KEY(i)   std::to_string(i)
    #define sKeyType    "string"
#endif

#if TVal == 0
    typedef int         valueType;
    #define TO_VAL(i)   i
    #define TO_SUM(i)   i
    #define sValueType  "int"
#elif TVal == 1
    typedef int64       valueType;
    #define TO_VAL(i)   i
    #define TO_SUM(i)   i
    #define sValueType  "int64"
#elif TVal == 2
    typedef std::string valueType;
    #define TO_VAL(i)   std::to_string(i)
    #define TO_SUM(i)   i.size()
    #define sValueType  "string"
#else
    typedef RankItem    valueType;
    #define TO_VAL(i)   i
    #define TO_SUM(i)   i.lScore
    #define sValueType  "RankItem"
#endif

emilib5::HashMap<std::string, std::string> show_name = {
//    {"stl_hash", "unordered_map"},

    {"emilib2", "emilib2"},
    {"emilib3", "emilib3"},
    {"emilib4", "emilib4"},
//    {"emilib5", "emilib5"},
//  {"emilib6", "emilib6"},
    {"emilib7", "emilib7"},

    {"martin", "martin flat"},
    {"phmap", "phmap flat"},

#if HOOD_HASH
    {"robin", "tessil robin"},
    {"hopsco", "tessil hopsco"},
    {"flat", "skarupk flat"},
#endif

//    {"byte", "skarupk byte"},
};

static int64 getTime()
{
#if _WIN32
    FILETIME ptime[4] = {0};
    GetThreadTimes(GetCurrentThread(), &ptime[0], &ptime[1], &ptime[2], &ptime[3]);
    return (ptime[2].dwLowDateTime + ptime[3].dwLowDateTime) / 10;
#elif __linux__ || __unix__
    struct rusage rup;
    getrusage(RUSAGE_SELF, &rup);
    long sec  = rup.ru_utime.tv_sec  + rup.ru_stime.tv_sec;
    long usec = rup.ru_utime.tv_usec + rup.ru_stime.tv_usec;
    return sec * 1000000 + usec;
#elif _WIN32
    return clock() * 1000ll;
#else
    return clock();
#endif
}

static std::map<std::string, long>  check_result;
static std::multimap <int64, std::string> func_time;
static std::map<std::string, int64> map_time;
static std::map<std::string, std::map<std::string, int64>> func_map_time;

#define COUNT_TIME(key)             map_time[show_name[key]] += getTime() - ts1; CHRCK_MAP_R()
#define AVE_TIME(ts, n)             1000 * (getTime() - ts) / (n)

static void check_mapfunc_result(const std::string& map_name, const std::string& func, long sum, long ts1, int size)
{
    if (check_result.find(func) == check_result.end()) {
        check_result[func] = sum;
    }
    else if (sum != check_result[func]) {
        printf("%s %s %ld != %ld\n", map_name.c_str(), func.c_str(), sum, check_result[func]);
    }

    auto& showname = show_name[map_name];
    auto timeuse = (getTime() - ts1);
    func_time.insert(std::pair<int64, std::string>(timeuse / 1000, showname));
    map_time[showname] += timeuse;
    func_map_time[func][showname] += timeuse;
}

static void set_func_time(std::map<std::string, std::map<std::string, int64_t>>& func_rank_time)
{
    for (auto& v : func_map_time) {
        for (auto& f : v.second) {
            func_rank_time[v.first][f.first] += f.second;
        }
    }
    func_map_time.clear();
}

static void dump_func(const std::string& func, const std::map<std::string, int64_t >& map_rtime)
{
    std::multimap <int64, std::string> functime;
    for (const auto& v : map_rtime)
        functime.insert(std::pair<int64_t, std::string>(v.second, v.first));

    puts(func.c_str());

    auto min = functime.begin()->first + 1;
    for (auto& v : functime)
        printf("   %-8ld     %-21s   %02d\n", v.first / 10000, v.second.c_str(), (int)((min * 100) / (v.first + 1)));
    putchar('\n');
}

static void dump_all(std::map<std::string, std::map<std::string, int64_t>>& func_rtime)
{
    for (const auto& v : func_rtime)
        dump_func(v.first, v.second);
}

template<class MAP>
void hash_iter(MAP& amap, const std::string& map_name, std::vector<keyType>& vList)
{
    if (show_name.find(map_name) != show_name.end())\
    {
        auto ts1 = getTime(); auto sum = 0;
        for (const auto& v : amap)
            sum += TO_SUM(v.second);
        for (auto it = amap.cbegin(); it != amap.cend(); ++it)
            sum ++;

        check_mapfunc_result(map_name, __FUNCTION__, sum, ts1, vList.size());
    }
}

template<class MAP>
void hash_reinsert(MAP& amap, const std::string& map_name, std::vector<keyType>& vList)
{
    if (show_name.find(map_name) != show_name.end())
    {
        size_t sum = 0;
        auto ts1 = getTime();
        for (const auto v : vList) {
            amap[v] = TO_VAL(1);
#if TVal < 2
            sum += amap[v];
#else
            sum += TO_SUM(amap[v]);
#endif

        }
        check_mapfunc_result(map_name, __FUNCTION__, sum, ts1, vList.size());
//        printf("    %12s    reinsert  %5ld ns, factor = %.2f\n", map_name.c_str(), AVE_TIME(ts1, vList.size()), amap.load_factor());
    }
}

template<class MAP>
void hash_insert2(MAP& amap, const std::string& map_name, std::vector<keyType>& vList)
{
    if (show_name.find(map_name) != show_name.end())
    {
        size_t sum = 0;
        auto ts1 = getTime();
        for (const auto v : vList) {
            amap.insert_unique(v, TO_VAL(0));
            sum ++;
        }

        check_mapfunc_result(map_name, "hash_emplace", sum, ts1, vList.size());
        printf("    %12s    %s  %5lu ns, factor = %.2f\n", "hash_insert2", map_name.c_str(), AVE_TIME(ts1, vList.size()), amap.load_factor());
    }
}

template<class MAP>
void hash_emplace(MAP& amap, const std::string& map_name, std::vector<keyType>& vList)
{
    if (show_name.find(map_name) != show_name.end())
    {
        size_t sum = 0;
        auto ts1 = getTime();
        for (const auto v : vList)
            sum += amap.emplace(v, TO_VAL(0)).second;

        check_mapfunc_result(map_name, __FUNCTION__, sum, ts1, vList.size());
        printf("    %12s    %s  %5lu ns, factor = %.2f\n", __FUNCTION__, map_name.c_str(), AVE_TIME(ts1, vList.size()), amap.load_factor());
    }
}

template<class MAP>
void hash_miss(MAP& amap, const std::string& map_name, std::vector<keyType>& vList)
{
    if (show_name.find(map_name) != show_name.end())\
    {
        auto n = vList.size();
        auto ts1 = getTime(); size_t sum = 0;
        for (size_t v = 1; v < 2*n; v++)
            sum += amap.count(TO_KEY(v));
        check_mapfunc_result(map_name, __FUNCTION__, sum, ts1, vList.size());
    }
}

template<class MAP>
void hash_erase(MAP& amap, const std::string& map_name, std::vector<keyType>& vList)
{
    if (show_name.find(map_name) != show_name.end())\
    {
        auto ts1 = getTime(); size_t sum = 0;
        for (const auto v : vList)
            sum += amap.erase(v);
        //func_time.insert(std::pair<int64, std::string>(AVE_TIME(ts1, vList.size() + sum % 2), show_name[map_name]));
        check_mapfunc_result(map_name, __FUNCTION__, sum, ts1, vList.size());
    }
}

template<class MAP>
void hash_find(MAP& amap, const std::string& map_name, std::vector<keyType>& vList)
{
    if (show_name.find(map_name) != show_name.end())
    {
        auto ts1 = getTime();
        size_t sum = 0;
        for (const auto v : vList) {
#ifdef KEY_INT
            sum += amap.count(v) + v;
#else
            sum += amap.count(v) + v.size();
#endif
        }
        check_mapfunc_result(map_name, __FUNCTION__, sum, ts1, vList.size());
    }
}

template<class MAP>
void hash_find2(MAP& amap, const std::string& map_name, std::vector<keyType>& vList)
{
    if (show_name.find(map_name) != show_name.end())
    {
        auto ts1 = getTime();
        size_t sum = 0;
        for (const auto v : vList)
            sum += amap.count(v);
        check_mapfunc_result(map_name, __FUNCTION__, sum, ts1, vList.size());
    }
}

template<class MAP>
void hash_clear(MAP& amap, const std::string& map_name, std::vector<keyType>& vList)
{
    if (show_name.find(map_name) != show_name.end())\
    {
        auto ts1 = getTime();
        size_t sum = amap.size();
        amap.clear();
        amap.clear();
        check_mapfunc_result(map_name, __FUNCTION__, sum, ts1, sum);
    }
}

template<class MAP>
void hash_copy(MAP& amap, const std::string& map_name, std::vector<keyType>& vList)
{
    if (show_name.find(map_name) != show_name.end())
    {
        size_t sum = 0;
        auto ts1 = getTime();
        MAP tmap = amap;
        amap = tmap;
        sum  = tmap.size();
        check_mapfunc_result(map_name, __FUNCTION__, sum, ts1, sum);
    }
}

//#define FUNC_RANK() for (auto& v : func_time) \ printf("   %-4ld     %-49s\n",  v.first, v.second.c_str()); putchar('\n'); func_time.clear();
#define FUNC_RANK()  func_time.clear();

#ifndef PACK
#define PACK 128
#endif
struct RankItem
{
    RankItem()
    {
        lScore = 0;
        lUid = iRank = 0;
        iUpdateTime = 0;
    }

    RankItem(int64 lUid1, int64 lScore1 = 0, int iTime = 0)
    {
        lUid   = lUid1;
        lScore = lScore1;
        iUpdateTime = iTime;
    }

/*:
    bool operator < (const RankItem& r) const
    {
        if (lScore != r.lScore)
            return lScore > r.lScore;
        else if (iUpdateTime != r.iUpdateTime)
            return iUpdateTime < r.iUpdateTime;
        return lUid < r.lUid;
    }

    bool operator == (const RankItem &r) const
    {
        if (lUid == r.lUid) {
            return true;
        }

        return false;
    }

    RankItem& operator += (const RankItem& r)
    {
        lScore += r.lScore;
        iRank   = r.iRank;
        //        iUpdateTime = time(0);
        return *this;
    }

    RankItem& operator += (int64 r)
    {
        lScore += r;
        return *this;
    }
**/
    int64 operator ()()
    {
        return lScore;
    }

    int64 lUid;
    int64 lScore;
    int  iUpdateTime;
    int  iRank;
#if PACK >= 24
    char data[(PACK - 24) / 8 * 8];
#else
    std::string data;
#endif
};

#if PACK >= 24
static_assert(sizeof(RankItem) == PACK, "PACK >=24");
#endif

static int ilog(int x, const int n = 2)
{
    int logn = 0;
    while (x / n) {
        logn ++;
        x /= n;
    }

    return logn;
}

#include <chrono>

static uint32_t get32rand()
{
    return (rand() ^ (rand() << 15u) ^ (rand() << 30u));
}

static int64 get64rand()
{
    return (((uint64_t)get32rand()) << 32) | get32rand();
}

static const std::array<char, 62> ALPHANUMERIC_CHARS = {
    '0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
    'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
    'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z'
};

static const std::int64_t SEED = 0;
static std::mt19937_64 generator(SEED);
std::uniform_int_distribution<std::size_t> rd_uniform(0, ALPHANUMERIC_CHARS.size() - 1);

static std::string get_random_alphanum_string(std::size_t size) {
    std::string str(size, '\0');
    for(std::size_t i = 0; i < size; i++) {
        str[i] = ALPHANUMERIC_CHARS[rd_uniform(generator)];
    }

    return str;
}

template<class RandomIt>
void shuffle(RandomIt first, RandomIt last)
{
    std::random_device rd;
    std::mt19937 g(rd());
    typedef typename std::iterator_traits<RandomIt>::difference_type diff_t;
    typedef std::uniform_int_distribution<diff_t> distr_t;
    typedef typename distr_t::param_type param_t;

    distr_t D;
    diff_t n = last - first;
    for (diff_t i = n-1; i > 0; --i) {
        using std::swap;
        swap(first[i], first[D(g, param_t(0, i))]);
    }
}

//https://en.wikipedia.org/wiki/Gamma_distribution#/media/File:Gamma_distribution_pdf.svg
//https://blog.csdn.net/luotuo44/article/details/33690179
static int buildTestData(int size, std::vector<keyType>& rankdata)
{
    rankdata.reserve(size);

#ifndef KEY_INT
    for (int i = 0; i < size; i++)
        rankdata.emplace_back(get_random_alphanum_string(rand() % 10 + 6));
    return 0;
#else

    sfc64 srng;
    auto flag = (int)rand() % 4 + 1;

    if (rand() % 100 > 20)
    {
        emilib6::HashSet<keyType> eset(size);
        for (int i = 0; ; i++) {
            auto key = TO_KEY(srng());
            if (eset.insert(key).second) {
                rankdata.emplace_back(key);
                if (rankdata.size() >= size)
                    break;
            }
        }
        flag = 0;
    }
    else
    {
        const int pow2 = 1 << ilog(size, 2);
        auto k = srng();
        for (int i = 1; i <= size; i ++) {
            k ++;
            if (flag == 2)
                k += (1 << 8) - 1;
            else if (flag == 3) {
                if (rand() % 256 == 0)
                    k += pow2 - rand() % 256;
            }
            else if (flag == 4) {
                if (rand() % 64 == 0)
                    k += pow2 / 8;
            }
            rankdata.emplace_back(k);
        }
    }

    printf("flag = %d\n", flag);
    return flag;
#endif
}

static int readFile(std::string fileName, int size)
{
    if (!freopen(fileName.c_str(), "rb", stdin)) {
        std::cerr << "Cannot open the File : " << fileName << std::endl;
        freopen(CONSOLE, "r", stdin);
        return -1;
    }

    auto ts1 = getTime();
    std::vector<int64> vList;
    vList.reserve(1038860);

    int64                uid, score, pid;
    if (size == 1) {
        while (scanf("%ld", &uid) == size) vList.emplace_back(uid);
    } else if (size == 2) {
        while (scanf("%ld %ld", &uid, &score) == size) vList.emplace_back(uid);
    } else if (size == 4) {
        int items = 0;
        while (scanf("%ld,%ld,%d,%ld", &uid, &pid, &items, &score) == size)
            vList.emplace_back(uid);
    }
    freopen(CONSOLE, "r", stdin);

    //    std::sort(vList.begin(), vList.end());
    int n = (int)vList.size();
    printf("\nread file %s  %ld ms, size = %zd\n", fileName.c_str(), getTime() - ts1, vList.size());

    emilib3::HashMap<int64, int> emap3(n);
    emilib2::HashMap<int64, int> emap2(n);
    emilib5::HashMap<int64, int> emap5(n);

#if _CPP11_HASH
    ska::flat_hash_map <int64, int, std::hash<int>> fmap;
    tsl::robin_map     <int64_t, int> rmap;
    tsl::hopscotch_map <int64_t, int> hmap;

    fmap.reserve(n); rmap.reserve(n); hmap.reserve(n);
#endif

    ts1 = getTime();    for (auto v : vList)        emap3[v] = 1; printf("    emilib3       insert  %4ld ms, size = %zd\n", getTime() - ts1, emap3.size());
    ts1 = getTime();    for (auto v : vList)        emap2[v] = 1; printf("    emilib2       insert  %4ld ms, size = %zd\n", getTime() - ts1, emap2.size());
    ts1 = getTime();    for (auto v : vList)        emap5[v] = 1; printf("    emilib3       insert  %4ld ms, size = %zd\n", getTime() - ts1, emap5.size());

#if _CPP11_HASH
    ts1 = getTime();    for (auto v : vList)        fmap[v] = 1;  printf("    flatmap       insert  %4ld ms, size = %zd\n", getTime() - ts1, fmap.size());
    ts1 = getTime();    for (auto v : vList)        rmap[v] = 1;  printf("    robinmap      insert  %4ld ms, size = %zd\n", getTime() - ts1, rmap.size());
    ts1 = getTime();    for (auto v : vList)        hmap[v] = 1;  printf("    oodflat      insert  %4ld ms, size = %zd\n", getTime() - ts1, hmap.size());
#endif

    printf("\n");
    return 0;
}

static int HashMapTest(int n, int max_loops = 1234567)
{
    emilib4::HashMap <keyType,int> emap5;
    emilib2::HashMap <keyType,int> emap2;

#if _CPP11_HASH
    robin_hood::unordered_flat_map <keyType, int> umap;
#else
    std::unordered_map<keyType,int> umap;
#endif

    const auto step = n % 2 + 1;
    emap2.reserve(n / 4); emap5.reserve(n / 2);
    umap.reserve(n/2);

    for (int i = 1; i < n * step; i += step) {
        auto ki = TO_KEY(i);
        emap5[ki] = umap[ki] = 0;
        emap2[ki] = 0;
    }

    int loops = 0;
    while (loops++ < max_loops) {
        const int type = rand() % 100;
        auto rid  = n ++;
        auto id   = TO_KEY(rid);
        if (type <= 50 || emap2.size() < 100) {
            emap2[id] += type; emap5[id] += type; umap[id]  += type;

            assert(emap2[id] == umap[id]); assert(emap5[id] == umap[id]);
            assert(emap2.size() == umap.size()); assert(emap5.size() == umap.size());
        }
        else if (type < 70) {
            if (n % 4 == 0)
                id = umap.begin()->first;
            else if (n % 8 == 0)
                id = emap2.begin()->first;
            else
                id = emap5.begin()->first;

            assert(emap2.size() == umap.size()); assert(emap5.size() == umap.size());
            if (umap.count(id) == 1) {
                emap5.erase(emap5.find(id));
                umap.erase(id);
                emap2.erase(id);

                assert(emap5.count(id) == umap.count(id));
                assert(emap2.count(id) == umap.count(id));
                assert (emap2.size() == umap.size());
                assert (emap5.size() == umap.size());
            }
        }
        else if (type < 90) {
            auto it = emap5.begin();
            for (int i = n % 32; i > 0; i--)
                it ++;
            id = it->first;
            umap.erase(id);
            emap2.erase(id);
            emap5.erase(it);
            assert(emap5.size() == umap.size());
            assert(emap2.count(id) == 0);
            assert(emap5.count(id) == umap.count(id));
            assert (emap2.size() == umap.size());
        }
        else if (type < 98) {
            if (emap5.count(id) == 0) {
                const auto vid = (int)rid;
                emap5.insert_unique(id, vid);
                assert(emap5.count(id) == 1);
                assert(emap2.count(id) == 0);

                //if (id == 1043)
                emap2.insert(id, vid);
                assert(emap2.count(id) == 1);

                umap[id] = emap2[id];
                assert(umap[id] == emap2[id]);
                assert(umap[id] == emap5[id]);
            } else {
                umap.erase(id);
                emap2.erase(id);
                emap5.erase(id);
#ifdef KEY_INT
                uint64_t sum1 = 0, sum2 = 0, sum3 = 0;
                for (auto v : umap)
                    sum1 += v.first * v.second;
                for (auto v : emap2)
                    sum2 += v.first * v.second;
                for (auto v : emap5)
                    sum3 += v.first * v.second;
                assert(sum1 == sum2);
                assert(sum1 == sum3);
#endif
            }
        }
        if (loops % 4096 == 0) {
            printf("%zd %d\r", emap2.size(), loops);
        }
    }

    printf("\n");
    return 0;
}

template<class MAP>
int benOneMap(MAP& hmap, const std::string& map_name, std::vector<keyType> vList)
{
    if (show_name.find(map_name) == show_name.end())
        return 80;

#if IR == 0
    const auto iRation = vList.size() / 1;
    hmap.reserve(iRation);
#else
    const auto iRation = vList.size() / IR;
    hmap.reserve(iRation);
#endif

    func_time.clear();

    hash_emplace(hmap, map_name, vList);
//  shuffle(vList.begin(), vList.end());
    hash_find(hmap, map_name, vList);

//    shuffle(vList.begin(), vList.end());
    hash_miss(hmap, map_name, vList);

#ifdef KEY_INT
    for (int v = 0; v < (int)vList.size() / 2; v ++)
        vList[v] ++;
#else
    sfc64 rng(vList.size());
    for (int v = 0; v < (int)vList.size() / 2; v ++)
        vList[v][0] = v % 10 + '0';
#endif

//    shuffle(vList.begin(), vList.end());
    hash_erase(hmap, map_name, vList);

    hash_find2(hmap, map_name, vList);
//    shuffle(vList.begin(), vList.end());
    hash_reinsert(hmap, map_name, vList);

    auto lf = (int)(hmap.load_factor() * 100);

#if UF
    hash_copy(hmap, map_name, vList);
    hash_iter(hmap, map_name, vList);
    hash_clear(hmap, map_name, vList);
#endif

    return lf;
}

static void benchMarkHashMap2(int n)
{
    if (n < 10000)
        n = 123456;
    printf("%s n = %d, keyType = %s, valueType = %s\n", __FUNCTION__, n, sKeyType, sValueType);

    typedef std::hash<keyType>        std_func;

#ifdef HOOD_HASH
    typedef robin_hood::hash<keyType> hood_func;
    using hash_func = hood_func;
#else
    using hash_func = std_func;
#endif

    auto iload = 0;
    float lf = 0.90f;

    check_result.clear();
    map_time.clear();
    func_map_time.clear();

    std::vector<keyType> vList;
    auto step = buildTestData(n, vList);

    { phmap::flat_hash_map <keyType, valueType, hash_func> emap; emap.max_load_factor(lf); benOneMap(emap, "phmap", vList); }

    { emilib2::HashMap <keyType, valueType, hash_func> emap; emap.max_load_factor(lf); benOneMap(emap, "emilib2", vList); }

    { emilib3::HashMap <keyType, valueType, hash_func> emap; emap.max_load_factor(lf); iload = benOneMap(emap, "emilib3", vList); }

    { emilib4::HashMap <keyType, valueType, hash_func> emap; emap.max_load_factor(lf); benOneMap(emap, "emilib4", vList); }

    { std::unordered_map<keyType, valueType, hash_func> umap; umap.max_load_factor(lf); benOneMap(umap, "stl_hash", vList); }

#if _CPP14_HASH
    { ska::bytell_hash_map <keyType, valueType, hash_func > bmap; bmap.max_load_factor(lf); benOneMap(bmap, "byte", vList); }
#endif

#if _CPP11_HASH
    { robin_hood::unordered_flat_map <keyType, valueType, hash_func> mmap; benOneMap(mmap, "martin", vList); }

    { ska::flat_hash_map <keyType, valueType, hash_func> fmap; fmap.max_load_factor(lf); benOneMap(fmap, "flat", vList); }

    { tsl::hopscotch_map <keyType, valueType, hash_func> hmap; hmap.max_load_factor(lf); benOneMap(hmap, "hopsco", vList); }

    { tsl::robin_map     <keyType, valueType, hash_func> rmap; rmap.max_load_factor(lf); benOneMap(rmap, "robin", vList); }
#endif

    {
//        emilib7::HashMap <keyType, valueType, hash_func> emap; emap.max_load_factor(lf); benOneMap(emap, "emilib7", vList);
    }

    static int tcase = 1;
    printf("\n %d ======== n = %d, flag = %d hash_map(%.2lf) ========\n", tcase, n, step, iload / 100.0);
    std::multimap <int64, std::string> time_map;
    for (auto& v : map_time)
        time_map.insert(std::pair<int64, std::string>(v.second, v.first));

    const auto last  = double(time_map.rbegin()->first);
    const auto first = double(time_map.begin()->first);
    if (first < 10 || last < 9)
        return;

    static std::map<std::string,int64> rank;
    static std::map<std::string,int64> rank_time;
    static std::map<std::string, std::map<std::string, int64>> func_rank_time;

    auto it0 = time_map.begin();
    auto it1 = *(it0++);
    auto it2 = *(it0++);
    auto it3 = *(it0++);

    constexpr auto base1 = 300000000;
    constexpr auto base2 =      20000;

    if (it1.first == it3.first) {
        rank[it1.second] += base1 / 3;
        rank[it2.second] += base1 / 3;
        rank[it3.second] += base1 / 3;
    } else if (it1.first == it2.first) {
        rank[it1.second] += base1 / 2;
        rank[it2.second] += base1 / 2;
        rank[it3.second] += 1;
    } else {
        rank[it1.second] += base1;
        if (it2.first == it3.first) {
            rank[it2.second] += base2 / 2;
            rank[it3.second] += base2 / 2;
        } else {
            rank[it2.second] += base2;
            rank[it3.second] += 1;
        }
    }

    set_func_time(func_rank_time);
    for (auto& v : time_map) {
        rank_time[v.second] += (int)(first * 100 / v.first);
        printf("%5ld   %13s   (%4.2lf %6.1lf%%)\n", v.first * 1000ll / n, v.second.c_str(), last * 1.0 / v.first, first * 100.0 / v.first);
    }

    if (tcase++ % 5 == 0) {
        printf("--------------------------------%s lf = %d--------------------------------\n", __FUNCTION__, iload);
        dump_all(func_rank_time);
        for (auto& v : rank)
            printf("%13s %10ld  %4.1lf %4.1lf %4ld\n", v.first.c_str(), v.second, v.second / (double)(base1), (v.second / (base2 / 2) % 1000) / 2.0, (v.second % (base2 / 2)));
        for (auto& v : rank_time)
            printf("%13s %4ld\n", v.first.c_str(), v.second / (tcase - 1));
#if _WIN32
        //        std::this_thread::sleep_for(std::chrono::milliseconds(5000));
#else
        usleep(1000*4000);
#endif
        printf("--------------------------------------------------------------------\n");
        return;
    }

    printf("=======================================================================\n\n");
}

void testSynax()
{
    emilib5::HashMap <std::string, std::string> mymap =
    {
        {"house","maison"},
        {"apple","pomme"},
        {"tree","arbre"},
        {"book","livre"},
        {"door","porte"},
        {"grapefruit","pamplemousse"}
    };

    std::vector<std::pair<std::string, std::string>> kv =
    {
        {"house2","maison"},
        {"apple2","pomme"},
        {"tree2","arbre"},
        {"book2","livre"},
        {"door2","porte"},
        {"grapefruit2","pamplemousse"}
    };

    for (auto& item : kv) {
        //mymap.insert_unique(item);
        mymap.emplace(item.first, item.second);
        mymap.insert(item);
        //mymap.emplace(std::move(item));
    }
    decltype(mymap) mymap2;
    mymap2.insert2(kv.begin(), kv.end());
    mymap2.insert2(kv.begin(), kv.end());
    mymap2 = mymap;
    //mymap2.reserve(mymap.bucket_count() * mymap.max_load_factor());

    auto nbuckets = mymap2.bucket_count();
    std::cout << "mymap has " << nbuckets << " buckets. and size = " << mymap2.size() << std::endl;
    for (unsigned i = 0; i < nbuckets; ++i) {
        std::cout << "bucket #" << i << " has " << mymap2.bucket_size(i) << " elements.\n";
    }

    for (auto& x : mymap2) {
        std::cout << "Element [" << x.first << ":" << x.second << "]";
        std::cout << " is in bucket #" << mymap2.bucket(x.first) << std::endl;
    }

    emilib6::HashSet<std::string> eset = {"abab",
        "cbacdcbc",
        "bcabc",
        "abc",
        "abba",
        "aaaaa",
        "bbbbaaab",
        "abada",
        "abadcea",
    };

    decltype(eset) eset2(std::move(eset));
    eset2.swap(eset);
    eset2 = eset;

    emilib3::HashMap <keyType, valueType> mymap_oper;
    emilib2::HashMap <keyType, valueType> mymap_insu;
    emilib5::HashMap <keyType, valueType> mymap_inse;

#if 0
    float load_factor = (36 + rand() % 64) / 100.0;
    mymap_oper.max_load_factor(load_factor);
    mymap_insu.max_load_factor(0.5);
    mymap_inse.max_load_factor(load_factor);
#endif

    std::vector<keyType> vList;
    int n = 123456 + get32rand() % (1 << 20);
    //n = (1 << 20) * load_factor;
    buildTestData(n, vList);

#if 1
    const int iRation = n / 4;// mymap_oper.size();
    mymap_oper.reserve(iRation);
    mymap_inse.reserve(iRation);
    mymap_insu.reserve(iRation);
#endif

    //auto ts1 = getTime();
    //    ts1 = getTime(); hash_insert(mymap_oper)    printf("insert []        %5ld ms, factor = %.2f\n", getTime() - ts1, mymap_oper.load_factor());
    //    ts1 = getTime(); hash_emplace(mymap_inse)   printf("insert nomal     %5ld ms, factor = %.2f\n", getTime() - ts1, mymap_inse.load_factor());
    //    ts1 = getTime(); hash_insert2(mymap_insu)   printf("insert_unique    %5ld ms, factor = %.2f\n", getTime() - ts1, mymap_insu.load_factor());

    //printf("testSynax load_factor = %.2f\n", load_factor);
}

int main(int argc, char* argv[])
{
    srand((unsigned)time(0));
    auto n = 1234567 + rand() % 1234567;

//    testSynax();

    printf("./test n load_factor (key=%s,value=%s)\n", sKeyType, sValueType);
    double load_factor = 0.1f;
#if LF
    load_factor = LF / 100.0;
#endif

    if (argc > 1 && argv[1][0] >= '0' && argv[1][0] <= '9')
        n = atoi(argv[1]);
    if (argc > 2)
        load_factor = atoi(argv[2]) / 100.0;

//    HashMapTest(n, 234567);

#ifdef TREAD
    //readFile("./uid_income.txt", 1);
    //readFile("./pid_income.txt", 1);
    //readFile("./uids.csv", 2);
    readFile("./item.log", 4);
#endif

    auto nows = time(0);

    while (true) {
#if INP
        char ccmd[257];
        printf(">> ");
        if (fgets(ccmd, 255, stdin)) {
            auto n = atol(ccmd);
            if (n == 0)
                n = (get32rand() >> 9) + 14567;
            else if (n < 0)
                break;
            if (load_factor > 0.4 && load_factor < 0.99) {
                int log2 = ilog(n, 2);
                n = int((1 << log2) * load_factor) + rand() % (1 << 10);
            }

            benchMarkHashMap(n);
        }
#else
        n = (get32rand() >> 10) + 12345;
        auto pow2 = 1 << ilog(n, 2);

#if LF > 100
        n = pow2 * (60 + rand()% 35) / 1000;
#endif

        if (load_factor > 0.4 && load_factor < 0.95)
            n = int(pow2 * load_factor) + rand() % (1 << 10);

//        if (rand() % 2 == 0)
            benchMarkHashMap2(n);
//        else
//            benchMarkHashMap(n);
#endif

#if TP || GCOV
        if (time(0) > nows + 20)
            break;
#endif

        if (rand() % 64 == 1)
            HashMapTest(n, rand() * rand() % 123456);
    }

    return 0;
}

