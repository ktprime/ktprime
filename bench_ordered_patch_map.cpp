#include <bitset>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <string>
#include <unordered_map>
#ifdef PATCHMAP
#include "ordered_patch_map.hpp"
#endif
//#include "doubling_patchmap.hpp"
//#include "unordered_map.hpp"
#ifdef SPARSE_PATCHMAP
#include "sparse_patchmap.hpp"
#endif
#ifdef FLATMAP
#include "tsl/flat_hash_map.hpp"
#endif
#ifdef BYTELL
#include "tsl/bytell_hash_map.hpp"
#endif
#ifdef SPARSEPP
#include "sparsepp/spp.h"
#endif
#ifdef SPARSEMAP
#include <sparsehash/sparse_hash_map>
#endif
#ifdef DENSEMAP
#include <sparsehash/dense_hash_map>
#endif
#ifdef KHASH
#include "khash.h"
#endif
#ifdef EMIMAP
#include "hash_table52.hpp"
#endif

#if HOODMAP
#include "tsl/robin_hood.h"
#endif

#if PHMAP
#include "phmap/phmap.h"
#endif

#include "tsl/robin_map.h"

#include <stdio.h>
#if defined (__has_include) && (__has_include(<x86intrin.h>))
#include <x86intrin.h>
#endif

#include <emmintrin.h>


//#include <proc/readproc.h>

//using wmath::ordered_patch_map;
#if 0
using wmath::reflect;
using wmath::rol;
using wmath::ror;
#endif

using std::cout;
using std::endl;
using std::bitset;
using std::setw;
using std::get;
using std::stoi;
using std::max;

#ifdef KHASH
KHASH_MAP_INIT_INT(32, uint32_t)
#endif

static int log2(uint64_t x)
{
	int logn = 0;
	while (x / 2) {
		logn++;
		x /= 2;
	}

	return logn;
}


uint32_t const clmul_mod(const uint32_t& i,const uint32_t& j){
    __m128i I{}; I[0]^=i;
    __m128i J{}; J[0]^=j;
    __m128i M{}; M[0]^=1073877003u;
    __m128i X = _mm_clmulepi64_si128(I,J,0);
    __m128i X0{};X0[0]^=X[0]&(~uint32_t(0));
    __m128i A = _mm_clmulepi64_si128(X0,M,0);
    __m128i A0{};A0[0]^=A[0]&(~uint32_t(0));
    __m128i B = _mm_clmulepi64_si128(A0,M,0);
    return A[0]^(A[0]>>32)^(B[0]>>32)^X[0]^(X[0]>>32);
  }

uint32_t gen_rand(uint32_t i){
  //return i;
  //return wmath::rol(i*uint32_t(3061963241ul),16)*uint32_t(3107070805ul);
  return clmul_mod(uint32_t(i*3061963241ul),uint32_t(3107070805ul));
}

/*
 * Author:  David Robert Nadeau
 * Site:    http://NadeauSoftware.com/
 * License: Creative Commons Attribution 3.0 Unported License
 *          http://creativecommons.org/licenses/by/3.0/deed.en_US
 */

#if defined(_WIN32)
#include <windows.h>
#include <psapi.h>

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
#include <unistd.h>
#include <sys/resource.h>

#if defined(__APPLE__) && defined(__MACH__)
#include <mach/mach.h>

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
#include <fcntl.h>
#include <procfs.h>

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
#include <stdio.h>

#endif

#else
#error "Cannot define getPeakRSS( ) or getCurrentRSS( ) for an unknown OS."
#endif





/**
 * Returns the peak (maximum so far) resident set size (physical
 * memory use) measured in bytes, or zero if the value cannot be
 * determined on this OS.
 */
size_t getPeakRSS( )
{
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
    return (size_t)info.PeakWorkingSetSize;

#elif (defined(_AIX) || defined(__TOS__AIX__)) || (defined(__sun__) || defined(__sun) || defined(sun) && (defined(__SVR4) || defined(__svr4__)))
    /* AIX and Solaris ------------------------------------------ */
    struct psinfo psinfo;
    int fd = -1;
    if ( (fd = open( "/proc/self/psinfo", O_RDONLY )) == -1 )
        return (size_t)0L;      /* Can't open? */
    if ( read( fd, &psinfo, sizeof(psinfo) ) != sizeof(psinfo) )
    {
        close( fd );
        return (size_t)0L;      /* Can't read? */
    }
    close( fd );
    return (size_t)(psinfo.pr_rssize * 1024L);

#elif defined(__unix__) || defined(__unix) || defined(unix) || (defined(__APPLE__) && defined(__MACH__))
    /* BSD, Linux, and OSX -------------------------------------- */
    struct rusage rusage;
    getrusage( RUSAGE_SELF, &rusage );
#if defined(__APPLE__) && defined(__MACH__)
    return (size_t)rusage.ru_maxrss;
#else
    return (size_t)(rusage.ru_maxrss * 1024L);
#endif

#else
    /* Unknown OS ----------------------------------------------- */
    return (size_t)0L;          /* Unsupported. */
#endif
}

/**
 * Returns the current resident set size (physical memory use) measured
 * in bytes, or zero if the value cannot be determined on this OS.
 */
size_t getCurrentRSS( )
{
#if defined(_WIN32)
    /* Windows -------------------------------------------------- */
    PROCESS_MEMORY_COUNTERS info;
    GetProcessMemoryInfo( GetCurrentProcess( ), &info, sizeof(info) );
    return (size_t)info.WorkingSetSize;

#elif defined(__APPLE__) && defined(__MACH__)
    /* OSX ------------------------------------------------------ */
    struct mach_task_basic_info info;
    mach_msg_type_number_t infoCount = MACH_TASK_BASIC_INFO_COUNT;
    if ( task_info( mach_task_self( ), MACH_TASK_BASIC_INFO,
        (task_info_t)&info, &infoCount ) != KERN_SUCCESS )
        return (size_t)0L;      /* Can't access? */
    return (size_t)info.resident_size;

#elif defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
    /* Linux ---------------------------------------------------- */
    long rss = 0L;
    FILE* fp = NULL;
    if ( (fp = fopen( "/proc/self/statm", "r" )) == NULL )
        return (size_t)0L;      /* Can't open? */
    if ( fscanf( fp, "%*s%ld", &rss ) != 1 )
    {
        fclose( fp );
        return (size_t)0L;      /* Can't read? */
    }
    fclose( fp );
    return (size_t)rss * (size_t)sysconf( _SC_PAGESIZE);

#else
    /* AIX, BSD, Solaris, and Unknown OS ------------------------ */
    return (size_t)0L;          /* Unsupported. */
#endif
}

int main(int argc, char** argv){
  std::ios_base::sync_with_stdio(false);
  if (argc<2) return 1;
  const size_t N = stoi(argv[1]);
  std::random_device urand;
  std::minstd_rand mr(urand());
  std::uniform_int_distribution<uint32_t> u32distr(1u<<30);
  double initial_memory;
  {
  //  struct proc_t usage;
 //   look_up_our_self(&usage);
//    initial_memory = usage.vsize;
	  initial_memory = getCurrentRSS();
  }
  double base_time;
  int sand = 0;
  {
    auto start = std::clock();
    for (size_t i=0;i!=N;++i){
      sand+=gen_rand(i);
    }
    base_time = std::clock()-start;
    base_time/=N;
  }
#ifdef PATCHMAP
  wmath::ordered_patch_map<uint32_t,uint32_t> test;
#endif
#ifdef SPARSE_PATCHMAP
  wmath::sparse_patchmap<uint32_t,uint32_t> test;
#endif
#ifdef WMATH_ROBIN_MAP
  wmath::robin_map<uint32_t,uint32_t> test;
#endif
#ifdef HOODMAP
  robin_hood::unordered_map<uint32_t,uint32_t> test;
#endif
#ifdef SPARSEPP
  spp::sparse_hash_map<uint32_t,uint32_t> test;
#endif
#ifdef UNORDERED_MAP
  std::unordered_map<uint32_t,uint32_t> test;
#endif
#ifdef MAP
  std::map<size_t,size_t> test;
#endif
#ifdef FLATMAP
  ska::flat_hash_map<uint32_t,uint32_t> test;
#endif
#ifdef BYTELL
  ska::bytell_hash_map<uint32_t,uint32_t> test;
#endif
#ifdef SPARSEMAP
  google::sparse_hash_map<uint32_t,uint32_t> test;
#endif
#ifdef DENSEMAP
  google::dense_hash_map<uint32_t,uint32_t> test;
  test.set_empty_key(0);
  test.set_deleted_key(~uint32_t(0));
#endif
#ifdef KHASH
  khash_t(32) *h = kh_init(32);
  int ret, is_missing;
  khiter_t k;
#elif EMIMAP
   emilib2::HashMap<uint32_t, uint32_t> test;
#elif HOODMAP
  robin_hood::unordered_map<uint32_t, uint32_t> test;
#elif ROBINMAP
   tsl::robin_map<uint32_t,uint32_t> test;
#elif PHMAP
  phmap::flat_hash_map<uint32_t, uint32_t> test;
#endif

  std::uniform_int_distribution<size_t> distr;
  double counter = 0;
  double acc_insert = 0;
  double acc_find = 0;
  double acc_delete = 0;
  double acc_not_find = 0;
  double typical_insert_time = 0;
  double typical_delete_time = 0;
  double typical_find_time = 0;
  double typical_not_find_time = 0;
  double typical_memory = 0;
  auto start = std::clock();
  for (size_t i=0;i!=N/4096;++i) {
    auto start = std::clock();
    for (size_t j=0;j!=4096;++j) {
      const uint32_t n = gen_rand(2*(i*4096+j));
#ifdef KHASH
      k = kh_put(32, h, n, &ret);
      kh_value(h, k) = n;
#else
      test[n]=n;
#endif
    }
    typical_insert_time+=(acc_insert+=(std::clock()-start-base_time))/(i+1);
//    struct proc_t usage;
//    look_up_our_self(&usage);

    typical_memory+=(getCurrentRSS() - initial_memory )/(i+1);
    mr.seed(i+7ul);
    std::uniform_int_distribution<size_t> distr(0,i);
    const size_t l0 = distr(mr);
    start = std::clock();
    for (size_t j=0;j!=4096;++j) {
      const uint32_t n = gen_rand(2*(l0*4096+j));
#ifdef KHASH
      k = kh_get(32, h, n);
      sand += k;
#else
      sand+=test.count(n);
#endif
    }
    typical_find_time+=(acc_find+=(std::clock()-start-base_time))/(i+1);
    const size_t l1 = distr(mr);
    start = std::clock();
    for (size_t j=0;j!=4096;++j) {
      const uint32_t n = gen_rand(2*(l1*4096+j)+1);
#ifdef KHASH
      k = kh_get(32, h, n);
      sand += k;
#else
      sand+=test.count(n);
#endif
    }
    typical_not_find_time+=(acc_not_find+=(std::clock()-start-base_time))/(i+1);
    const size_t l2 = distr(mr);
    start = std::clock();
    for (size_t j=0;j!=4096;++j){
      const uint32_t n = gen_rand(2*(l2*4096+j));
#ifdef KHASH
      k = kh_get(32, h, n);
      if (k!=kh_end(h)) kh_del(32, h, k);
#else
      test.erase(n);
#endif
    }
    typical_delete_time+=(acc_delete+=(std::clock()-start-base_time))/(i+1);
    for (size_t j=0;j!=4096;++j){
      const uint32_t n = gen_rand(2*(l2*4096+j));
#ifdef KHASH
      k = kh_put(32, h, n, &ret);
      kh_value(h, k) = n;
#else
      test[n]=n;
#endif
    }
    if (0 == (i & (i - 1)) && i > 20) {
      cout <<  
		  i << "   "
		  <<  typical_memory/(i*4096)        << " "
           << typical_insert_time/(i*4096)   << " "
           << typical_delete_time/(i*4096)   << " "
           << typical_find_time/(i*4096)     << " "
           << typical_not_find_time/(i*4096) << endl;
    }
  }
  cout << typical_memory/N        << " "
       << typical_insert_time/N   << " "
       << typical_delete_time/N   << " "
       << typical_find_time/N     << " "
       << typical_not_find_time/N << endl;
#ifdef KHASH
  kh_destroy(32, h);
#endif
//#ifdef PATCHMAP
//  std::cerr << test.size() << " " << test.test_size() << endl;
//#endif
#ifdef UNORDERED_MAP
  std::cerr << test.size() << endl;
#endif
  return sand;
}
