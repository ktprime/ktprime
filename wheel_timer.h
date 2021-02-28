#pragma once

#include <stdint.h>
#include <vector>

#include "net/quic/quartc/quartc_task_runner_interface.h"
typedef net::QuartcTaskRunnerInterface::Task* QTask;
typedef net::QuartcTaskRunnerInterface QuartcTaskRunnerInterface;

#if 1
    #include "hash_table5.hpp"
    #define HASH_MAP emhash5::HashMap
#else
    #include "robin_hood.h"
    #define HASH_MAP robin_hood::unordered_map
#endif



// timer queue implemented with hashed hierarchical wheel.
//
// inspired by linux kernel, see links below
// https://git.kernel.org/pub/scm/linux/kernel/git/stable/linux-stable.git/tree/kernel/timer.c?h=v3.2.98
//
// We model timers as the number of ticks until the next
// due event. We allow 32-bits of space to track this
// due interval, and break that into 4 regions of 8 bits.
// Each region indexes into a bucket of 256 lists.
//
// complexity:
//      AddTimer   CancelTimer   PerTick
//       O(1)       O(1)          O(1)
//

#define WHEEL_BUCKETS 2
enum TIME_WHEEL
{
    TVN_BITS = 6,                   // time vector level shift bits
    TVR_BITS = 10,                  // timer vector shift bits
    TVN_SIZE = (1 << TVN_BITS),     // wheel slots of level vector
    TVR_SIZE = (1 << TVR_BITS),     // wheel slots of vector
    TVN_MASK = (TVN_SIZE - 1),      //
    TVR_MASK = (TVR_SIZE - 1),      // near vector size
    TIME_UNIT = 1,                  // centisecond, i.e. 1/1000 second

    TIME_INDEX1 = 1 << (TVR_BITS + 1 * TVN_BITS),
    TIME_INDEX2 = 1 << (TVR_BITS + 2 * TVN_BITS),
    TIME_INDEX3 = 1 << (TVR_BITS + 3 * TVN_BITS),

    ALL_BITS  = TVR_BITS + WHEEL_BUCKETS * TVN_BITS,

    MAX_TVAL  = (uint64_t)((1ULL << ALL_BITS) - 1),
    MAX_SECOS = (uint32_t)(MAX_TVAL / (1000 / TIME_UNIT)),
    MAX_MINUS = (uint32_t)(MAX_SECOS / 60),
    MAX_HOURS = (uint32_t)(MAX_SECOS / 3600),
};

#if __cplusplus > 201702 || __GUNC__
static_assert(ALL_BITS  < 32);
static_assert(MAX_MINUS >= 10 && MAX_HOURS < 24);
static_assert(WHEEL_BUCKETS >= 1 && WHEEL_BUCKETS <= 4);//one hour - one day
#endif

class WheelTimer
{
public:
    static constexpr int INVALID_NODE_TIMERID = -1;
    static constexpr int INVALID_NODE_SLOTID  = -2;
    static constexpr int FREE_LIST_CAPACITY   = 1024;
    static constexpr int ALLOC_SIZE = 1024;

    struct TimerNode
    {
        int64_t id;         // timer id
        int64_t expire;     // jiffies
        int     slot;       // near index
        QTask   cb;
    };

    typedef std::vector<TimerNode*> TimerList;

public:
    WheelTimer(int64_t current_ms);
    ~WheelTimer();

    int64_t Schedule(int64_t deadline_ms, QTask cb);
    bool UpdateExpire(int64_t deadline_ms, int64_t timer_id);
    bool Cancel(const int64_t timer_id);
    bool Erase(const int64_t timer_id);
    int Update(int64_t now_ms);
    int64_t nextTimer(int slots);
    int Size() const  { return ref_.size(); }
    int64_t MaxId() const { return counter_;  }

private:
    int tick();
    void addTimerNode(TimerNode* node);
    int execute();
    bool cascade(int bucket);

    void clearAll();
    TimerNode* allocNode();
    void freeNode(TimerNode* node);
    void freeList(TimerList& list);
    int64_t nextCounter() { return ++counter_; }

private:
    int64_t counter_ = 0;
    int64_t current_ = 0;
    int64_t jiffies_ = 0;

    HASH_MAP<int64_t, TimerNode*> ref_;
    TimerList free_list_;
    TimerList near_[TVR_SIZE];
    TimerList buckets_[WHEEL_BUCKETS][TVN_SIZE];
    std::vector<TimerNode*> alloc_list_;
    int alloc_size_;
};
