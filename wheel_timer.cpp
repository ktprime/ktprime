#include "wheel_timer.h"

//#define FAST_REMOVE_TIMER  1

WheelTimer::WheelTimer(int64_t current)
    : current_(current / TIME_UNIT)
{
    ref_.reserve(8);            // reserve a little space
    free_list_.reserve(FREE_LIST_CAPACITY);
    jiffies_ = current_;
    alloc_size_ = 0;
}

WheelTimer::~WheelTimer()
{
    clearAll();
}

void WheelTimer::clearAll()
{
    for (auto ptr : alloc_list_)
        free(ptr);
    alloc_list_.clear();
    alloc_size_ = 0;
    free_list_.clear();
    ref_.clear();
}

WheelTimer::TimerNode* WheelTimer::allocNode()
{
    TimerNode* node;
    if (free_list_.size() > 0)
    {
        node = free_list_.back();
        free_list_.pop_back();
    }
    else
    {
       if (alloc_size_ == 0)
       {
	   alloc_size_ = ALLOC_SIZE;
	   alloc_list_.emplace_back((TimerNode*) malloc(alloc_size_ * sizeof(TimerNode)) );
       }
       node = alloc_list_.back() + (-- alloc_size_);
    }
    return node;
}

void WheelTimer::freeList(TimerList& list)
{
    free_list_.insert(free_list_.end(), list.begin(), list.end());
    list.clear();
}

void WheelTimer::freeNode(TimerNode* node)
{
    free_list_.emplace_back(node);
}

// Do lazy cancellation
bool WheelTimer::Erase(const int64_t timer_id)
{
    const auto itnode = ref_.find(timer_id);
    if (itnode == ref_.end())
        return false;

    auto node = itnode->second;
    ref_._erase(itnode);
    node->expire = 0;
    return true;
}

// Do lazy cancellation, so we can effectively use vector as container of timer nodes
bool WheelTimer::Cancel(const int64_t timer_id)
{
    auto nodeptr = ref_.try_get(timer_id);
    if (!nodeptr)
        return false;

    auto node = *nodeptr;
    if (node->expire <= 0)
        return false;

    node->expire = -node->expire;

#ifdef FAST_REMOVE_TIMER
    ref_._erase(timer_id);
    if (node->slot != INVALID_NODE_SLOTID) {
        auto& near = near_[node->slot];
        if (near.size() == 1) {
            assert(near[0]->id == timer_id);
            near.clear();
            return true;
        }
        for (auto it = near.begin(); it != near.end(); it++) {
            if ((*it)->id == timer_id) {
                near.erase(it);
                return true;
            }
        }
    }
#endif

    return true;
}

bool WheelTimer::UpdateExpire(int64_t deadline_ms, int64_t timer_id)
{
    auto nodeptr = ref_.try_get(timer_id);
    if (nodeptr) {
        auto& expire = (*nodeptr)->expire;
        if (std::abs(expire) <= deadline_ms / TIME_UNIT) {
            expire = deadline_ms / TIME_UNIT;
            return true;
        }
        ref_.erase(timer_id);
        expire = 0;
    }

    return false;
}

#define TIMER_SLOT(j, N) ((j >> (TVR_BITS + (N) * TVN_BITS)) & TVN_MASK)

void WheelTimer::addTimerNode(TimerNode* node)
{
    node->slot   = INVALID_NODE_SLOTID;
    auto expires = node->expire;
    const uint64_t idx = (uint64_t)(expires - jiffies_);
    TimerList* list;
    if (idx < TVR_SIZE) // [0, 0x100)
    {
        int i = expires & TVR_MASK;
        list = &near_[i];
        node->slot = i;
    }
    else if (idx < TIME_INDEX1) // [0x100, 0x4000)
    {
        int i = TIMER_SLOT(expires, 0);
        list = &buckets_[0][i];
    }
#if WHEEL_BUCKETS >= 2
    else if (idx < TIME_INDEX2) // [0x4000, 0x100000)
    {
        int i = TIMER_SLOT(expires, 1);
        list = &buckets_[1][i];
    }
#endif
#if WHEEL_BUCKETS >= 3
    else if (idx < TIME_INDEX3) // [0x100000, 0x4000000)
    {
        int i = TIMER_SLOT(expires, 2);
        list = &buckets_[2][i];
    }
#endif
    else if ((int64_t)idx < 0)
    {
        // Can happen if you add a timer with expires == jiffies,
        // or you set a timer to go off in the past
        node->expire = jiffies_;
        int i = jiffies_ & TVR_MASK;
        list = &near_[i];
        node->slot = i;
    }
    else
    {
        // If the timeout is larger than MAX_TVAL on 64-bit
        // architectures then we use the maximum timeout
        if (idx > MAX_TVAL)
        {
             node->expire = expires = MAX_TVAL + jiffies_;
        }
        constexpr int slots = WHEEL_BUCKETS - 1;
        int i = TIMER_SLOT(expires, slots);
        list = &buckets_[slots][i];
    }
    // add to linked list
    list->emplace_back(node);
}

int64_t WheelTimer::Schedule(int64_t deadline_ms, QTask cb)
{
    auto node    = allocNode();
    node->cb     = cb;
    node->expire = deadline_ms / TIME_UNIT;
    node->id     = nextCounter();
    addTimerNode(node);
    ref_.insert_unique(node->id, node);
    //ref_[node->id] = node;
    return node->id;
}

// cascade all the timers at bucket of index up one level
bool WheelTimer::cascade(int bucket)
{
    if (bucket >= WHEEL_BUCKETS)
        return false;

    const int index = TIMER_SLOT(jiffies_, bucket);
    TimerList list = std::move(buckets_[bucket][index]);
    for (auto node : list)
    {
        if (node->expire > 0)
            addTimerNode(node);
        else {
            ref_.erase(node->id);
            freeNode(node);
        }
    }
    return index == 0;
}

// cascades all vectors and executes all expired timer
int WheelTimer::tick()
{
    const auto index = jiffies_ & TVR_MASK;
    if (index == 0 && cascade(0) && cascade(1)) // cascade timers
    {
#if WHEEL_BUCKETS > 2
        if (cascade(2) && cascade(3));
#endif
    }

    int fired = execute();
    jiffies_++;
    return fired;
}

int WheelTimer::execute()
{
    const int index = jiffies_ & TVR_MASK;
    auto& near   = near_[index];
    auto actives = near.size();
    if (actives == 0)
        return 0;

    int fired = 0, moves = 0;
    TimerList expired = std::move(near);

    for (auto node : expired)
    {
        if (node->expire > jiffies_) {
            addTimerNode(node); moves ++;
            //assert(ref_.count(node->id) == 1);
        }
        else {
            //assert(node->expire == jiffies_);
            ref_.erase(node->id);
            if (node->expire > 0) {
                //erase it then call Run, it can be removed by Cancel in Run Callback
                node->cb->Run(); fired++;
                //node->expire = 0;
            }
            freeNode(node);
        }
    }

#ifdef NDEBUG0
    if (moves > 2)
    printf("index = %d, list = %zd, timers = %u, fired/moves = %d/%d, free = %zd\n",
        index, expired.size(), ref_.size(), fired, moves, free_list_.size());
#endif

    return fired;
}

int WheelTimer::Update(int64_t now)
{
    now /= TIME_UNIT;
    const auto ticks = (int)(now - current_);
    if (ticks <= 0)
        return 0;

    current_ = now;
    int fired = 0;
    for (int i = 0; i < ticks; i++)
        fired += tick();
    return fired;
}

int64_t WheelTimer::nextTimer(int slots)
{
    int64_t next_expire = jiffies_ + TVR_MASK;
    if (ref_.bucket_count() < 50) {
        for (const auto& v : ref_) {
            auto expire = v.second->expire;
            if (expire < next_expire && expire > 0)
                next_expire = expire;
        }
    }
    else {
        for (int i = 1; i < TVR_MASK && i < slots; i++) {
            const int slots = (jiffies_ + i) & TVR_MASK;
            if (near_[slots].size() == 0)
                continue;
#ifndef NDEBUG
            //assert(near_[slots][0]->expire - jiffies_ >= i);
#endif
            for (auto node : near_[slots]) {
                if (node->expire > 0)
                    return jiffies_ + i;
            }
        }
    }

    return next_expire;
}
