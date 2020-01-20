// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan suspendable queue.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \author David Weese <david.weese AT fu-berlin.de>
 */

#pragma once

#include <condition_variable>
#include <mutex>
#include <thread>
#include <vector>

#include <seqan3/core/platform.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>
#include <seqan3/std/span>

namespace seqan3::contrib
{

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class ConcurrentQueue
// ----------------------------------------------------------------------------

template <typename TSpec = void>
struct Tag{};

template <typename TSpec = Tag<void>>
struct Suspendable;

template <typename TValue, typename TSpec = Suspendable<>>
class ConcurrentQueue;

struct Limit_;
using Limit = Tag<Limit_>;

template <typename TValue, typename TSpec>
class ConcurrentQueue<TValue, Suspendable<TSpec> >
{
public:
    typedef std::vector<TValue>          TString;
    typedef typename TString::size_type  TSize;

    size_t                  readerCount;
    size_t                  writerCount;

    TString                 data;
    TSize                   occupied;
    TSize                   back;
    TSize                   front;

    std::mutex              cs;
    std::condition_variable more;

    bool                    virgin;

    ConcurrentQueue():
        readerCount(0),
        writerCount(0),
        occupied(0),
        back(0),
        front(0),
        virgin(true)
    {}

    ~ConcurrentQueue()
    {
        assert(writerCount == 0u);

        // wait for all pending readers to finish
        while (readerCount != 0u)
        {}
    }
};

template <typename TValue>
class ConcurrentQueue<TValue, Suspendable<Limit> >:
    public ConcurrentQueue<TValue, Suspendable<> >
{
public:
    typedef ConcurrentQueue<TValue, Suspendable<> > TBase;
    typedef typename TBase::TString                 TString;
    typedef typename TBase::TSize                   TSize;

    std::condition_variable less;

    ConcurrentQueue(TSize maxSize):
        TBase()
    {
        this->data.resize(maxSize);
        // reserve(this->data, maxSize, Exact());
        // _setLength(this->data, maxSize);
    }

    ConcurrentQueue(ConcurrentQueue const & other):
        TBase((TBase const &)other)
    {}
};

template <typename TValue, typename TSpec>
inline void
lockReading(ConcurrentQueue<TValue, Suspendable<TSpec> > &)
{}

template <typename TValue, typename TSpec>
inline void
unlockReading(ConcurrentQueue<TValue, Suspendable<TSpec> > & me)
{
    {
        std::lock_guard<std::mutex> lock(me.cs);
        if (--me.readerCount != 0u)
            return;
    }
    me.less.notify_all();  // publish the condition that reader count is 0.
}

template <typename TValue, typename TSpec>
inline void
lockWriting(ConcurrentQueue<TValue, Suspendable<TSpec> > &)
{}

template <typename TValue, typename TSpec>
inline void
unlockWriting(ConcurrentQueue<TValue, Suspendable<TSpec> > & me)
{
    {
        std::lock_guard<std::mutex> lk(me.cs);
        if (--me.writerCount != 0u)
            return;
    }
    me.more.notify_all();  // publish the condition, that writer count is 0.
}

template <typename TValue, typename TSize, typename TSpec>
inline void
setReaderCount(ConcurrentQueue<TValue, Suspendable<TSpec> > & me, TSize readerCount)
{
    std::unique_lock<std::mutex> lock(me.cs);
    me.readerCount = readerCount;
}

template <typename TValue, typename TSize, typename TSpec>
inline void
setWriterCount(ConcurrentQueue<TValue, Suspendable<TSpec> > & me, TSize writerCount)
{
    std::unique_lock<std::mutex> lock(me.cs);
    me.writerCount = writerCount;
}

template <typename TValue, typename TSize1, typename TSize2, typename TSpec>
inline void
setReaderWriterCount(ConcurrentQueue<TValue, Suspendable<TSpec> > & me, TSize1 readerCount, TSize2 writerCount)
{
    std::unique_lock<std::mutex> lock(me.cs);
    me.readerCount = readerCount;
    me.writerCount = writerCount;
}

template <typename TValue, typename TSize, typename TSpec>
inline bool
waitForMinSize(ConcurrentQueue<TValue, Suspendable<TSpec> > & me,
               TSize minSize)
{
    std::unique_lock<std::mutex> lock(me.cs);
    while (me.occupied < minSize && me.writerCount > 0u)
        me.more.wait(lock);
    return me.occupied >= minSize;
}

template <typename TValue, typename TSpec>
inline bool
empty(ConcurrentQueue<TValue, Suspendable<TSpec> > const & me)
{
    return me.occupied == 0;
}

template <typename TValue, typename TSpec>
inline typename ConcurrentQueue<TValue, Suspendable<TSpec> >::SizeType
length(ConcurrentQueue<TValue, Suspendable<TSpec> > const & me)
{
    return me.occupied;
}

template <typename TValue, typename TSpec>
inline bool
_popFront(TValue & result, ConcurrentQueue<TValue, Suspendable<TSpec> > & me,
          std::unique_lock<std::mutex> & lk)
{
    typedef ConcurrentQueue<TValue, Suspendable<TSpec> >    TQueue;
    typedef typename TQueue::TString                        TString;
    typedef typename TString::size_type                     TSize;

    TSize cap = me.data.size();

    while (me.occupied == 0u && me.writerCount > 0u)
        me.more.wait(lk);

    if (me.occupied == 0u)
        return false;

    assert(me.occupied > 0u);

    // extract value and destruct it in the data string
    // TIter it = me.data.begin() + me.front;
    result = std::ranges::iter_move(std::ranges::next(me.data.begin(), me.front));
    // std::swap(result, *it);
    // valueDestruct(it);

    me.front = (me.front + 1) % cap;
    me.occupied--;

    /* now: either me.occupied > 0 and me.nextout is the index
       of the next occupied slot in the buffer, or
       me.occupied == 0 and me.nextout is the index of the next
       (empty) slot that will be filled by a producer (such as
       me.nextout == me.nextin) */

    return true;
}

template <typename TValue, typename TSpec>
inline bool
_popBack(TValue & result,
         ConcurrentQueue<TValue, Suspendable<TSpec> > & me,
         std::unique_lock<std::mutex> & lk)
{
    typedef ConcurrentQueue<TValue, Suspendable<TSpec> > TQueue;
    typedef typename TQueue::TString                     TString;
    typedef typename TString::size_type                  TSize;

    TSize cap = me.data.size();

    while (me.occupied == 0u && me.writerCount > 0u)
        me.more.wait(lk);

    if (me.occupied == 0u)
        return false;

    assert(me.occupied > 0u);

    me.back = (me.back + cap - 1) % cap;

    // extract value and destruct it in the data string
    // TIter it = me.data.begin() + me.back;
    result = std::ranges::iter_move(std::ranges::next(me.data.begin(), me.back));
    // std::swap(result, *it);
    // valueDestruct(it);

    me.occupied--;

    /* now: either me.occupied > 0 and me.nextout is the index
       of the next occupied slot in the buffer, or
       me.occupied == 0 and me.nextout is the index of the next
       (empty) slot that will be filled by a producer (such as
       me.nextout == me.nextin) */

    return true;
}

template <typename TValue, typename TSpec>
inline bool
popFront(TValue & result, ConcurrentQueue<TValue, Suspendable<TSpec> > & me)
{
    std::unique_lock<std::mutex> lock(me.cs);
    return _popFront(result, me, lock);
}

template <typename TValue>
inline bool
popFront(TValue & result, ConcurrentQueue<TValue, Suspendable<Limit> > & me)
{
    {
        std::unique_lock<std::mutex> lk(me.cs);
        if (!_popFront(result, me, lk))
            return false;
    }
    me.less.notify_all();
    return true;
}

template <typename TValue, typename TSpec>
inline bool
popBack(TValue & result, ConcurrentQueue<TValue, Suspendable<TSpec> > & me)
{
    std::unique_lock<std::mutex> lk(me.cs);
    return _popBack(result, me, lk);
}

template <typename TValue>
inline bool
popBack(TValue & result, ConcurrentQueue<TValue, Suspendable<Limit> > & me)
{
    {
        std::unique_lock<std::mutex> lk(me.cs);
        if (!_popBack(result, me, lk))
            return false;
    }
    me.less.notify_all();
    return true;
}


template <typename TValue, typename TValue2, typename TSpec, typename TExpand>
inline bool
appendValue(ConcurrentQueue<TValue, Suspendable<TSpec> > & me,
            TValue2 && val,
            [[maybe_unused]] Tag<TExpand> expandTag)
{
    typedef ConcurrentQueue<TValue, Suspendable<TSpec> > TQueue;
    typedef typename TQueue::TString                     TString;
    typedef typename TString::size_type                  TSize;

    {
        std::lock_guard<std::mutex> lock(me.cs);
        TSize cap = me.data.size();

        if (me.occupied >= cap)
        {
            // increase capacity
            // _setLength(me.data, cap);
            // reserve(me.data, cap + 1, expandTag);
            me.data.resize(cap + 1);
            TSize delta = me.data.size() - cap;
            assert(delta == 1);

            // create a gap of delta many values between tail and head
            // Why?
            // _clearSpace(me.data, delta, me.back, me.back, expandTag);
            std::ranges::move_backward(std::span{me.data.data() + me.front, me.data.data() + cap},
                                       me.data.data() + me.data.size());
            if (me.occupied != 0 && me.back <= me.front)
                me.front += delta;

            cap += delta;
        }

        // valueConstruct(begin(me.data, Standard()) + me.back, val);
        *std::ranges::next(me.data.begin(), me.back) = std::forward<TValue2>(val);
        me.back = (me.back + 1) % cap;

        ++me.occupied;
    }

    /* now: either me.occupied < BSIZE and me.nextin is the index
     of the next empty slot in the buffer, or
     me.occupied == BSIZE and me.nextin is the index of the
     next (occupied) slot that will be emptied by a consumer
     (such as me.nextin == me.nextout) */
    me.more.notify_all();
    return true;
}

template <typename TValue, typename TValue2, typename TSpec, typename TExpand>
inline bool
appendValue(ConcurrentQueue<TValue, Suspendable<Limit> > & me,
            TValue2 && val,
            Tag<TExpand> expandTag);

template <typename TValue, typename TValue2>
inline bool
appendValue(ConcurrentQueue<TValue, Suspendable<Limit> > & me,
            TValue2 && val,
            Limit)
{
    typedef ConcurrentQueue<TValue, Suspendable<Limit> > TQueue;
    typedef typename TQueue::TString                     TString;
    typedef typename TString::size_type                  TSize;

    {
        std::unique_lock<std::mutex> lock(me.cs);
        TSize cap = me.data.size();

        while (me.occupied >= cap && me.readerCount > 0u)
            me.less.wait(lock);

        if (me.occupied >= cap)
            return false;

        assert(me.occupied < cap);

        // valueConstruct(begin(me.data, Standard()) + me.back, val);
        *std::ranges::next(me.data.begin(), me.back) = std::forward<TValue2>(val);
        me.back = (me.back + 1) % cap;
        me.occupied++;
    }

    /* now: either me.occupied < BSIZE and me.nextin is the index
     of the next empty slot in the buffer, or
     me.occupied == BSIZE and me.nextin is the index of the
     next (occupied) slot that will be emptied by a consumer
     (such as me.nextin == me.nextout) */
    me.more.notify_all();
    return true;
}

template <typename TValue, typename TValue2, typename TSpec>
inline bool
appendValue(ConcurrentQueue<TValue, Suspendable<TSpec> > & me,
            TValue2 && val)
{
    return appendValue(me, std::forward<TValue2>(val), TSpec{});
}

} // namespace seqan3::contrib
