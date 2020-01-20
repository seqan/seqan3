// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides helper structs from SeqAn2 for the bgzf_ostream.
 * \author David Weese <david.weese AT fu-berlin.de>
 * \author Rene Rahn <rene.rahn AT fu-berlin.de
 */

#pragma once

#include <mutex>

#include <seqan3/contrib/parallel/suspendable_queue.hpp>

namespace seqan3::contrib
{

// ============================================================================
// Classes
// ============================================================================

template <typename TValue>
struct ResourcePool
{
    typedef ConcurrentQueue<TValue *, Suspendable<> > TStack;
    typedef typename TStack::TSize                    TSize;

    TStack recycled;

    ResourcePool(TSize maxSize)
    {
        setWriterCount(recycled, 1);
        for (; maxSize != 0; --maxSize)
            appendValue(recycled, (TValue *)NULL);
    }

    ~ResourcePool()
    {
        unlockWriting(recycled);
        TValue *ptr = NULL;
        unsigned count = 0;
        while (popBack(ptr, recycled))
        {
            if (ptr != NULL)
                count++;
            delete ptr;
        }
    }
};

// ----------------------------------------------------------------------------
// Struct SerializerItem
// ----------------------------------------------------------------------------

template <typename TValue>
struct SerializerItem
{
    TValue          val;
    SerializerItem  *next;
    bool            ready;
};

// ----------------------------------------------------------------------------
// Class Serializer
// ----------------------------------------------------------------------------

template <typename TValue, typename TWorker>
class Serializer
{
public:
    typedef SerializerItem<TValue>          TItem;
    typedef TItem *                         TItemPtr;
    typedef ResourcePool<TItem>             TPool;
    typedef size_t                          TSize;

    std::mutex          cs;
    TWorker             worker;
    TItemPtr            first;
    TItemPtr            last;
    TPool               pool;
    bool                stop;

    Serializer() :
        first(NULL),
        last(NULL),
        stop(false)
    {}

    template <typename TArg>
    explicit
    Serializer(TArg &arg, TSize maxItems = 1024) :
        worker(arg),
        first(NULL),
        last(NULL),
        pool(maxItems),
        stop(false)
    {}

    ~Serializer()
    {
        while (first != NULL)
        {
            TItemPtr item = first;
            first = first->next;
            delete item;
        }
    }

    operator bool()
    {
        return !stop;
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function aquireValue()
// ----------------------------------------------------------------------------

template <typename TValue>
inline TValue *
aquireValue(ResourcePool<TValue> & me)
{
    TValue *ptr = NULL;
    if (!popBack(ptr, me.recycled))
        return NULL;

    if (ptr == NULL)
        ptr = new TValue;

    return ptr;
}

// ----------------------------------------------------------------------------
// Function releaseValue()
// ----------------------------------------------------------------------------

template <typename TValue>
inline void
releaseValue(ResourcePool<TValue> & me, TValue *ptr)
{
    appendValue(me.recycled, ptr);
}


// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TValue, typename TWorker>
inline void
clear(Serializer<TValue, TWorker> & me)
{
    me.stop = false;
    while (me.first != NULL)
    {
        TValue *item = me.first;
        me.first = me.first->next;
        releaseValue(me.recycled, item);
    }
    me.last = NULL;
}

// ----------------------------------------------------------------------------
// Function aquireValue()
// ----------------------------------------------------------------------------

// this function is not thread-safe as it would make
// not much sense to order a stream by the random
// order of executition behind a mutex
template <typename TValue, typename TWorker>
inline TValue *
aquireValue(Serializer<TValue, TWorker> & me)
{
    typedef SerializerItem<TValue> TItem;

    TItem *item = aquireValue(me.pool);
    item->next = NULL;
    item->ready = false;

    // add item to the end of our linked list
    {
        std::lock_guard<std::mutex> lock(me.cs);
        if (me.first == NULL)
            me.first = item;
        else
            me.last->next = item;
        me.last = item;
    }
    return &item->val;
}

// ----------------------------------------------------------------------------
// Function releaseValue()
// ----------------------------------------------------------------------------

template <typename TValue, typename TWorker>
inline bool
releaseValue(Serializer<TValue, TWorker> & me, TValue *ptr)
{
    typedef SerializerItem<TValue> TItem;

    TItem *item = reinterpret_cast<TItem *>(ptr);
    assert(!item->ready);

    // changing me.first or the ready flag must be done synchronized (me.mutex)
    // the thread who changed me.first->ready to be true has to write it.

    // change our ready flag and test if me.first->ready became true
    {
        std::lock_guard<std::mutex> lock(me.cs);
        item->ready = true;
        if (item != me.first)
            return true;
    }

    // ok, if we are here it seems that we are responsible for writing the buffer

    assert(me.first != NULL);

    bool success;
    do
    {
        // process item
        success = me.worker(item->val);

        // remove item from linked list
        {
            std::lock_guard<std::mutex> lock(me.cs);
            me.first = item->next;

            // recycle released items
            releaseValue(me.pool, item);

            // can we leave?
            item = me.first;
            if (item == NULL || !item->ready)
                return success;
        }

        // we continue to process the next buffer
    }
    while (success);

    return false;
}

}  // namespace seqan3::contrib
