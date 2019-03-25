// Copyright 2011 Google Inc. All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#pragma once

#include <atomic>
#include <iostream>
#include <iterator>
#include <stddef.h>

#include <seqan3/core/platform.hpp>

namespace seqan3::contrib
{

template <typename Queue>
class queue_back_iter
:
    public std::iterator<std::output_iterator_tag, void, void, void, void>
{
  public:
    typedef typename Queue::value_type value_type;

    queue_back_iter(Queue& q) : q_(&q) { }
    queue_back_iter() : q_(static_cast<Queue*>(NULL)) { }

    queue_back_iter& operator *() { return *this; }
    queue_back_iter& operator ++() { return *this; }
    queue_back_iter& operator ++(int) { return *this; }
    queue_back_iter& operator =(const value_type& value);

    bool operator ==(const queue_back_iter& y) { return q_ == y.q_; }
    bool operator !=(const queue_back_iter& y) { return q_ != y.q_; }

  private:
    Queue* q_;
};

template <typename Queue>
class queue_front_iter
:
    public std::iterator<std::input_iterator_tag, void, void, void, void>
{
  public:
    typedef typename Queue::value_type value_type;

    class value
    {
      public:
        value(value_type v) : v_(v) { }
        value_type operator *() const { return v_; }
      private:
        value_type v_;
    };

    queue_front_iter(Queue& q) : q_(&q) { if ( q_ ) next(); }
    queue_front_iter() : q_(static_cast<Queue*>(NULL)) { }

    const value_type& operator *() const { return v_; }
    const value_type* operator ->() const { return &v_; }
    queue_front_iter& operator ++() { next(); return *this; }
    value operator ++(int) { value t = v_; next(); return t; }

    bool operator ==(const queue_front_iter& y)
    { return q_ == y.q_; }
    bool operator !=(const queue_front_iter& y)
    { return q_ != y.q_; }

  private:
    void next();

    Queue* q_;
    value_type v_;
};

enum class queue_op_status
{
    success = 0,
    empty,
    full,
    closed,
    busy
};

#if 0
template <typename Value>
class queue_common
{
  public:
    typedef Value& reference;
    typedef const Value& const_reference;
    typedef Value value_type;

    virtual void close() = 0;
    virtual bool is_closed() = 0;
    virtual bool is_empty() = 0;

  protected:
    virtual ~queue_common();
};

template <typename Value>
queue_common<Value>::~queue_common() = default;
#endif

template <typename Queue>
class generic_queue_back
{
  public:
    typedef typename Queue::value_type value_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;

    typedef queue_back_iter<generic_queue_back> iterator;
    typedef const queue_back_iter<generic_queue_back> const_iterator;

    //FIX generic_queue_back() = default;
    generic_queue_back(Queue& queue) : queue_(&queue) { }
    generic_queue_back(Queue* queue) : queue_(queue) { }
    generic_queue_back(const generic_queue_back& other)
        = default;
    generic_queue_back& operator =(const generic_queue_back& other)
        = default;

    void close() { queue_->close(); }
    bool is_closed() { return queue_->is_closed(); }
    bool is_empty() { return queue_->is_empty(); }

    iterator begin() { return iterator(*this); }
    iterator end() { return iterator(); }
    const iterator cbegin() { return const_iterator(*this); }
    const iterator cend() { return const_iterator(); }

    void push(const value_type& x)
        { queue_->push(x); }
    queue_op_status wait_push(const value_type& x)
        { return queue_->wait_push(x); }
    queue_op_status try_push(const value_type& x)
        { return queue_->try_push(x); }
    queue_op_status nonblocking_push(const value_type& x)
        { return queue_->nonblocking_push(x); }
    void push(value_type&& x)
        { queue_->push( std::move(x) ); }
    queue_op_status wait_push(value_type&& x)
        { return queue_->wait_push( std::move(x) ); }
    queue_op_status try_push(value_type&& x)
        { return queue_->try_push( std::move(x) ); }
    queue_op_status nonblocking_push(value_type&& x)
        { return queue_->nonblocking_push( std::move(x) ); }

    bool has_queue() { return queue_ != NULL; }

  protected:
    Queue* queue_;
};

template <typename Queue>
class generic_queue_front
{
  public:
    typedef typename Queue::value_type value_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;

    typedef queue_front_iter<generic_queue_front> iterator;
    typedef queue_front_iter<generic_queue_front> const_iterator;

    //FIX generic_queue_front() = default;
    generic_queue_front(Queue& queue) : queue_(&queue) { }
    generic_queue_front(Queue* queue) : queue_(queue) { }
    generic_queue_front(const generic_queue_front& other)
        = default;
    generic_queue_front& operator =(const generic_queue_front& other)
        = default;

    void close() { queue_->close(); }
    bool is_closed() { return queue_->is_closed(); }
    bool is_empty() { return queue_->is_empty(); }

    iterator begin() { return iterator(*this); }
    iterator end() { return iterator(); }
    const iterator cbegin() { return const_iterator(*this); }
    const iterator cend() { return const_iterator(); }

    value_type value_pop()
        { return queue_->value_pop(); }
    queue_op_status wait_pop(value_type& x)
        { return queue_->wait_pop(x); }
    queue_op_status try_pop(value_type& x)
        { return queue_->try_pop(x); }
    queue_op_status nonblocking_pop(value_type& x)
        { return queue_->nonblocking_pop(x); }

    bool has_queue() { return queue_ != NULL; }

  protected:
    Queue* queue_;
};

template <typename Queue>
queue_back_iter<Queue>&
queue_back_iter<Queue>::operator =(const value_type& value)
{
    queue_op_status s = q_->wait_push(value);
    if ( s != queue_op_status::success ) {
        q_ = NULL;
        throw s;
    }
    return *this;
}

template <typename Queue>
void
queue_front_iter<Queue>::next()
{
    queue_op_status s = q_->wait_pop(v_);
    if ( s == queue_op_status::closed )
        q_ = NULL;
}

template <typename Value>
class queue_base
{
  public:
    typedef Value value_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;

    virtual ~queue_base() { }

    virtual void close() = 0;
    virtual bool is_closed() = 0;
    virtual bool is_empty() = 0;

    virtual void push(const Value& x) = 0;
    virtual queue_op_status wait_push(const Value& x) = 0;
    virtual queue_op_status try_push(const Value& x) = 0;
    virtual queue_op_status nonblocking_push(const Value& x) = 0;

    virtual void push(Value&& x) = 0;
    virtual queue_op_status wait_push(Value&& x) = 0;
    virtual queue_op_status try_push(Value&& x) = 0;
    virtual queue_op_status nonblocking_push(Value&& x) = 0;

    virtual Value value_pop() = 0;
    virtual queue_op_status wait_pop(Value&) = 0;
    virtual queue_op_status try_pop(Value&) = 0;
    virtual queue_op_status nonblocking_pop(Value&) = 0;
};

//TODO(crowl): Use template aliases for queue_back and queue_front?

template <typename Value>
class queue_back
: public generic_queue_back< queue_base<Value> >
{
  public:
    queue_back() = default;
    queue_back(queue_base<Value>& queue)
        : generic_queue_back< queue_base<Value> >(queue) { }
    queue_back(queue_base<Value>* queue)
        : generic_queue_back< queue_base<Value> >(queue) { }
    queue_back(const queue_back<Value>& other)
        : generic_queue_back< queue_base<Value> >(other.queue_) { }
};

template <typename Value>
class queue_front
: public generic_queue_front< queue_base<Value> >
{
  public:
    queue_front() = default;
    queue_front(queue_base<Value>& queue)
        : generic_queue_front< queue_base<Value> >(queue) { }
    queue_front(queue_base<Value>* queue)
        : generic_queue_front< queue_base<Value> >(queue) { }
    queue_front(const queue_front<Value>& other)
        : generic_queue_front< queue_base<Value> >(other.queue_) { }
};

template <typename Queue>
class queue_wrapper
:
    public virtual queue_base <typename Queue::value_type>
{
    Queue* ptr;

  public:
    typedef typename Queue::value_type value_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;

    queue_wrapper(Queue* arg)
    : ptr(arg)
    { }

    queue_wrapper(Queue& arg)
    : ptr(&arg)
    { }

    virtual ~queue_wrapper()
    { }

    virtual void close()
    { ptr->close(); }

    virtual bool is_closed()
    { return ptr->is_closed(); }

    virtual bool is_empty()
    { return ptr->is_empty(); }

    virtual void push(const value_type& x)
    { ptr->push(x); }

    virtual queue_op_status wait_push(const value_type& x)
    { return ptr->wait_push(x); }

    virtual queue_op_status try_push(const value_type& x)
    { return ptr->try_push(x); }

    virtual queue_op_status nonblocking_push(const value_type& x)
    { return ptr->nonblocking_push(x); }

    virtual void push(value_type&& x)
    { ptr->push(std::move(x)); }

    virtual queue_op_status wait_push(value_type&& x)
    { return ptr->wait_push(std::move(x)); }

    virtual queue_op_status try_push(value_type&& x)
    { return ptr->try_push(std::move(x)); }

    virtual queue_op_status nonblocking_push(value_type&& x)
    { return ptr->nonblocking_push(std::move(x)); }

    virtual value_type value_pop()
    { return ptr->value_pop(); }

    virtual queue_op_status wait_pop(value_type& x)
    { return ptr->wait_pop(x); }

    virtual queue_op_status try_pop(value_type& x)
    { return ptr->try_pop(x); }

    virtual queue_op_status nonblocking_pop(value_type& x)
    { return ptr->nonblocking_pop(x); }

    queue_back<value_type> back()
    { return queue_back<value_type>(this); }

    queue_front<value_type> front()
    { return queue_front<value_type>(this); }
};

template <typename Value>
class queue_counted
:
    public queue_base<Value>
{
  public:
    queue_counted() : bk_(0), ft_(0) { }
    virtual ~queue_counted() { }

    void inc_back() { bk_++; }
    void inc_front() { ft_++; }
    bool dec_back() { return --bk_ == 0; }
    bool dec_front() { return --ft_ == 0; }
    bool no_back() { return bk_ == 0; }
    bool no_front() { return ft_ == 0; }

  private:
    std::atomic<int> bk_;
    std::atomic<int> ft_;
};

template <typename Value>
class shared_queue_back
{
  public:
    typedef Value value_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;

    typedef queue_back_iter<shared_queue_back> iterator;
    typedef const queue_back_iter<shared_queue_back> const_iterator;

    //FIX shared_queue_back()
    //FIX     : queue_(NULL) { }
    shared_queue_back(queue_counted<value_type>* queue)
        : queue_(queue) { queue->inc_back(); }
    shared_queue_back(const shared_queue_back& other)
        : queue_(other.queue_) { queue_->inc_back(); }
    shared_queue_back(shared_queue_back&& other)
        : queue_(other.queue_) { other.queue_ = NULL; }

  private:
    void release()
    {
        if ( queue_ != NULL && queue_->dec_back() ) {
            queue_->close();
            if ( queue_->no_front() ) {
                delete queue_;
            }
        }
    }

  public:
    ~shared_queue_back() { release(); }

    shared_queue_back& operator =(const shared_queue_back& other)
    {
        if ( this != &other ) {
            release();
            queue_ = other->queue_;
            if ( queue_ != NULL )
                queue_->inc_back();
        }
        return *this;
    }
    shared_queue_back& operator =(shared_queue_back&& other)
    {
        if ( this != &other ) {
            release();
            queue_ = other->queue_;
            other->queue_ == NULL;
        }
        return *this;
    }

    void close() { queue_->close(); }
    bool is_closed() { return queue_->is_closed(); }
    bool is_empty() { return queue_->is_empty(); }

    iterator begin() { return iterator(*this); }
    iterator end() { return iterator(); }
    const iterator cbegin() { return const_iterator(*this); }
    const iterator cend() { return const_iterator(); }

    void push(const value_type& x)
        { queue_->push(x); }
    queue_op_status wait_push(const value_type& x)
        { return queue_->wait_push(x); }
    queue_op_status try_push(const value_type& x)
        { return queue_->try_push(x); }
    queue_op_status nonblocking_push(const value_type& x)
        { return queue_->nonblocking_push(x); }
    void push(value_type&& x)
        { queue_->push( std::move(x) ); }
    queue_op_status wait_push(value_type&& x)
        { return queue_->wait_push( std::move(x) ); }
    queue_op_status try_push(value_type&& x)
        { return queue_->try_push( std::move(x) ); }
    queue_op_status nonblocking_push(value_type&& x)
        { return queue_->nonblocking_push( std::move(x) ); }

  private:
    queue_counted<value_type>* queue_;
};

template <typename Value>
class shared_queue_front
{
  public:
    typedef Value value_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;

    typedef queue_front_iter<shared_queue_front> iterator;
    typedef queue_front_iter<shared_queue_front> const_iterator;

    //FIX shared_queue_front()
    //FIX     : queue_(NULL) { }
    shared_queue_front(queue_counted<value_type>* queue)
        : queue_(queue) { queue->inc_front(); }
    shared_queue_front(const shared_queue_front& other)
        : queue_(other.queue_) { queue_->inc_front(); }
    shared_queue_front(shared_queue_front&& other)
        : queue_(other.queue_) { other.queue_ = NULL; }

  private:
    void release()
    {
        if ( queue_ != NULL && queue_->dec_front() ) {
            queue_->close();
            if ( queue_->no_back() ) {
                delete queue_;
            }
        }
    }

  public:
    ~shared_queue_front() { release(); }

    shared_queue_front& operator =(const shared_queue_front& other)
    {
        if ( this != &other ) {
            release();
            queue_ = other->queue_;
            if ( queue_ != NULL )
                queue_->inc_front();
        }
        return *this;
    }

    shared_queue_front& operator =(shared_queue_front&& other)
    {
        if ( this != &other ) {
            release();
            queue_ = other->queue_;
            other->queue_ = NULL;
        }
        return *this;
    }

    void close() { queue_->close(); }
    bool is_closed() { return queue_->is_closed(); }
    bool is_empty() { return queue_->is_empty(); }

    iterator begin() { return iterator(*this); }
    iterator end() { return iterator(); }
    const iterator cbegin() { return const_iterator(*this); }
    const iterator cend() { return const_iterator(); }

    value_type value_pop()
        { return queue_->value_pop(); }
    queue_op_status wait_pop(value_type& x)
        { return queue_->wait_pop(x); }
    queue_op_status try_pop(value_type& x)
        { return queue_->try_pop(x); }
    queue_op_status nonblocking_pop(value_type& x)
        { return queue_->nonblocking_pop(x); }

  private:
    queue_counted<value_type>* queue_;
};

template <typename Queue>
class queue_owner
:
    public queue_counted <typename Queue::value_type>
{
    Queue* ptr;

  public:
    typedef typename Queue::value_type value_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;

    queue_owner(const queue_owner&) = delete;
    queue_owner(Queue* arg) : ptr(arg) { }

    virtual ~queue_owner() { delete ptr; }

    queue_back<value_type> back()
        { return queue_back<value_type>(this); }
    queue_front<value_type> front()
        { return queue_front<value_type>(this); }

    virtual void close() { ptr->close(); }
    virtual bool is_closed() { return ptr->is_closed(); }
    virtual bool is_empty() { return ptr->is_empty(); }

    virtual void push(const value_type& x)
        { ptr->push(x); }
    virtual queue_op_status wait_push(const value_type& x)
        { return ptr->wait_push(x); }
    virtual queue_op_status try_push(const value_type& x)
        { return ptr->try_push(x); }
    virtual queue_op_status nonblocking_push(const value_type& x)
        { return ptr->nonblocking_push(x); }

    virtual void push(value_type&& x)
        { ptr->push(std::move(x)); }
    virtual queue_op_status wait_push(value_type&& x)
        { return ptr->wait_push(std::move(x)); }
    virtual queue_op_status try_push(value_type&& x)
        { return ptr->try_push(std::move(x)); }
    virtual queue_op_status nonblocking_push(value_type&& x)
        { return ptr->nonblocking_push(std::move(x)); }

    virtual value_type value_pop()
        { return ptr->value_pop(); }
    virtual queue_op_status wait_pop(value_type& x)
        { return ptr->wait_pop(x); }
    virtual queue_op_status try_pop(value_type& x)
        { return ptr->try_pop(x); }
    virtual queue_op_status nonblocking_pop(value_type& x)
        { return ptr->nonblocking_pop(x); }
};


template <typename Queue>
class queue_object
:
    public queue_counted <typename Queue::value_type>
{
    Queue obj_;

  public:
    typedef typename Queue::value_type value_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;

    queue_object(const queue_object&) = delete;
    template <typename ... Args>
    queue_object(Args ... args) : obj_(args...) { }

    virtual ~queue_object() { }

    operator queue_back<value_type>() //TODO(crowl) Really?
        { return queue_back<value_type>(this); }
    queue_back<value_type> back()
        { return queue_back<value_type>(this); }
    queue_front<value_type> front()
        { return queue_front<value_type>(this); }

    virtual void close() { obj_.close(); }
    virtual bool is_closed() { return obj_.is_closed(); }
    virtual bool is_empty() { return obj_.is_empty(); }

    virtual void push(const value_type& x)
        { obj_.push(x); }
    virtual queue_op_status wait_push(const value_type& x)
        { return obj_.wait_push(x); }
    virtual queue_op_status try_push(const value_type& x)
        { return obj_.try_push(x); }
    virtual queue_op_status nonblocking_push(const value_type& x)
        { return obj_.nonblocking_push(x); }

    virtual void push(value_type&& x)
        { obj_.push(std::move(x)); }
    virtual queue_op_status wait_push(value_type&& x)
        { return obj_.wait_push(std::move(x)); }
    virtual queue_op_status try_push(value_type&& x)
        { return obj_.try_push(std::move(x)); }
    virtual queue_op_status nonblocking_push(value_type&& x)
        { return obj_.nonblocking_push(std::move(x)); }

    virtual value_type value_pop()
        { return obj_.value_pop(); }
    virtual queue_op_status wait_pop(value_type& x)
        { return obj_.wait_pop(x); }
    virtual queue_op_status try_pop(value_type& x)
        { return obj_.try_pop(x); }
    virtual queue_op_status nonblocking_pop(value_type& x)
        { return obj_.nonblocking_pop(x); }
};

template <typename Queue, typename ... Args>
std::pair< shared_queue_back<typename Queue::value_type>,
           shared_queue_front<typename Queue::value_type> >
share_queue_ends(Args ... args)
{
  typedef typename Queue::value_type elemtype;
  auto  q =  new queue_object<Queue>(args...) ;
  return std::make_pair(shared_queue_back<elemtype>(q),
                        shared_queue_front<elemtype>(q));
}

} // namespace seqan3::contrib
