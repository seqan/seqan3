// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================
// Authors: Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
//          Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
// ============================================================================
// Implementation of the (const) iterator given a container reference.
// ============================================================================

#pragma once

#include <seqan3/container/concepts.hpp>

#include <iomanip>
#include <limits>
#include <type_traits>

namespace seqan3::detail
{

template <typename container_type, bool is_const=false>
    requires container_concept<container_type> //requires random_access_range_concept<container_type>
struct ra_iterator
{

private:
    container_type & host{*(container_type*)0};
    typename container_type::size_type pos{0};
public:
    typedef typename container_type::size_type size_type;
    typedef typename container_type::difference_type difference_type;
    typedef typename container_type::value_type value_type;
    using reference = std::conditional_t<is_const, value_type const &, value_type &>;
    typedef typename container_type::const_reference const_reference;
    using pointer = std::conditional_t<is_const, value_type const *, value_type *>;

    typedef typename container_type::iterator iterator;
    typedef typename container_type::const_iterator const_iterator;

    //! Note: iterator can only be initialized with container reference
    //ra_iterator() : host(*(container_type*)0) {};
    //! set position index to max in order to indicate that iterator exceeds container size
    ra_iterator(container_type & _host, bool const _at_end = false) : host{_host}, pos{0}
    {
        if (_at_end)
            pos = std::numeric_limits<size_type>::max();
    }

    //! copy constructor
    ra_iterator(ra_iterator const & _it) : host{_it.host}, pos{_it.pos} {}

    //! copy via assignment
    ra_iterator& operator=(ra_iterator const & rhs)
    {
        //assert(host == rhs.host);
        // TODO: remove above assert, because with assert: this is not possible:
        // it = it + 2, but we defined it += 2; Hence, w would expect same behaviour.
        pos = rhs.pos;
        return *this;
    }

    //! move operator, TODO: needed for this light-weight class?
    /*ra_iterator& operator=(ra_iterator && rhs)
    {
        //assert(host == rhs.host);
        swap(*this, rhs);
        return *this;
    }*/

    ~ra_iterator() = default;

    //! return element currently pointed at
    reference operator*()
    {
        return host[pos];
    }

    //! two iterators are equal if their absolute positions are the same
    bool operator==(ra_iterator const & rhs) const
    {
        return this->pos == rhs.pos;
    }

    bool operator!=(ra_iterator const & rhs) const
    {
        return this->pos != rhs.pos;
    }

    bool operator<(ra_iterator const & rhs) const
    {
        return this->pos < rhs.pos;
    }

    bool operator>(ra_iterator const & rhs) const
    {
        return this->pos > rhs.pos;
    }

    bool operator<=(ra_iterator const & rhs) const
    {
        return this->pos <= rhs.pos;
    }

    bool operator>=(ra_iterator const & rhs) const
    {
        return this->pos >= rhs.pos;
    }

    //! prefix increment
    ra_iterator& operator++()
    {
        ++pos;
        return (*this);
    }

    //! postfix increment
    ra_iterator operator++(int)
    {
        ra_iterator cpy{*this};
        ++pos;
        if (pos >= host.size())
            pos = std::numeric_limits<size_type>::max();
        return cpy;
    }

    //! prefix decrement
    ra_iterator& operator--()
    {
        --pos;
        return *this;
    }

    //! postfix decrement
    ra_iterator operator--(int)
    {
        ra_iterator cpy{*this};
        --pos;
        return cpy;
    }

    //! forward iterator and assignment
    ra_iterator& operator+=(size_type skip)
    {
        pos += skip;
        return *this;
    }

    //! forward iterator    ra_iterator operator+(size_type) const; //optional
    ra_iterator operator+(size_type skip) const
    {
        ra_iterator cpy{*this};
        cpy.pos = cpy.pos + skip;
        return cpy;
    }

    friend ra_iterator operator+(size_type skip , const ra_iterator& _it)
    {
        ra_iterator cpy{_it};
        cpy.pos = cpy.pos + skip;
        return cpy;
    }

    ra_iterator& operator-=(size_type skip)
    {
        pos -= skip;
        return *this;
    }

    ra_iterator operator-(size_type skip) const
    {
        ra_iterator cpy{*this};
        cpy.pos -= skip;
        return cpy;
    }

    difference_type operator-(ra_iterator lhs) const
    {
        return pos - lhs.pos;
    }

    pointer operator->()
    {
        return &(*this);
    }

    reference operator[](size_type const n)
    {
        return host[pos + n];
    }

};

} // namespace seqan3::detail
