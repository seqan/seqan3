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
// Authors: Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
// ============================================================================
// Implementation of the (const) iterator given a container reference.
// ============================================================================

#pragma once

#include <iomanip>
#include <limits>
#include <type_traits>

#include <seqan3/container/concepts.hpp>

namespace seqan3::detail
{

template <typename container_type, bool is_const=false>
    requires random_access_range_concept<container_type>
struct ra_iterator
{

private:
    //! Note: an iterator can only be initialized with container (_host) reference
    container_type & _host{*(container_type*)0};
    //! set to numeric_limits::max to indicate that iterator exceeds container size
    typename container_type::size_type _pos{0};
public:
    typedef typename container_type::size_type size_type;
    typedef typename container_type::difference_type difference_type;
    typedef typename container_type::value_type value_type;
    using reference = std::conditional_t<is_const, value_type const &, value_type &>;
    typedef typename container_type::const_reference const_reference;
    using pointer = std::conditional_t<is_const, value_type const *, value_type *>;

    //! Explicitly not possible: ra_iterator() : _host(*(container_type*)0) {};
    ra_iterator(container_type & _host, bool const _at_end = false) : _host{_host}, _pos{0}
    {
        if (_at_end)
            _pos = std::numeric_limits<size_type>::max();
    }

    //! copy constructor
    ra_iterator(ra_iterator const & _it) : _host{_it._host}, _pos{_it._pos} {}

    //! copy via assignment
    ra_iterator& operator=(ra_iterator const & rhs)
    {
        //assert(_host == rhs._host);
        // TODO: remove above assert, because with assert: this is not possible:
        // it = it + 2, but we allow it += 2; Hence, we would expect same behaviour.
        _pos = rhs._pos;
        return *this;
    }

    ~ra_iterator() = default;

    //! return element currently pointed at
    reference operator*() const
    {
        return _host[_pos];
    }

    //! two iterators are equal if their absolute positions are the same
    bool operator==(ra_iterator const & rhs) const
    {
        return this->_pos == rhs._pos;
    }

    bool operator!=(ra_iterator const & rhs) const
    {
        return this->_pos != rhs._pos;
    }

    bool operator<(ra_iterator const & rhs) const
    {
        return this->_pos < rhs._pos;
    }

    bool operator>(ra_iterator const & rhs) const
    {
        return this->_pos > rhs._pos;
    }

    bool operator<=(ra_iterator const & rhs) const
    {
        return this->_pos <= rhs._pos;
    }

    bool operator>=(ra_iterator const & rhs) const
    {
        return this->_pos >= rhs._pos;
    }

    //! pre-increment, return updated iterator
    ra_iterator& operator++()
    {
        ++_pos;
        return (*this);
    }

    //! post-increment, return old iterator
    ra_iterator operator++(int)
    {
        ra_iterator cpy{*this};
        ++_pos;
        if (_pos >= _host.size())
            _pos = std::numeric_limits<size_type>::max();
        return cpy;
    }

    //! pre-decrement, return updated iterator
    ra_iterator& operator--()
    {
        --_pos;
        return *this;
    }

    //! post-decrement, return the old iterator
    ra_iterator operator--(int)
    {
        ra_iterator cpy{*this};
        --_pos;
        return cpy;
    }

    //! forward iterator and assignment
    ra_iterator& operator+=(size_type skip)
    {
        _pos += skip;
        if (_pos >= _host.size())
            _pos = std::numeric_limits<size_type>::max();
        return *this;
    }

    //! forward iterator    ra_iterator operator+(size_type) const; //optional
    ra_iterator operator+(size_type skip) const
    {
        ra_iterator cpy{*this};
        cpy._pos = cpy._pos + skip;
        if (cpy._pos >= _host.size())
            cpy._pos = std::numeric_limits<size_type>::max();
        return cpy;
    }

    friend ra_iterator operator+(size_type skip , const ra_iterator& _it)
    {
        ra_iterator cpy{_it};
        cpy._pos = cpy._pos + skip;
        if (cpy._pos >= _it._host.size())
            cpy._pos = std::numeric_limits<size_type>::max();
        return cpy;
    }

    ra_iterator& operator-=(size_type skip)
    {
        _pos -= skip;
        return *this;
    }

    ra_iterator operator-(size_type skip) const
    {
        ra_iterator cpy{*this};
        cpy._pos -= skip;
        return cpy;
    }

    difference_type operator-(ra_iterator lhs) const
    {
        return _pos - lhs._pos;
    }

    pointer operator->() const
    {
        return &(*this);
    }

    reference operator[](size_type const n) const
    {
        return _host[_pos + n];
    }

};

} // namespace seqan3::detail
