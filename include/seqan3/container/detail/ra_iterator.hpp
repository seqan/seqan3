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

#include <range/v3/all.hpp>
#include <seqan3/container/concepts.hpp>

namespace seqan3::detail
{

template <typename container_type, bool is_const=false>
    requires random_access_range_concept<container_type> && sized_range_concept<container_type>
class ra_iterator
{

private:
    // by default iterator is end iterator
    container_type & host{*(container_type*)0};
    using position_type = uint8_t;
    position_type pos{std::numeric_limits<uint8_t>::max()};
public:
    using difference_type = int8_t;
    using value_type = typename container_type::value_type;
    using reference = std::conditional_t<is_const, value_type const &, value_type &>;
    using const_reference = typename container_type::const_reference;
    using pointer = std::conditional_t<is_const, value_type const *, value_type *>;
    using category_type = std::random_access_iterator_tag;

    ra_iterator() {}

    ra_iterator(container_type & host, bool const _at_end = false) : host{host}, pos{0}
    {
        if (_at_end)
            pos = std::numeric_limits<difference_type>::max();
    }

    // copy constructor
    ra_iterator(ra_iterator const & _it) : host{_it.host}, pos{_it.pos} {}

    // copy via assignment
    ra_iterator& operator=(ra_iterator const & rhs)
    {
        pos = rhs.pos;
        return *this;
    }

    ~ra_iterator() = default;

    // return element currently pointed at
    reference operator*() const
    {
        return host[pos];
    }

    // two iterators are equal if their absolute positions are the same
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

    // pre-increment, return updated iterator
    ra_iterator& operator++()
    {
        ++pos;
        return (*this);
    }

    // post-increment, return old iterator
    ra_iterator operator++(int)
    {
        ra_iterator cpy{*this};
        ++pos;
        if (pos >= ranges::size(host))
            pos = std::numeric_limits<difference_type>::max();
        return cpy;
    }

    // pre-decrement, return updated iterator
    ra_iterator& operator--()
    {
        --pos;
        return *this;
    }

    // post-decrement, return the old iterator
    ra_iterator operator--(int)
    {
        ra_iterator cpy{*this};
        --pos;
        return cpy;
    }

    // forward iterator and assignment
    ra_iterator& operator+=(difference_type skip)
    {
        pos += skip;
        if (pos >= ranges::size(host))
            pos = std::numeric_limits<difference_type>::max();
        return *this;
    }

    // forward iterator
    ra_iterator operator+(difference_type skip) const
    {
        ra_iterator cpy{*this};
        cpy.pos = cpy.pos + skip;
        if (cpy.pos >= ranges::size(host))
            cpy.pos = std::numeric_limits<difference_type>::max();
        return cpy;
    }

    friend ra_iterator operator+(difference_type skip , const ra_iterator& _it)
    {
        ra_iterator cpy{_it};
        cpy.pos = cpy.pos + skip;
        if (cpy.pos >= ranges::size(_it.host))
            cpy.pos = std::numeric_limits<difference_type>::max();
        return cpy;
    }

    ra_iterator& operator-=(difference_type skip)
    {
        pos -= skip;
        return *this;
    }

    ra_iterator operator-(difference_type skip) const
    {
        ra_iterator cpy{*this};
        cpy.pos -= skip;
        return cpy;
    }

    difference_type operator-(ra_iterator lhs) const
    {
        return (difference_type)pos - (difference_type)lhs.pos;
    }

    pointer operator->() const
    {
        return &(*this);
    }

    reference operator[](value_type const n) const
    {
        return host[pos + n];
    }

};

} // namespace seqan3::detail

static_assert(seqan3::random_access_iterator_concept<seqan3::detail:ra_iterator>);
