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

/*!\file container/detail/random_access_iterator.hpp
 * \brief Random access iterator for const and non const containers.
 * \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 * \ingroup iterator
 */

#pragma once

#include <iomanip>
#include <limits>
#include <type_traits>

#include <range/v3/all.hpp>

#include <seqan3/core/concept/iterator.hpp>
#include <seqan3/range/concept.hpp>

namespace seqan3::detail
{

// todo: replace is_const flag by container value function
template <typename container_type, bool is_const=false>
    requires random_access_range_concept<container_type> && sized_range_concept<container_type>
class random_access_iterator
{

private:

    // TODO: query container's diff type with to be provided metafunctions
    using difference_type_container = int8_t;
    using size_type_container = uint8_t;

    typename std::add_pointer<container_type>::type host{nullptr};

    using position_type = size_type_container;
    position_type pos{static_cast<position_type>(0)};

public:
    //! value type for distances between iterators
    using difference_type = difference_type_container;
    //! value type obtained by dereferencing the iterator
    using value_type = typename container_type::value_type;

    using reference = std::conditional_t<is_const, value_type const &, value_type &>;
    using const_reference = typename container_type::const_reference;
    using pointer = std::conditional_t<is_const, value_type const *, value_type *>;
    //using iterator_category = std::random_access_iterator_tag;

    constexpr random_access_iterator() = default;

    //! construct by host, default position pointer with 0
    constexpr random_access_iterator(container_type & host) : host{&host} {}

    //! construct by host and explicit position
    constexpr random_access_iterator(container_type & host, position_type pos) : host{&host}, pos{pos} {}

    //! copy constructor
    constexpr random_access_iterator(random_access_iterator const &) = default;

    //! copy construction via assignment
    constexpr random_access_iterator & operator=(random_access_iterator const &) = default;

    //! move constructor
    constexpr random_access_iterator (random_access_iterator &&) = default;

    //! move assignment
    constexpr random_access_iterator & operator=(random_access_iterator &&) = default;

    //! use default deconstructor
    ~random_access_iterator() = default;

    //! dereference operator returns element currently pointed at
    value_type operator*() const
    {
        return (*host)[pos];
    }

    //! two iterators are equal if their absolute positions are the same
    bool operator==(random_access_iterator const & rhs) const
    {
        return this->pos == rhs.pos;
    }

    //! iterator inequality comparison refers to their positions
    bool operator!=(random_access_iterator const & rhs) const
    {
        return this->pos != rhs.pos;
    }

    //! iterator comparison refers to their positions
    bool operator<(random_access_iterator const & rhs) const
    {
        return this->pos < rhs.pos;
    }

    //! iterator comparison refers to their positions
    bool operator>(random_access_iterator const & rhs) const
    {
        return this->pos > rhs.pos;
    }

    //! iterator comparison refers to their positions
    bool operator<=(random_access_iterator const & rhs) const
    {
        return this->pos <= rhs.pos;
    }

    //! iterator comparison refers to their positions
    bool operator>=(random_access_iterator const & rhs) const
    {
        return this->pos >= rhs.pos;
    }

    //! pre-increment, return updated iterator
    random_access_iterator & operator++()
    {
        ++pos;
        return (*this);
    }

    //! post-increment, return previous iterator state
    random_access_iterator operator++(int)
    {
        random_access_iterator cpy{*this};
        ++pos;
        return cpy;
    }

    //! pre-decrement, return updated iterator
    random_access_iterator& operator--()
    {
        --pos;
        return *this;
    }

    //! post-decrement, return previous iterator state
    random_access_iterator operator--(int)
    {
        random_access_iterator cpy{*this};
        --pos;
        return cpy;
    }

    //! forward this iterator
    random_access_iterator& operator+=(difference_type skip)
    {
        pos += skip;
        return *this;
    }

    //! forward copy of this iterator
    random_access_iterator operator+(difference_type skip) const
    {
        return random_access_iterator{*this->host, static_cast<position_type>(pos + skip)};
    }

    //! non-member operator+ delegates to non-friend operator+
    friend random_access_iterator operator+(difference_type skip , const random_access_iterator& it)
    {
        return it + skip;
    }

    //! decrement iterator by skip
    random_access_iterator& operator-=(difference_type skip)
    {
        pos -= skip;
        return *this;
    }

    //! return decremented copy of this iterator
    random_access_iterator operator-(difference_type skip) const
    {
        return random_access_iterator{*this->host, static_cast<position_type>(pos - skip)};
    }

    //! non-member operator- delegates to non-friend operator-
    friend random_access_iterator operator-(difference_type skip , const random_access_iterator& it)
    {
        return it - skip;
    }

    //! return offset between this and remote iterator's position
    difference_type operator-(random_access_iterator lhs) const
    {
        return static_cast<difference_type>(pos) - static_cast<difference_type>(lhs.pos);
    }

    //! return pointer to this iterator
    pointer operator->() const
    {
        return &this->host[pos]; //&(*this);
    }

    //! return underlying container value currently pointed at
    reference operator[](position_type const n) const
    {
        return (*host)[pos + n];
    }
};

} // namespace seqan3::detail

//static_assert(seqan3::random_access_iterator_concept<seqan3::detail:random_access_iterator>);
//static_assert(seqan3::readable_concept<seqan3::detail:random_access_iterator<std::vector<int>>>);
