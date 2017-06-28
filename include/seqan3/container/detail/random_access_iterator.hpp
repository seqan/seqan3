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
#include <cassert>

#include <range/v3/all.hpp>

#include <seqan3/core/concept/iterator.hpp>
#include <seqan3/range/concept.hpp>

namespace seqan3::detail
{

template <typename container_type>
    requires random_access_range_concept<container_type> && sized_range_concept<container_type>
class random_access_iterator
{

private:
    //! \brief Iterator stores pointer to underlying container structure
    typename std::add_pointer_t<container_type> host{nullptr};
    using position_type =  ranges::v3::size_type_t<container_type>;
    //! \brief Store position index for container
    position_type pos{static_cast<position_type>(0)};

public:
    //! \brief Type for distances between iterators
    using difference_type = ranges::v3::difference_type_t<container_type>;
    //! \brief Value type of container elements
    using value_type = ranges::v3::value_type_t<container_type>;
    //! \brief Reference types defined by container
    using reference = typename container_type::reference;
    using const_reference = typename container_type::const_reference;
    //! \brief Pointer type is pointer of container element type
    using pointer = value_type *;
    using iterator_category = std::random_access_iterator_tag;

    // \brief Default constructor
    constexpr random_access_iterator() = default;

    //! \brief Construct by host, default position pointer with 0
    constexpr random_access_iterator(container_type & host) : host{&host} {}

    //! \brief Construct by host and explicit position
    constexpr random_access_iterator(container_type & host, position_type pos) : host{&host}, pos{pos} {}

    //! \brief Copy constructor
    constexpr random_access_iterator(random_access_iterator const &) = default;

    //! \brief Copy construction via assignment
    constexpr random_access_iterator & operator=(random_access_iterator const &) = default;

    //! \brief Move constructor
    constexpr random_access_iterator (random_access_iterator &&) = default;

    //! \brief Move assignment
    constexpr random_access_iterator & operator=(random_access_iterator &&) = default;

    //! \brief Use default deconstructor
    ~random_access_iterator() = default;

    //! \brief Dereference operator returns element currently pointed at
    reference operator*()
    {
        return (reference)(*host)[pos];
    }

    //! \brief Two iterators are equal if their absolute positions are the same
    bool operator==(random_access_iterator const & rhs) const
    {
        return this->pos == rhs.pos;
    }

    //! \brief Iterator inequality comparison refers to their positions
    bool operator!=(random_access_iterator const & rhs) const
    {
        return this->pos != rhs.pos;
    }

    //! \brief Iterator comparison refers to their positions
    bool operator<(random_access_iterator const & rhs) const
    {
        return static_cast<bool>(this->pos < rhs.pos);
    }

    //! \brief Iterator comparison refers to their positions
    bool operator>(random_access_iterator const & rhs) const
    {
        return this->pos > rhs.pos;
    }

    //! \brief Iterator comparison refers to their positions
    bool operator<=(random_access_iterator const & rhs) const
    {
        return this->pos <= rhs.pos;
    }

    //! \brief Iterator comparison refers to their positions
    bool operator>=(random_access_iterator const & rhs) const
    {
        return this->pos >= rhs.pos;
    }

    //! \brief Pre-increment, return updated iterator
    random_access_iterator & operator++()
    {
        ++pos;
        return (*this);
    }

    //! \brief Post-increment, return previous iterator state
    random_access_iterator operator++(int)
    {
        random_access_iterator cpy{*this};
        ++pos;
        return cpy;
    }

    //! \brief Pre-decrement, return updated iterator
    random_access_iterator& operator--()
    {
        --pos;
        return *this;
    }

    //! \brief Post-decrement, return previous iterator state
    random_access_iterator operator--(int)
    {
        random_access_iterator cpy{*this};
        --pos;
        return cpy;
    }

    //! \brief Forward this iterator
    random_access_iterator& operator+=(difference_type skip)
    {
        pos += skip;
        return *this;
    }

    //! \brief Forward copy of this iterator
    random_access_iterator operator+(difference_type skip) const
    {
        return random_access_iterator{*this->host, static_cast<position_type>(pos + skip)};
    }

    //! \brief Non-member operator+ delegates to non-friend operator+
    friend random_access_iterator operator+(difference_type skip , const random_access_iterator& it)
    {
        return it + skip;
    }

    //! \brief Decrement iterator by skip
    random_access_iterator& operator-=(difference_type skip)
    {
        pos -= skip;
        return *this;
    }

    //! \brief Return decremented copy of this iterator
    random_access_iterator operator-(difference_type skip) const
    {
        return random_access_iterator{*this->host, static_cast<position_type>(this->pos - skip)};
    }

    //! \brief Non-member operator- delegates to non-friend operator-
    friend random_access_iterator operator-(difference_type skip , const random_access_iterator& it)
    {
        return it - skip;
    }

    //! \brief Return offset between this and remote iterator's position
    difference_type operator-(random_access_iterator lhs) const
    {
        return static_cast<difference_type>(this->pos - lhs.pos);
    }

    //! \brief Return pointer to this iterator
    pointer operator->() const
    {
        return &this->host[pos]; //&(*this);
    }

    //! \brief Return underlying container value currently pointed at
    reference operator[](position_type const n) const
    {
        return (reference)(*host)[pos + n];
    }
};

} // namespace seqan3::detail

static_assert(static_cast<bool>(ranges::concepts::models<ranges::concepts::RandomAccessIterator, seqan3::detail::random_access_iterator<std::vector<int>>>()));
