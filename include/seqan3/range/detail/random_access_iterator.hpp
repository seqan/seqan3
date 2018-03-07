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

/*!\file
 * \brief Provides the seqan3::detail::random_access_iterator class.
 * \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 * \ingroup range
 */

#pragma once

#include <cassert>
#include <type_traits>

#include <range/v3/utility/iterator_traits.hpp>

#include <seqan3/core/concept/iterator.hpp>
#include <seqan3/range/concept.hpp>

namespace seqan3::detail
{

/*!\brief Implementation of a random access iterator on an input container pointer.
 * \tparam container_type The data structure on which the iterator operates, e.g. `std::vector<int>`.
 *
 * No iterator operation will modify the container. Arithmetic and boolean
 * operations are applied to the iterator positions, not the corresponding values
 * of their containers.
 *
 * The iterator makes certain assumptions on the `container_type`, but does not formally require
 * it to satisfy the seqan3::random_access_range_concept, because the iterator itself is
 * a requirement for this.
 */
template <typename container_type>
class random_access_iterator
{

  private:
    //!\brief Iterator stores pointer to underlying container structure.
    typename std::add_pointer_t<container_type> host{ nullptr };
    //!\brief Use container's size_type as a position.
    using position_type = ranges::v3::size_type_t<container_type>;
    //!\brief Store position index for container.
    position_type pos{ static_cast<position_type>(0) };

    //!\brief This friend declaration is required to allow non-const to const-construction.
    template <typename t>
    requires std::is_same_v<std::remove_const_t<container_type>, t> &&
      std::is_same_v<container_type, std::add_const_t<t>> friend class random_access_iterator;

  public:
    //!\brief Type for distances between iterators.
    using difference_type = ranges::v3::difference_type_t<container_type>;
    //!\brief Value type of container elements.
    using value_type = ranges::v3::value_type_t<container_type>;
    //!\brief Use reference type defined by container.
    using reference = std::conditional_t<std::is_const_v<container_type>,
                                         typename container_type::const_reference,
                                         typename container_type::reference>;
    //!\brief Use const reference type provided by container.
    using const_reference = typename container_type::const_reference;
    //!\brief Pointer type is pointer of container element type.
    using pointer = value_type *;
    //!\brief Tag this class as a random access iterator.
    using iterator_category = std::random_access_iterator_tag;

    /*!\name Constructors/Destructors
     * \{
    */
    // \brief Default constructor.
    constexpr random_access_iterator() = default;

    //!\brief Construct by host, default position pointer with 0.
    explicit constexpr random_access_iterator(container_type & host) noexcept
      : host{ &host }
    {
    }

    //!\brief Construct by host and explicit position.
    constexpr random_access_iterator(container_type & host, position_type const pos) noexcept
      : host{ &host }
      , pos{ pos }
    {
    }

    //!\brief Copy constructor.
    constexpr random_access_iterator(random_access_iterator const &) = default;

    //!\brief Copy construction via assignment.
    constexpr random_access_iterator & operator=(random_access_iterator const &) = default;

    //!\brief Move constructor.
    constexpr random_access_iterator(random_access_iterator &&) = default;

    //!\brief Move assignment.
    constexpr random_access_iterator & operator=(random_access_iterator &&) = default;

    //!\brief Use default deconstructor.
    ~random_access_iterator() = default;

    //!\brief Constructor for const version from non-const version.
    template <typename t>
    //!\cond
    requires std::is_same_v<std::remove_const_t<container_type>, t> &&
      std::is_same_v<container_type, std::add_const_t<t>>
      //!\endcond
      constexpr random_access_iterator(random_access_iterator<t> const & rhs) noexcept
      : host{ rhs.host }
      , pos{ rhs.pos }
    {
    }
    //!\}

    /*!\name Comparison operators
     * \brief Compares only the absolute position of two iterators.
     * \{
     */
    constexpr bool operator==(random_access_iterator const & rhs) const noexcept { return pos == rhs.pos; }

    constexpr bool operator!=(random_access_iterator const & rhs) const noexcept { return pos != rhs.pos; }

    constexpr bool operator<(random_access_iterator const & rhs) const noexcept
    {
        return static_cast<bool>(pos < rhs.pos);
    }

    constexpr bool operator>(random_access_iterator const & rhs) const noexcept { return pos > rhs.pos; }

    constexpr bool operator<=(random_access_iterator const & rhs) const noexcept { return pos <= rhs.pos; }

    constexpr bool operator>=(random_access_iterator const & rhs) const noexcept { return pos >= rhs.pos; }
    //!\}

    /*!\name Arithmetic operators
     * \{
    */
    //!\brief Pre-increment, return updated iterator.
    constexpr random_access_iterator & operator++() noexcept
    {
        ++pos;
        return (*this);
    }

    //!\brief Post-increment, return previous iterator state.
    constexpr random_access_iterator operator++(int)noexcept
    {
        random_access_iterator cpy{ *this };
        ++pos;
        return cpy;
    }

    //!\brief Pre-decrement, return updated iterator.
    constexpr random_access_iterator & operator--() noexcept
    {
        --pos;
        return *this;
    }

    //!\brief Post-decrement, return previous iterator state.
    constexpr random_access_iterator operator--(int)noexcept
    {
        random_access_iterator cpy{ *this };
        --pos;
        return cpy;
    }

    //!\brief Forward this iterator.
    constexpr random_access_iterator & operator+=(difference_type const skip) noexcept
    {
        pos += skip;
        return *this;
    }

    //!\brief Forward copy of this iterator.
    constexpr random_access_iterator operator+(difference_type const skip) const noexcept
    {
        return random_access_iterator{ *host, static_cast<position_type>(pos + skip) };
    }

    //!\brief Non-member operator+ delegates to non-friend operator+.
    friend random_access_iterator operator+(difference_type const skip, random_access_iterator const & it) noexcept
    {
        return it + skip;
    }

    //!\brief Decrement iterator by skip.
    constexpr random_access_iterator & operator-=(difference_type const skip) noexcept
    {
        pos -= skip;
        return *this;
    }

    //!\brief Return decremented copy of this iterator.
    constexpr random_access_iterator operator-(difference_type const skip) const noexcept
    {
        return random_access_iterator{ *host, static_cast<position_type>(pos - skip) };
    }

    //!\brief Non-member operator- delegates to non-friend operator-.
    constexpr friend random_access_iterator operator-(difference_type const skip,
                                                      random_access_iterator const & it) noexcept
    {
        return it - skip;
    }

    //!\brief Return offset between this and remote iterator's position.
    constexpr difference_type operator-(random_access_iterator const lhs) const noexcept
    {
        return static_cast<difference_type>(pos - lhs.pos);
    }
    //!\}

    /*!\name Reference/Dereference operators
     * \{
    */
    //!\brief Dereference operator returns element currently pointed at.
    constexpr reference operator*() noexcept(noexcept((*host)[pos])) { return (*host)[pos]; }

    //!\brief Return pointer to this iterator.
    constexpr pointer operator->() const noexcept(noexcept((&host)[pos])) { return &host[pos]; }

    //!\brief Return underlying container value currently pointed at.
    constexpr reference operator[](position_type const n) const noexcept(noexcept((*host)[pos + n]))
    {
        return (*host)[pos + n];
    }
    //!\}
};

} // namespace seqan3::detail

static_assert(static_cast<bool>(ranges::concepts::models<ranges::concepts::RandomAccessIterator,
                                                         seqan3::detail::random_access_iterator<std::vector<int>>>()));
