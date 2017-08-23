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

/*!\brief A CRTP base template for creating random access iterators.
 * \tparam container_type The data structure on which the iterator operates, e.g. `std::vector<int>`.
 * \tparam derived_type The derived type.
 * \implements seqan3::random_access_iterator_concept
 * \ingroup range
 *
 * The iterator makes certain assumptions on the `container_type`, but does not formally require
 * it to satisfy the seqan3::random_access_range_concept, because the iterator itself may be
 * a requirement for this.
 */
template <typename container_type, typename derived_type>
class random_access_iterator_base
{
protected:
    //!\privatesection
    //!\brief Iterator stores pointer to underlying container structure.
    typename std::add_pointer_t<container_type> host{nullptr};
    //!\brief Use container's size_type as a position.
    using position_type =  ranges::v3::size_type_t<container_type>;
    //!\brief Store position index for container.
    position_type pos{static_cast<position_type>(0)};

    //!\brief This friend declaration is required to allow non-const to const-construction.
    template <typename t, typename derived2_type>
        requires !std::is_const_v<t> && std::is_const_v<container_type>
    friend class random_access_iterator_base;

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
    //!\brief Default constructor.
    constexpr random_access_iterator_base() = default;
    //!\brief Copy constructor.
    constexpr random_access_iterator_base(random_access_iterator_base const &) = default;
    //!\brief Copy construction via assignment.
    constexpr random_access_iterator_base & operator=(random_access_iterator_base const &) = default;
    //!\brief Move constructor.
    constexpr random_access_iterator_base (random_access_iterator_base &&) = default;
    //!\brief Move assignment.
    constexpr random_access_iterator_base & operator=(random_access_iterator_base &&) = default;
    //!\brief Use default deconstructor.
    ~random_access_iterator_base() = default;

    //!\brief Construct by host, default position pointer with 0.
    explicit constexpr random_access_iterator_base(container_type & host) noexcept : host{&host} {}
    //!\brief Construct by host and explicit position.
    constexpr random_access_iterator_base(container_type & host, position_type const pos) noexcept :
        host{&host}, pos{pos}
    {}

// TODO make this inheritable
//     //!\brief Constructor for const version from non-const version.
//     template <typename t>
//     //!\cond
//         requires std::is_same_v<t, typename iterator_on_const_container<derived_type>::type> &&
//                  !std::is_same_v<std::remove_const_t<container_type>, container_type>
//     //!\endcond
//     constexpr random_access_iterator_base(t const & rhs) noexcept :
//         host{rhs.host}, pos{rhs.pos}
//     {}
    //!\}

    /*!\name Comparison operators
     * \brief seqan3::detail::random_access_iterator_base operators are used unless specialised in derived type.
     * \{
     */
    constexpr bool operator==(derived_type const & rhs) const noexcept
    {
        return pos == rhs.pos;
    }

    constexpr bool operator!=(derived_type const & rhs) const noexcept
    {
        return pos != rhs.pos;
    }

    constexpr bool operator<(derived_type const & rhs) const noexcept
    {
        return static_cast<bool>(pos < rhs.pos);
    }

    constexpr bool operator>(derived_type const & rhs) const noexcept
    {
        return pos > rhs.pos;
    }

    constexpr bool operator<=(derived_type const & rhs) const noexcept
    {
        return pos <= rhs.pos;
    }

    constexpr bool operator>=(derived_type const & rhs) const noexcept
    {
        return pos >= rhs.pos;
    }
    //!\}

    /*!\name Arithmetic operators
     * \brief seqan3::detail::random_access_iterator_base operators are used unless specialised in derived type.
     * \{
    */
    //!\brief Pre-increment, return updated iterator.
    constexpr derived_type & operator++() noexcept
    {
        ++pos;
        return *(static_cast<derived_type*>(this));
    }

    //!\brief Post-increment, return previous iterator state.
    constexpr derived_type operator++(int) noexcept
    {
        derived_type cpy{*(static_cast<derived_type*>(this))};
        ++pos;
        return cpy;
    }

    //!\brief Pre-decrement, return updated iterator.
    constexpr derived_type & operator--() noexcept
    {
        --pos;
        return *(static_cast<derived_type*>(this));
    }

    //!\brief Post-decrement, return previous iterator state.
    constexpr derived_type operator--(int) noexcept
    {
        derived_type cpy{*(static_cast<derived_type*>(this))};
        --pos;
        return cpy;
    }

    //!\brief Forward this iterator.
    constexpr derived_type & operator+=(difference_type const skip) noexcept
    {
        pos += skip;
        return *(static_cast<derived_type*>(this));
    }

    //!\brief Forward copy of this iterator.
    constexpr derived_type operator+(difference_type const skip) const noexcept
    {
        derived_type cpy{*(static_cast<derived_type const *>(this))};
        return cpy += skip;
    }

    //!\brief Non-member operator+ delegates to non-friend operator+.
    constexpr friend derived_type operator+(difference_type const skip , derived_type const & it) noexcept
    {
        return it + skip;
    }

    //!\brief Decrement iterator by skip.
    constexpr derived_type & operator-=(difference_type const skip) noexcept
    {
        pos -= skip;
        return *(static_cast<derived_type*>(this));
    }

    //!\brief Return decremented copy of this iterator.
    constexpr derived_type operator-(difference_type const skip) const noexcept
    {
        derived_type cpy{*(static_cast<derived_type const *>(this))};
        return cpy -= skip;
    }

    //!\brief Non-member operator- delegates to non-friend operator-.
    constexpr friend derived_type operator-(difference_type const skip, derived_type const & it) noexcept
    {
        return it - skip;
    }

    //!\brief Return offset between this and remote iterator's position.
    constexpr difference_type operator-(derived_type const lhs) const noexcept
    {
        return static_cast<difference_type>(pos - lhs.pos);
    }
    //!\}

    /*!\name Reference/Dereference operators
     * \brief seqan3::detail::random_access_iterator_base operators are used unless specialised in derived type.
     * \{
    */
    //!\brief Dereference operator returns element currently pointed at.
    constexpr reference operator*() const noexcept(noexcept((*host)[pos]))
    {
        return (*host)[pos];
    }

    //!\brief Return pointer to this iterator.
    constexpr pointer operator->() const noexcept(noexcept((&host)[pos]))
    {
        return &host[pos];
    }

    //!\brief Return underlying container value currently pointed at.
    constexpr reference operator[](position_type const n) const noexcept(noexcept((*host)[pos+n]))
    {
        return (*host)[pos + n];
    }
    //!\}
};

/*!\brief A generic random access iterator that delegates most operations to the range.
 * \tparam container_type The data structure on which the iterator operates, e.g. `std::vector<int>`.
 * \ingroup range
 *
 * The iterator makes certain assumptions on the `container_type`, but does not formally require
 * it to satisfy the seqan3::random_access_range_concept, because the iterator itself may be
 * a requirement for this.
 */
template <typename container_type>
class random_access_iterator :
    public random_access_iterator_base<container_type, random_access_iterator<container_type>>
{
private:
    //!\brief Shortcut for the base class.
    using base = random_access_iterator_base<container_type, random_access_iterator<container_type>>;
    //!\brief Import from base class.
    using typename base::position_type;

public:
    /*!\name Member types
     * \brief Make the parent's member types visible.
     * \{
     */
    using typename base::difference_type;
    using typename base::value_type;
    using typename base::reference;
    using typename base::const_reference;
    using typename base::pointer;
    using typename base::iterator_category;
    //!\}

    //!\brief Import the parent's constructors.
    using base::base;
    //!\brief Befriend the parent.
    friend base;

    //!\brief This friend declaration is required to allow non-const to const-construction.
    template <typename t>
        requires std::is_const_v<t> && !std::is_const_v<container_type>
    friend class random_access_iterator;

    //!\brief Constructor for const version from non-const version.
    template <typename t>
    //!\cond
        requires !std::is_const_v<t> && std::is_const_v<container_type>
    //!\endcond
    constexpr random_access_iterator(random_access_iterator<t> const & rhs) noexcept :
        random_access_iterator{*rhs.host, rhs.pos}
    {}
};

} // namespace seqan3::detail

static_assert(seqan3::random_access_iterator_concept<seqan3::detail::random_access_iterator<std::vector<int>>>);
