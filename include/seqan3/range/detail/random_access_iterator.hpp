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
 */

#pragma once

#include <cassert>
#include <type_traits>

#include <range/v3/utility/iterator_traits.hpp>
#include <range/v3/range_traits.hpp>

#include <seqan3/std/concept/range.hpp>
#include <seqan3/std/concept/iterator.hpp>

namespace seqan3::detail
{

/*!\brief A CRTP base template for creating random access iterators.
 * \tparam range_type The data structure on which the iterator operates, e.g. `std::vector<int>`.
 * \tparam derived_t_template The derived type template, see below.
 * \implements seqan3::random_access_iterator_concept
 * \ingroup range
 *
 * The iterator makes certain assumptions on the `range_type`, but does not formally require
 * it to satisfy the seqan3::random_access_range_concept, because the iterator itself may be
 * a requirement for this.
 *
 * Actual functionality of this iterator is realised via the host range's [] operator and member type definitions, i.e.
 * you need to implement these in the range before you can make use of this iterator.
 *
 * Since the CRTP parameter is in fact a template template, CRTP instantiation looks a little different, e.g.:
 * ```cpp
 * template <typename range_type>
 * class random_access_iterator : public random_access_iterator_base<range_type, random_access_iterator>
 * {
 * //...
 * };
 * ```
 * and not:
 * ```cpp
 * template <typename range_type>
 * class random_access_iterator : public random_access_iterator_base<range_type, random_access_iterator<range_type>>
 * {
 * //...
 * };
 * ```
 */
template <typename range_type, template <typename...> typename derived_t_template>
class random_access_iterator_base
{
protected:
    //!\brief Iterator stores pointer to underlying container structure.
    typename std::add_pointer_t<range_type> host{nullptr};
    //!\brief Use container's size_type as a position.
    using position_type =  ranges::v3::size_type_t<range_type>;
    //!\brief Store position index for container.
    position_type pos{static_cast<position_type>(0)};

    //!\brief This friend declaration is required to allow non-const to const-construction.
    template <typename range_type2, template <typename...> typename derived_t_template2>
    //!\cond
        requires std::is_const_v<range_type> && !std::is_const_v<range_type2> &&
                 std::is_same_v<std::remove_const_t<range_type>, range_type2> &&
                 std::is_same_v<derived_t_template2, derived_t_template>
    //!\endcond
    friend class random_access_iterator_base;

    //!\brief Because this is CRTP, we know the full derived type:
    using derived_t = derived_t_template<range_type>;

public:
    //!\brief Type for distances between iterators.
    using difference_type = typename range_type::difference_type; // TODO should be range_ but is broken in ranges
    //!\brief Value type of container elements.
    using value_type = typename range_type::value_type;
    //!\brief Use reference type defined by container.
    using reference = std::conditional_t<std::is_const_v<range_type>,
                                         typename range_type::const_reference,
                                         typename range_type::reference>;
    //!\brief Use const reference type provided by container.
    using const_reference = typename range_type::const_reference; //TODO: there is no metafunction for this, yet :o
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
    explicit constexpr random_access_iterator_base(range_type & host) noexcept : host{&host} {}
    //!\brief Construct by host and explicit position.
    constexpr random_access_iterator_base(range_type & host, position_type const pos) noexcept :
        host{&host}, pos{pos}
    {}

    //!\brief Constructor for const version from non-const version.
    template <typename range_type2>
    //!\cond
        requires std::is_const_v<range_type> && !std::is_const_v<range_type2> &&
                 std::is_same_v<std::remove_const_t<range_type>, range_type2>
    //!\endcond
    constexpr random_access_iterator_base(random_access_iterator_base<range_type2, derived_t_template> const & rhs) noexcept :
        host{rhs.host}, pos{rhs.pos}
    {}
    //!\}

    /*!\name Comparison operators
     * \brief Compare iterators by position. seqan3::detail::random_access_iterator_base operators are used unless
     * specialised in derived type.
     * \{
     */
    template <typename range_type2>
    //!\cond
        requires std::is_same_v<std::remove_const_t<range_type>, std::remove_const_t<range_type2>>
    //!\endcond
    constexpr bool operator==(random_access_iterator_base<range_type2, derived_t_template> const & rhs) const noexcept
    {
        return pos == rhs.pos;
    }

    template <typename range_type2>
    //!\cond
        requires std::is_same_v<std::remove_const_t<range_type>, std::remove_const_t<range_type2>>
    //!\endcond
    constexpr bool operator!=(random_access_iterator_base<range_type2, derived_t_template> const & rhs) const noexcept
    {
        return pos != rhs.pos;
    }

    template <typename range_type2>
    //!\cond
        requires std::is_same_v<std::remove_const_t<range_type>, std::remove_const_t<range_type2>>
    //!\endcond
    constexpr bool operator<(random_access_iterator_base<range_type2, derived_t_template> const & rhs) const noexcept
    {
        return static_cast<bool>(pos < rhs.pos);
    }

    template <typename range_type2>
    //!\cond
        requires std::is_same_v<std::remove_const_t<range_type>, std::remove_const_t<range_type2>>
    //!\endcond
    constexpr bool operator>(random_access_iterator_base<range_type2, derived_t_template> const & rhs) const noexcept
    {
        return pos > rhs.pos;
    }

    template <typename range_type2>
    //!\cond
        requires std::is_same_v<std::remove_const_t<range_type>, std::remove_const_t<range_type2>>
    //!\endcond
    constexpr bool operator<=(random_access_iterator_base<range_type2, derived_t_template> const & rhs) const noexcept
    {
        return pos <= rhs.pos;
    }

    template <typename range_type2>
    //!\cond
        requires std::is_same_v<std::remove_const_t<range_type>, std::remove_const_t<range_type2>>
    //!\endcond
    constexpr bool operator>=(random_access_iterator_base<range_type2, derived_t_template> const & rhs) const noexcept
    {
        return pos >= rhs.pos;
    }
    //!\}

    /*!\name Arithmetic operators
     * \brief seqan3::detail::random_access_iterator_base operators are used unless specialised in derived type.
     * \{
    */
    //!\brief Pre-increment, return updated iterator.
    constexpr derived_t & operator++() noexcept
    {
        ++pos;
        return *this_derived();
    }

    //!\brief Post-increment, return previous iterator state.
    constexpr derived_t operator++(int) noexcept
    {
        derived_t cpy{*this_derived()};
        ++pos;
        return cpy;
    }

    //!\brief Pre-decrement, return updated iterator.
    constexpr derived_t & operator--() noexcept
    {
        --pos;
        return *this_derived();
    }

    //!\brief Post-decrement, return previous iterator state.
    constexpr derived_t operator--(int) noexcept
    {
        derived_t cpy{*this_derived()};
        --pos;
        return cpy;
    }

    //!\brief Forward this iterator.
    constexpr derived_t & operator+=(difference_type const skip) noexcept
    {
        pos += skip;
        return *this_derived();
    }

    //!\brief Forward copy of this iterator.
    constexpr derived_t operator+(difference_type const skip) const noexcept
    {
        derived_t cpy{*this_derived()};
        return cpy += skip;
    }

    //!\brief Non-member operator+ delegates to non-friend operator+.
    constexpr friend derived_t operator+(difference_type const skip, derived_t const & it) noexcept
    {
        return it + skip;
    }

    //!\brief Decrement iterator by skip.
    constexpr derived_t & operator-=(difference_type const skip) noexcept
    {
        pos -= skip;
        return *this_derived();
    }

    //!\brief Return decremented copy of this iterator.
    constexpr derived_t operator-(difference_type const skip) const noexcept
    {
        derived_t cpy{*this_derived()};
        return cpy -= skip;
    }

    //!\brief Non-member operator- delegates to non-friend operator-.
    constexpr friend derived_t operator-(difference_type const skip, derived_t const & it) noexcept
    {
        return it - skip;
    }

    //!\brief Return offset between this and remote iterator's position.
    constexpr difference_type operator-(derived_t const lhs) const noexcept
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

private:

    //!\brief Cast this to derived type.
    derived_t* this_derived()
    {
        return static_cast<derived_t*>(this);
    }

    //!\copydoc this_derived
    derived_t const * this_derived() const
    {
        return static_cast<derived_t const *>(this);
    }
};

/*!\brief A generic random access iterator that delegates most operations to the range.
 * \tparam range_type The data structure on which the iterator operates, e.g. `std::vector<int>`.
 * \ingroup range
 *
 * The iterator makes certain assumptions on the `range_type`, but does not formally require
 * it to satisfy the seqan3::random_access_range_concept, because the iterator itself may be
 * a requirement for this.
 */
template <typename range_type>
class random_access_iterator :
    public random_access_iterator_base<range_type, random_access_iterator>
{
private:
    //!\brief Shortcut for the base class.
    using base = random_access_iterator_base<range_type, random_access_iterator>;
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
};

} // namespace seqan3::detail
