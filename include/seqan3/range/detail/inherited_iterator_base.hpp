// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the seqan3::detail::inherited_iterator_base template.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <cassert>
#include <type_traits>

#include <range/v3/utility/iterator_traits.hpp>
#include <range/v3/range_traits.hpp>

#include <seqan3/std/ranges>
#include <seqan3/std/iterator>

namespace seqan3::detail
{

/*!\brief A CRTP base template for creating iterators that inherit from other iterators.
 * \tparam derived_t The CRTP specialisation.
 * \tparam base_t    The type to inherit from; must satisfy std::Iterator.
 * \implements std::Iterator
 * \ingroup range
 *
 * \details
 *
 * This template enables you to inherit from another iterator and just overload those functions
 * that you wish to change.
 *
 * ### Example
 *
 * You may want to have an iterator that skips those elements in an int-vector that are not divisable by two:
 *
 * \snippet test/unit/range/detail/inherited_iterator_base_test.cpp inherited_iterator_base desired
 *
 * You could define it like this:
 *
 * \snippet test/unit/range/detail/inherited_iterator_base_test.cpp inherited_iterator_base def
 */
template <typename derived_t, std::Iterator base_t>
class inherited_iterator_base : public base_t
{
public:
    /*!\name Associated types
     * \brief All are derived from the base_t.
     * \{
     */
    using difference_type       = typename std::iterator_traits<base_t>::difference_type;
    using value_type            = typename std::iterator_traits<base_t>::value_type;
    using reference             = typename std::iterator_traits<base_t>::reference;
    using pointer               = typename std::iterator_traits<base_t>::pointer;
    using iterator_category     = typename std::iterator_traits<base_t>::iterator_category;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    inherited_iterator_base() = default;
    constexpr inherited_iterator_base(inherited_iterator_base const & rhs) = default;
    constexpr inherited_iterator_base(inherited_iterator_base && rhs) = default;
    constexpr inherited_iterator_base & operator=(inherited_iterator_base const & rhs) = default;
    constexpr inherited_iterator_base & operator=(inherited_iterator_base && rhs) = default;
    ~inherited_iterator_base() = default;

    inherited_iterator_base(base_t it) : base_t{it} {}
    //!\}

    /*!\name Comparison operators
     * \brief Unless specialised in derived_type, all operators perform base_t's operator and cast to derived_t.
     * \{
     */
    constexpr bool operator==(derived_t const & rhs) const noexcept(noexcept(base_t{} == base_t{}))
    //!\cond
        requires std::EqualityComparable<base_t>
    //!\endcond
    {
        return *this_to_base() == static_cast<base_t>(rhs);
    }

    constexpr bool operator!=(derived_t const & rhs) const noexcept(noexcept(base_t{} == base_t{}))
    //!\cond
        requires std::EqualityComparable<base_t>
    //!\endcond
    {
        return !(*this == rhs);
    }

    constexpr bool operator<(derived_t const & rhs) const noexcept(noexcept(base_t{} < base_t{}))
    //!\cond
        requires std::StrictTotallyOrdered<base_t>
    //!\endcond
    {
        return *this_to_base() < static_cast<base_t>(rhs);
    }

    constexpr bool operator>(derived_t const & rhs) const noexcept(noexcept(base_t{} > base_t{}))
    //!\cond
        requires std::StrictTotallyOrdered<base_t>
    //!\endcond
    {
        return *this_to_base() > static_cast<base_t>(rhs);
    }

    constexpr bool operator<=(derived_t const & rhs) const noexcept(noexcept(base_t{} < base_t{}))
    //!\cond
        requires std::StrictTotallyOrdered<base_t>
    //!\endcond
    {
        return !(*this > rhs);
    }

    constexpr bool operator>=(derived_t const & rhs) const noexcept(noexcept(base_t{} < base_t{}))
    //!\cond
        requires std::StrictTotallyOrdered<base_t>
    //!\endcond
    {
        return !(*this < rhs);
    }
    //!\}

    /*!\name Arithmetic operators
     * \brief Unless specialised in derived_type, all operators perform base_t's operator and cast to derived_t.
     * \{
    */
    //!\brief Pre-increment, return updated iterator.
    constexpr derived_t & operator++() noexcept(noexcept(++base_t{}))
    //!\cond
        requires std::InputIterator<base_t>
    //!\endcond
    {
        ++(*this_to_base());
        return *this_derived();
    }

    //!\brief Post-increment, return previous iterator state.
    constexpr derived_t operator++(int) noexcept(noexcept(++base_t{}))
    //!\cond
        requires std::InputIterator<base_t>
    //!\endcond
    {
        inherited_iterator_base cpy{*this};
        ++(*this_derived());
        return static_cast<derived_t>(cpy);
    }

    //!\brief Pre-decrement, return updated iterator.
    constexpr derived_t & operator--() noexcept(noexcept(--base_t{}))
    //!\cond
        requires std::BidirectionalIterator<base_t>
    //!\endcond
    {
        --(*this_to_base());
        return *this_derived();
    }

    //!\brief Post-decrement, return previous iterator state.
    constexpr derived_t operator--(int) noexcept(noexcept(--base_t{}))
    //!\cond
        requires std::BidirectionalIterator<base_t>
    //!\endcond
    {
        inherited_iterator_base cpy{*this};
        --(*this);
        return static_cast<derived_t>(cpy);
    }

    //!\brief Move iterator to the right.
    constexpr derived_t & operator+=(difference_type const skip) noexcept(noexcept(base_t{} += skip))
    //!\cond
        requires std::RandomAccessIterator<base_t>
    //!\endcond
    {
        *this_to_base() += skip;
        return *this_derived();
    }

    //!\brief Return a an iterator moved to the right.
    constexpr derived_t operator+(difference_type const skip) const noexcept(noexcept(base_t{} += skip))
    //!\cond
        requires std::RandomAccessIterator<base_t>
    //!\endcond
    {
        derived_t cpy{*this_derived()};
        return cpy += skip;
    }

    //!\brief Non-member operator+ delegates to non-friend operator+.
    constexpr friend derived_t operator+(difference_type const skip, derived_t const & it) noexcept(noexcept(base_t{} += skip))
    //!\cond
        requires std::RandomAccessIterator<base_t>
    //!\endcond
    {
        return it + skip;
    }

    //!\brief Decrement iterator by skip.
    constexpr derived_t & operator-=(difference_type const skip) noexcept(noexcept(base_t{} -= skip))
    //!\cond
        requires std::RandomAccessIterator<base_t>
    //!\endcond
    {
        *this_to_base() -= skip;
        return *this_derived();
    }

    //!\brief Return decremented copy of this iterator.
    constexpr derived_t operator-(difference_type const skip) const noexcept(noexcept(base_t{} -= skip))
    //!\cond
        requires std::RandomAccessIterator<base_t>
    //!\endcond
    {
        derived_t cpy{*this_derived()};
        return cpy -= skip;
    }

    //!\brief Non-member operator- delegates to non-friend operator-.
    constexpr friend derived_t operator-(difference_type const skip, derived_t const & it) noexcept(noexcept(base_t{} -= skip))
    //!\cond
        requires std::RandomAccessIterator<base_t>
    //!\endcond
    {
        return it - skip;
    }

    //!\brief Return offset between this and remote iterator's position.
    constexpr difference_type operator-(derived_t const rhs) const noexcept
    //!\cond
        requires std::RandomAccessIterator<base_t>
    //!\endcond
    {
        assert(static_cast<base_t>(rhs) > *this_to_base());
        return static_cast<difference_type>(*this_to_base() - static_cast<base_t>(rhs));
    }
    //!\}

    /*!\name Reference/Dereference operators
     * \{
    */
    //!\brief Dereference operator returns element currently pointed at.
    constexpr reference operator*() const noexcept(noexcept(*base_t{}))
    //!\cond
        requires std::InputIterator<base_t>
    //!\endcond
    {
        return **this_to_base();
    }

    //!\brief Return pointer to this iterator.
    constexpr pointer operator->() const noexcept(noexcept(*base_t{}))
    //!\cond
        requires std::InputIterator<base_t>
    //!\endcond
    {
        return &*this_to_base();
    }

    //!\brief Return underlying container value currently pointed at.
    constexpr reference operator[](std::make_signed_t<difference_type> const n) const noexcept(noexcept(base_t{}[0]))
    //!\cond
        requires std::RandomAccessIterator<base_t>
    //!\endcond
    {
        return this_to_base()[n];
    }
    //!\}

private:

    //!\brief Cast this to derived type.
    derived_t * this_derived()
    {
        return static_cast<derived_t*>(this);
    }

    //!\copydoc this_derived
    derived_t const * this_derived() const
    {
        return static_cast<derived_t const *>(this);
    }

    //!\brief Cast this to base type.
    base_t * this_to_base()
    {
        return static_cast<base_t*>(this);
    }

    //!\copydoc this_to_base
    base_t const * this_to_base() const
    {
        return static_cast<base_t const *>(this);
    }
};

} // namespace seqan3::detail
