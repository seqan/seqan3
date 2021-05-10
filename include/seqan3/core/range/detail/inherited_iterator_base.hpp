// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the seqan3::detail::inherited_iterator_base template.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <cassert>
#include <seqan3/std/iterator>
#include <type_traits>

#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/core/detail/iterator_traits.hpp>

namespace seqan3::detail
{

/*!\brief A CRTP base template for creating iterators that inherit from other iterators.
 * \tparam derived_t The CRTP specialisation.
 * \tparam base_t    The type to inherit from; must satisfy std::input_or_output_iterator.
 * \implements std::input_or_output_iterator
 * \ingroup range
 *
 * \details
 *
 * This template enables you to inherit from another iterator and just overload those functions
 * that you wish to change.
 *
 * Note that many of this class's members assume that the derived type is constructible from the base type.
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
template <typename derived_t, std::input_or_output_iterator base_t>
class inherited_iterator_base : public std::conditional_t<std::is_pointer_v<base_t> || !std::semiregular<base_t>,
                                                          empty_type,
                                                          base_t>,
                                public maybe_inherited_iterator_category<base_t>
{
private:
    //!\brief Whether this iterator inherits or wraps.
    static constexpr bool wrap_base = std::is_pointer_v<base_t> || !std::semiregular<base_t>;
public:
    /*!\name Associated types
     * \brief All are derived from the base_t.
     * \{
     */

    //!\brief The difference type.
    using difference_type = std::iter_difference_t<base_t>;
    //!\brief The value type.
    using value_type = std::iter_value_t<base_t>;
    //!\brief The reference type.
    using reference = std::iter_reference_t<base_t>;
    //!\brief The pointer type.
    using pointer = detail::iter_pointer_t<base_t>;
#if SEQAN3_DOXYGEN_ONLY(1)0
    //!\brief The iterator category tag.
    using iterator_category = maybe_present;
#endif // SEQAN3_DOXYGEN_ONLY(1)0
    //!\brief The iterator concept tag.
    using iterator_concept = detail::iterator_concept_tag_t<base_t>;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \brief The exception specification is explicitly "inherited" to also work for pointers as base.
     * \{
     */
    constexpr inherited_iterator_base()
        noexcept(std::is_nothrow_default_constructible_v<base_t>)                       = default; //!< Defaulted.
    constexpr inherited_iterator_base(inherited_iterator_base const & rhs)
        noexcept(std::is_nothrow_copy_constructible_v<base_t>)                          = default; //!< Defaulted.
    constexpr inherited_iterator_base(inherited_iterator_base && rhs)
        noexcept(std::is_nothrow_move_constructible_v<base_t>)                          = default; //!< Defaulted.
    constexpr inherited_iterator_base & operator=(inherited_iterator_base const & rhs)
        noexcept(std::is_nothrow_copy_assignable_v<base_t>)                             = default; //!< Defaulted.
    constexpr inherited_iterator_base & operator=(inherited_iterator_base && rhs)
        noexcept(std::is_nothrow_move_assignable_v<base_t>)                             = default; //!< Defaulted.
    ~inherited_iterator_base()
        noexcept(std::is_nothrow_destructible_v<base_t>)                                = default; //!< Defaulted.

    //!\brief Delegate to base class if inheriting from non-pointer iterator.
    constexpr inherited_iterator_base(base_t it) noexcept(std::is_nothrow_move_constructible_v<base_t>)
    //!\cond
        requires (!wrap_base)
    //!\endcond
        : base_t{std::move(it)}
    {}

    //!\brief Initialise member if deriving from pointer.
    constexpr inherited_iterator_base(base_t it) noexcept
    //!\cond
        requires wrap_base
    //!\endcond
        : member{std::move(it)}
    {}
    //!\}

    //!\brief Get a copy of the base.
    constexpr base_t base() const &
    //!\cond
        requires std::copyable<base_t>
    //!\endcond
    {
        return as_base();
    }

    //!\brief Returns an [rvalue](https://en.cppreference.com/w/cpp/language/value_category) of the base.
    constexpr base_t base() &&
    {
        return std::move(as_base());
    }

    /*!\name Comparison operators
     * \brief Unless specialised in derived_type, all operators perform base_t's operator and cast to derived_t.
     * \{
     */

    //!\brief Checks whether `*this` is equal to `rhs`.
    constexpr bool operator==(derived_t const & rhs) const
        noexcept(noexcept(std::declval<base_t &>() == std::declval<base_t &>()))
    //!\cond
        requires std::equality_comparable<base_t>
    //!\endcond
    {
        return base() == rhs.base();
    }

    //!\brief Checks whether `*this` is not equal to `rhs`.
    constexpr bool operator!=(derived_t const & rhs) const
        noexcept(noexcept(std::declval<base_t &>() == std::declval<base_t &>()))
    //!\cond
        requires std::equality_comparable<base_t>
    //!\endcond
    {
        return !(*this == rhs);
    }

    //!\brief Checks whether `*this` is less than `rhs`.
    constexpr bool operator<(derived_t const & rhs) const
        noexcept(noexcept(std::declval<base_t &>() < std::declval<base_t &>()))
    //!\cond
        requires std::totally_ordered<base_t>
    //!\endcond
    {
        return base() < rhs.base();
    }

    //!\brief Checks whether `*this` is greater than `rhs`.
    constexpr bool operator>(derived_t const & rhs) const
        noexcept(noexcept(std::declval<base_t &>() > std::declval<base_t &>()))
    //!\cond
        requires std::totally_ordered<base_t>
    //!\endcond
    {
        return base() > rhs.base();
    }

    //!\brief Checks whether `*this` is less than or equal to `rhs`.
    constexpr bool operator<=(derived_t const & rhs) const
        noexcept(noexcept(std::declval<base_t &>() > std::declval<base_t &>()))
    //!\cond
        requires std::totally_ordered<base_t>
    //!\endcond
    {
        return !(*this > rhs);
    }

    //!\brief Checks whether `*this` is greater than or equal to `rhs`.
    constexpr bool operator>=(derived_t const & rhs) const
        noexcept(noexcept(std::declval<base_t &>() < std::declval<base_t &>()))
    //!\cond
        requires std::totally_ordered<base_t>
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
    //!\cond
    template <typename base_t_ = base_t>
    //!\endcond
    constexpr derived_t & operator++() noexcept(noexcept(++std::declval<base_t &>()))
    //!\cond
        requires requires (base_t_ i) { ++i; }
    //!\endcond
    {
        ++as_base();
        return *this_derived();
    }

    //!\brief Post-increment of input iterators that do not return a copy of themselves but `void` or a proxy type.
    //!\cond
    template <typename base_t_ = base_t>
    //!\endcond
    constexpr auto operator++(int) noexcept(noexcept(std::declval<base_t &>()++))
    //!\cond
        requires requires (base_t_ i) { i++; requires !std::same_as<decltype(i++), base_t_>; }
    //!\endcond
    {
        return as_base()++;
    }

    //!\brief Post-increment, return previous iterator state.
    //!\cond
    template <typename base_t_ = base_t>
    //!\endcond
    constexpr derived_t operator++(int) noexcept(noexcept(std::declval<base_t &>()++) &&
                                                 noexcept(derived_t(std::declval<base_t &>())))
    //!\cond
        requires requires (base_t_ i) { i++; SEQAN3_RETURN_TYPE_CONSTRAINT(i++, std::same_as, base_t_); } &&
                 std::constructible_from<derived_t, base_t_>
    //!\endcond
    {
        return derived_t{as_base()++};
    }

    //!\brief Pre-decrement, return updated iterator.
    //!\cond
    template <typename base_t_ = base_t>
    //!\endcond
    constexpr derived_t & operator--() noexcept(noexcept(--std::declval<base_t &>()))
    //!\cond
        requires requires (base_t_ i) { --i; }
    //!\endcond
    {
        --as_base();
        return *this_derived();
    }

    //!\brief Post-decrement, return previous iterator state.
    //!\cond
    template <typename base_t_ = base_t>
    //!\endcond
    constexpr derived_t operator--(int) noexcept(noexcept(std::declval<base_t &>()--) &&
                                                 noexcept(derived_t{std::declval<base_t &>()}))
    //!\cond
        requires requires (base_t_ i) { i--; } && std::constructible_from<derived_t, base_t_>
    //!\endcond
    {
        return derived_t{as_base()--};
    }

    //!\brief Move iterator to the right.
    //!\cond
    template <typename base_t_ = base_t>
    //!\endcond
    constexpr derived_t & operator+=(difference_type const skip) noexcept(noexcept(std::declval<base_t &>() += skip))
    //!\cond
        requires requires (base_t_ i, difference_type const n) { i += n; }
    //!\endcond
    {
        as_base() += skip;
        return *this_derived();
    }

    //!\brief Returns an iterator which is advanced by `skip` positions.
    //!\cond
    template <typename base_t_ = base_t>
    //!\endcond
    constexpr derived_t operator+(difference_type const skip) const
        noexcept(noexcept(std::declval<base_t &>() + skip) && noexcept(derived_t{std::declval<base_t &>()}))
    //!\cond
        requires requires (base_t_ const i, difference_type const n) { i + n; } &&
                 std::constructible_from<derived_t, base_t_>
    //!\endcond
    {
        return derived_t{as_base() + skip};
    }

    //!\brief Non-member operator+ delegates to non-friend operator+.
    constexpr friend derived_t operator+(difference_type const skip, derived_t const & it)
        noexcept(noexcept(skip + std::declval<base_t const &>()))
    //!\cond
        requires requires (base_t const i, difference_type const n) { n + i; } &&
                 std::constructible_from<derived_t, base_t>
    //!\endcond
    {
        return it + skip;
    }

    //!\brief Decrement iterator by skip.
    //!\cond
    template <typename base_t_ = base_t>
    //!\endcond
    constexpr derived_t & operator-=(difference_type const skip) noexcept(noexcept(std::declval<base_t &>() -= skip))
    //!\cond
        requires requires (base_t_ i, difference_type const n) { i -= n; }
    //!\endcond
    {
        as_base() -= skip;
        return *this_derived();
    }

    //!\brief Return decremented copy of this iterator.
    //!\cond
    template <typename base_t_ = base_t>
    //!\endcond
    constexpr derived_t operator-(difference_type const skip) const
        noexcept(noexcept(std::declval<base_t const &>() - skip) && noexcept(derived_t(std::declval<base_t &>())))
    //!\cond
        requires requires (base_t_ i, difference_type const n) { i - n; } &&
                 std::constructible_from<derived_t, base_t_>
    //!\endcond
    {
        return derived_t{as_base() - skip};
    }

    //!\brief Return offset between this and remote iterator's position.
    constexpr difference_type operator-(derived_t const & rhs) const
        noexcept(noexcept(std::declval<base_t &>() - std::declval<base_t &>()))
    //!\cond
        requires std::sized_sentinel_for<base_t, base_t>
    //!\endcond
    {
        return as_base() - rhs.as_base();
    }
    //!\}

    /*!\name Reference/Dereference operators
     * \{
     */
    //!\brief Dereference operator returns element currently pointed at.
    constexpr reference operator*() noexcept(noexcept(*std::declval<base_t &>()))
    //!\cond
        requires std::indirectly_readable<base_t>
    //!\endcond
    {
        return *as_base();
    }

    //!\brief Dereference operator returns element currently pointed at.
    constexpr decltype(auto) operator*() const noexcept(noexcept(*std::declval<base_t const &>()))
    //!\cond
        requires std::indirectly_readable<base_t>
    //!\endcond
    {
        return *as_base();
    }

    //!\brief Return pointer to this iterator.
    constexpr pointer operator->() noexcept(noexcept(*std::declval<base_t &>()))
    //!\cond
        requires std::input_iterator<base_t>
    //!\endcond
    {
        return &as_base();
    }

    //!\brief Return pointer to this iterator.
    constexpr decltype(auto) operator->() const noexcept(noexcept(*std::declval<base_t const &>()))
    //!\cond
        requires std::input_iterator<base_t>
    //!\endcond
    {
        return &as_base();
    }

    //!\brief Return underlying container value currently pointed at.
    //!\cond
    template <typename base_t_ = base_t>
    //!\endcond
    constexpr decltype(auto) operator[](std::make_signed_t<difference_type> const n)
        noexcept(noexcept(std::declval<base_t &>()[0]))
    //!\cond
        requires requires (base_t_ i, difference_type const n) { i[n]; }
    //!\endcond
    {
        return as_base()[n];
    }

    //!\brief Return underlying container value currently pointed at.
    //!\cond
    template <typename base_t_ = base_t>
    //!\endcond
    constexpr decltype(auto) operator[](std::make_signed_t<difference_type> const n) const
        noexcept(noexcept(std::declval<base_t const &>()[0]))
    //!\cond
        requires requires (base_t_ const i, difference_type const n) { i[n]; }
    //!\endcond
    {
        return as_base()[n];
    }
    //!\}

private:
    //!\brief If the base is a pointer, we wrap it instead of inheriting.
    std::conditional_t<wrap_base, base_t, empty_type> member;

    //!\brief Befriend the derived type so it can access the private members.
    friend derived_t;

    //!\brief Cast this to base type.
    constexpr base_t & as_base() & noexcept
    {
        if constexpr (wrap_base)
            return member;
        else
            return *this;
    }

    //!\copydoc as_base
    constexpr base_t const & as_base() const & noexcept
    {
        if constexpr (wrap_base)
            return member;
        else
            return *this;
    }

    //!\brief Cast this to derived type.
    constexpr derived_t * this_derived()
    {
        return static_cast<derived_t*>(this);
    }

    //!\copydoc this_derived
    constexpr derived_t const * this_derived() const
    {
        return static_cast<derived_t const *>(this);
    }
};

} // namespace seqan3::detail
