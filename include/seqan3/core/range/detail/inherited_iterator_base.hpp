// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides the seqan3::detail::inherited_iterator_base template.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <cassert>
#include <iterator>
#include <type_traits>

#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/core/detail/iterator_traits.hpp>

namespace seqan3::detail
{

/*!\brief A CRTP base template for creating iterators that inherit from other iterators.
 * \tparam derived_t The CRTP specialisation.
 * \tparam base_t    The type to inherit from; must satisfy std::input_or_output_iterator.
 * \implements std::input_or_output_iterator
 * \ingroup core_range
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
 * \snippet test/unit/core/range/detail/inherited_iterator_base_test.cpp inherited_iterator_base desired
 *
 * You could define it like this:
 *
 * \snippet test/unit/core/range/detail/inherited_iterator_base_test.cpp inherited_iterator_base def
 */
template <typename derived_t, std::input_or_output_iterator base_t>
class inherited_iterator_base :
    public std::conditional_t<std::is_pointer_v<base_t> || !std::semiregular<base_t>, empty_type, base_t>,
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
#if SEQAN3_DOXYGEN_ONLY(1) 0
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
        noexcept(std::is_nothrow_default_constructible_v<base_t>) = default; //!< Defaulted.
    constexpr inherited_iterator_base(inherited_iterator_base const & rhs)
        noexcept(std::is_nothrow_copy_constructible_v<base_t>) = default; //!< Defaulted.
    constexpr inherited_iterator_base(inherited_iterator_base && rhs)
        noexcept(std::is_nothrow_move_constructible_v<base_t>) = default; //!< Defaulted.
    constexpr inherited_iterator_base & operator=(inherited_iterator_base const & rhs)
        noexcept(std::is_nothrow_copy_assignable_v<base_t>) = default; //!< Defaulted.
    constexpr inherited_iterator_base & operator=(inherited_iterator_base && rhs)
        noexcept(std::is_nothrow_move_assignable_v<base_t>) = default;                     //!< Defaulted.
    ~inherited_iterator_base() noexcept(std::is_nothrow_destructible_v<base_t>) = default; //!< Defaulted.

    //!\brief Delegate to base class if inheriting from non-pointer iterator.
    constexpr inherited_iterator_base(base_t it) noexcept(std::is_nothrow_move_constructible_v<base_t>)
        requires (!wrap_base)
        : base_t{std::move(it)}
    {}

    //!\brief Initialise member if deriving from pointer.
    constexpr inherited_iterator_base(base_t it) noexcept
        requires wrap_base
        : member{std::move(it)}
    {}
    //!\}

    //!\brief Get a const reference to the base.
    constexpr base_t const & base() const & noexcept
    {
        return as_base();
    }

    //!\brief Get a reference to the base.
    constexpr base_t & base() & noexcept
    {
        return as_base();
    }

    //!\brief Returns an [rvalue](https://en.cppreference.com/w/cpp/language/value_category) of the base.
    constexpr base_t base() && noexcept
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
        requires std::equality_comparable<base_t>
    {
        return base() == rhs.base();
    }

    //!\brief Checks whether `*this` is not equal to `rhs`.
    constexpr bool operator!=(derived_t const & rhs) const
        noexcept(noexcept(std::declval<base_t &>() == std::declval<base_t &>()))
        requires std::equality_comparable<base_t>
    {
        return !(*this == rhs);
    }

    //!\brief Checks whether `*this` is less than `rhs`.
    constexpr bool operator<(derived_t const & rhs) const
        noexcept(noexcept(std::declval<base_t &>() < std::declval<base_t &>()))
        requires std::totally_ordered<base_t>
    {
        return base() < rhs.base();
    }

    //!\brief Checks whether `*this` is greater than `rhs`.
    constexpr bool operator>(derived_t const & rhs) const
        noexcept(noexcept(std::declval<base_t &>() > std::declval<base_t &>()))
        requires std::totally_ordered<base_t>
    {
        return base() > rhs.base();
    }

    //!\brief Checks whether `*this` is less than or equal to `rhs`.
    constexpr bool operator<=(derived_t const & rhs) const
        noexcept(noexcept(std::declval<base_t &>() > std::declval<base_t &>()))
        requires std::totally_ordered<base_t>
    {
        return !(*this > rhs);
    }

    //!\brief Checks whether `*this` is greater than or equal to `rhs`.
    constexpr bool operator>=(derived_t const & rhs) const
        noexcept(noexcept(std::declval<base_t &>() < std::declval<base_t &>()))
        requires std::totally_ordered<base_t>
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
        requires requires (base_t_ i) { ++i; }
    {
        ++as_base();
        return *this_derived();
    }

    //!\brief Post-increment of input iterators that do not return a copy of themselves but `void` or a proxy type.
    //!\cond
    template <typename base_t_ = base_t>
    //!\endcond
    constexpr auto operator++(int) noexcept(noexcept(std::declval<base_t &>()++))
        requires requires (base_t_ i) {
            i++;
            requires !std::same_as<decltype(i++), base_t_>;
        }
    {
        return as_base()++;
    }

    //!\brief Post-increment, return previous iterator state.
    //!\cond
    template <typename base_t_ = base_t>
    //!\endcond
    constexpr derived_t operator++(int)
        noexcept(noexcept(std::declval<base_t &>()++) && std::is_nothrow_copy_constructible_v<derived_t>)
        requires requires (base_t_ i) {
            i++;
            { i++ } -> std::same_as<base_t_>;
        } && std::copy_constructible<derived_t>
    {
        derived_t tmp = *this_derived();
        ++as_base();
        return tmp;
    }

    //!\brief Pre-decrement, return updated iterator.
    //!\cond
    template <typename base_t_ = base_t>
    //!\endcond
    constexpr derived_t & operator--() noexcept(noexcept(--std::declval<base_t &>()))
        requires requires (base_t_ i) { --i; }
    {
        --as_base();
        return *this_derived();
    }

    //!\brief Post-decrement, return previous iterator state.
    //!\cond
    template <typename base_t_ = base_t>
    //!\endcond
    constexpr derived_t operator--(int)
        noexcept(noexcept(std::declval<base_t &>()--) && std::is_nothrow_copy_constructible_v<derived_t>)
        requires requires (base_t_ i) { i--; } && std::copy_constructible<derived_t>
    {
        derived_t tmp = *this_derived();
        --as_base();
        return tmp;
    }

    //!\brief Move iterator to the right.
    //!\cond
    template <typename base_t_ = base_t>
    //!\endcond
    constexpr derived_t & operator+=(difference_type const skip) noexcept(noexcept(std::declval<base_t &>() += skip))
        requires requires (base_t_ i, difference_type const n) { i += n; }
    {
        as_base() += skip;
        return *this_derived();
    }

    //!\brief Returns an iterator which is advanced by `skip` positions.
    //!\cond
    template <typename base_t_ = base_t>
    //!\endcond
    constexpr derived_t operator+(difference_type const skip) const
        noexcept(noexcept(std::declval<base_t &>() + skip) && std::is_nothrow_copy_constructible_v<derived_t>)
        requires requires (base_t_ const i, difference_type const n) { i + n; } && std::copy_constructible<derived_t>
    {
        derived_t tmp = *this_derived();
        tmp.as_base() += skip;
        return tmp;
    }

    //!\brief Non-member operator+ delegates to non-friend operator+.
    //!\cond
    // Make this function a function template. This means that constructible_from will be evaluated with the complete
    // derived_t (which is instantiated in the second pass). If this wasn't a function template, the concept would have
    // to be evaluated in the first pass, where derived_t is incomplete.
    template <typename base_t_ = base_t>
    //!\endcond
    constexpr friend derived_t operator+(difference_type const skip, derived_t const & it)
        noexcept(noexcept(skip + std::declval<base_t const &>()))
        requires requires (base_t const i, difference_type const n) { n + i; }
              && std::constructible_from<derived_t, base_t>
    {
        return it + skip;
    }

    //!\brief Decrement iterator by skip.
    //!\cond
    template <typename base_t_ = base_t>
    //!\endcond
    constexpr derived_t & operator-=(difference_type const skip) noexcept(noexcept(std::declval<base_t &>() -= skip))
        requires requires (base_t_ i, difference_type const n) { i -= n; }
    {
        as_base() -= skip;
        return *this_derived();
    }

    //!\brief Return decremented copy of this iterator.
    //!\cond
    template <typename base_t_ = base_t>
    //!\endcond
    constexpr derived_t operator-(difference_type const skip) const
        noexcept(noexcept(std::declval<base_t const &>() - skip) && std::is_nothrow_copy_constructible_v<derived_t>)
        requires requires (base_t_ i, difference_type const n) { i - n; } && std::copy_constructible<derived_t>
    {
        derived_t tmp = *this_derived();
        tmp.as_base() -= skip;
        return tmp;
    }

    //!\brief Return offset between this and remote iterator's position.
    constexpr difference_type operator-(derived_t const & rhs) const
        noexcept(noexcept(std::declval<base_t &>() - std::declval<base_t &>()))
        requires std::sized_sentinel_for<base_t, base_t>
    {
        return as_base() - rhs.as_base();
    }
    //!\}

    /*!\name Reference/Dereference operators
     * \{
     */
    //!\brief Dereference operator returns element currently pointed at.
    constexpr reference operator*() noexcept(noexcept(*std::declval<base_t &>()))
        requires std::indirectly_readable<base_t>
    {
        return *as_base();
    }

    //!\brief Dereference operator returns element currently pointed at.
    constexpr decltype(auto) operator*() const noexcept(noexcept(*std::declval<base_t const &>()))
        requires std::indirectly_readable<base_t>
    {
        return *as_base();
    }

    //!\brief Return pointer to this iterator.
    constexpr pointer operator->() noexcept(noexcept(*std::declval<base_t &>()))
        requires std::input_iterator<base_t>
    {
        return &as_base();
    }

    //!\brief Return pointer to this iterator.
    constexpr decltype(auto) operator->() const noexcept(noexcept(*std::declval<base_t const &>()))
        requires std::input_iterator<base_t>
    {
        return &as_base();
    }

    //!\brief Return underlying container value currently pointed at.
    //!\cond
    template <typename base_t_ = base_t>
    //!\endcond
    constexpr decltype(auto) operator[](std::make_signed_t<difference_type> const n)
        noexcept(noexcept(std::declval<base_t &>()[0]))
        requires requires (base_t_ i, difference_type const n) { i[n]; }
    {
        return as_base()[n];
    }

    //!\brief Return underlying container value currently pointed at.
    //!\cond
    template <typename base_t_ = base_t>
    //!\endcond
    constexpr decltype(auto) operator[](std::make_signed_t<difference_type> const n) const
        noexcept(noexcept(std::declval<base_t const &>()[0]))
        requires requires (base_t_ const i, difference_type const n) { i[n]; }
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

    //!\copydoc seqan3::detail::inherited_iterator_base::as_base
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
        return static_cast<derived_t *>(this);
    }

    //!\copydoc seqan3::detail::inherited_iterator_base::this_derived
    constexpr derived_t const * this_derived() const
    {
        return static_cast<derived_t const *>(this);
    }
};

} // namespace seqan3::detail
