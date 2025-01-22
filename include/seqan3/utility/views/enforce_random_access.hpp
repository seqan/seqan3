// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \brief Provides seqan3::views::enforce_random_access.
 */

#pragma once

#include <iterator>
#include <ranges>
#include <type_traits>

#include <seqan3/core/range/detail/adaptor_base.hpp>
#include <seqan3/core/range/detail/inherited_iterator_base.hpp>
#include <seqan3/utility/range/concept.hpp>

namespace seqan3::detail
{

// ============================================================================
//  view_pseudo_random_access
// ============================================================================

/*!\brief View to force random access range iterator for seqan3::pseudo_random_access_range.
 * \tparam urng_t The underlying range type; must model std::ranges::view and seqan3::pseudo_random_access_range.
 * \ingroup utility_views
 *
 * \details
 *
 * Wraps the iterator of a seqan3::pseudo_random_access_range and overwrites the iterator category to be
 * std::random_access_iterator_tag. Thus, the resulting range can be used in algorithms or other contexts that require
 * random access, although the time complexity still depends on the underlying range, that is it is not guaranteed to
 * be constant.
 */
template <std::ranges::view urng_t>
    requires pseudo_random_access_range<urng_t>
class view_enforce_random_access : public std::ranges::view_interface<view_enforce_random_access<urng_t>>
{
private:
    // Iterator declaration.
    template <typename underlying_iter_t>
    class basic_iterator;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr view_enforce_random_access() = default;                                               //!< Defaulted.
    constexpr view_enforce_random_access(view_enforce_random_access const &) = default;             //!< Defaulted.
    constexpr view_enforce_random_access(view_enforce_random_access &&) = default;                  //!< Defaulted.
    constexpr view_enforce_random_access & operator=(view_enforce_random_access const &) = default; //!< Defaulted.
    constexpr view_enforce_random_access & operator=(view_enforce_random_access &&) = default;      //!< Defaulted.
    ~view_enforce_random_access() = default;                                                        //!< Defaulted.

    //!\brief Construction from the underlying view.
    explicit constexpr view_enforce_random_access(urng_t && range) : urng{std::move(range)}
    {}

    //!\brief Construction from the underlying viewable range.
    template <typename viewable_rng_t>
        requires (!std::same_as<std::remove_cvref_t<viewable_rng_t>, view_enforce_random_access>)
              && std::ranges::viewable_range<viewable_rng_t>
              && std::constructible_from<urng_t, std::ranges::ref_view<std::remove_reference_t<viewable_rng_t>>>
    explicit constexpr view_enforce_random_access(viewable_rng_t && range) :
        view_enforce_random_access{std::views::all(std::forward<viewable_rng_t>(range))}
    {}
    //!\}

    /*!\name Iterators
     * \{
     */
    //!\brief Returns the iterator to the begin of the range.
    constexpr auto begin() noexcept
    {
        return basic_iterator<decltype(std::ranges::begin(urng))>{std::ranges::begin(urng)};
    }

    //!\copydoc seqan3::detail::view_enforce_random_access::begin
    constexpr auto begin() const noexcept
        requires const_iterable_range<urng_t>
    {
        return basic_iterator<decltype(std::ranges::begin(urng))>{std::ranges::begin(urng)};
    }

    /*!\brief Returns the sentinel to the end of the range.
     *
     * \details
     *
     * If the underlying range is a common range this functions returns
     * seqan3::detail::view_enforce_random_access::basic_iterator initialised with the end of the
     * underlying range. Otherwise it returns the sentinel of the underlying range.
     */
    constexpr auto end() noexcept
    {
        if constexpr (std::ranges::common_range<urng_t>)
            return basic_iterator<decltype(std::ranges::end(urng))>{std::ranges::end(urng)};
        else
            return urng.end();
    }

    //!\copydoc seqan3::detail::view_enforce_random_access::end
    constexpr auto end() const noexcept
        requires const_iterable_range<urng_t>
    {
        if constexpr (std::ranges::common_range<urng_t>)
            return basic_iterator<decltype(std::ranges::end(urng))>{std::ranges::end(urng)};
        else
            return std::ranges::cend(urng);
    }
    //!\}

    urng_t urng; //!< The underlying range.
};

/*!\brief Iterator wrapper for the underlying range iterator enforcing std::random_access_iterator_tag.
 * \tparam underlying_iter_t The type of the underlying range iterator.
 *
 * \details
 *
 * This class inherits all properties of the underlying range iterator and overwrites the iterator category to be
 * std::random_access_range_tag.
 */
template <std::ranges::view urng_t>
    requires pseudo_random_access_range<urng_t>
template <typename underlying_iter_t>
class view_enforce_random_access<urng_t>::basic_iterator :
    public inherited_iterator_base<basic_iterator<underlying_iter_t>, underlying_iter_t>
{
private:
    //!\brief The type of the base class.
    using base_t = inherited_iterator_base<basic_iterator<underlying_iter_t>, underlying_iter_t>;

public:
    //!\brief The new iterator category.
    using iterator_category = std::random_access_iterator_tag;
    //!\brief The new iterator concept.
    using iterator_concept = iterator_category;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    // Importing base's constructors.
    using base_t::base_t;
    //!\brief Defaulted.
    constexpr basic_iterator() = default;
    //!\brief Defaulted.
    constexpr basic_iterator(basic_iterator const &) = default;
    //!\brief Defaulted.
    constexpr basic_iterator(basic_iterator &&) = default;
    //!\brief Defaulted.
    constexpr basic_iterator & operator=(basic_iterator const &) = default;
    //!\brief Defaulted.
    constexpr basic_iterator & operator=(basic_iterator &&) = default;
    //!\brief Defaulted.
    ~basic_iterator() = default;
    //!\}

    /*!\name Comparison operators
     * \brief Comparison with sentinel of underlying range.
     * \{
     */
    // Importing base's equality operators
    using base_t::operator==;
    using base_t::operator!=;
    //!\brief Tests if iterator is at the end.
    friend constexpr bool operator==(basic_iterator const & lhs, std::ranges::sentinel_t<urng_t> const & rhs)
        noexcept(noexcept(std::declval<underlying_iter_t const &>()
                          == std::declval<std::ranges::sentinel_t<urng_t> const &>()))
    {
        return lhs.base() == rhs;
    }

    //!\brief Tests if iterator is at the end.
    friend constexpr bool operator==(std::ranges::sentinel_t<urng_t> const & lhs, basic_iterator const & rhs)
        noexcept(noexcept(std::declval<underlying_iter_t const &>()
                          == std::declval<std::ranges::sentinel_t<urng_t> const &>()))
    {
        return rhs == lhs;
    }

    //!\brief Tests if iterator is not at the end.
    friend constexpr bool operator!=(basic_iterator const & lhs, std::ranges::sentinel_t<urng_t> const & rhs)
        noexcept(noexcept(std::declval<underlying_iter_t const &>()
                          != std::declval<std::ranges::sentinel_t<urng_t> const &>()))
    {
        return !(lhs == rhs);
    }

    //!\brief Tests if iterator is not at the end.
    friend constexpr bool operator!=(std::ranges::sentinel_t<urng_t> const & lhs, basic_iterator const & rhs)
        noexcept(noexcept(std::declval<underlying_iter_t const &>()
                          != std::declval<std::ranges::sentinel_t<urng_t> const &>()))
    {
        return rhs != lhs;
    }
    //!\}

    /*!\name Arithmetic operators
     * \{
     */
    // Import operator from base.
    using base_t::operator-;

    //!\brief Computes the distance betwen this iterator and the sentinel of the underlying range.
    constexpr typename base_t::difference_type operator-(std::ranges::sentinel_t<urng_t> const & rhs) const
        noexcept(noexcept(std::declval<underlying_iter_t const &>()
                          - std::declval<std::ranges::sentinel_t<urng_t> const &>()))
        requires std::sized_sentinel_for<std::ranges::sentinel_t<urng_t>, underlying_iter_t>
    {
        return this->base() - rhs;
    }

    //!\brief Computes the distance betwen this iterator and the sentinel of the underlying range.
    constexpr friend typename base_t::difference_type operator-(std::ranges::sentinel_t<urng_t> const & lhs,
                                                                basic_iterator const & rhs)
        noexcept(noexcept(std::declval<std::ranges::sentinel_t<urng_t> const &>()
                          - std::declval<underlying_iter_t const &>()))
        requires std::sized_sentinel_for<std::ranges::sentinel_t<urng_t>, underlying_iter_t>
    {
        return lhs - rhs.base();
    }
    //!\}
};

/*!\name Type deduction guides
 * \relates seqan3::detail::view_enforce_random_access
 * \{
 */
//!\brief A deduction guide for the view class template.
template <std::ranges::viewable_range rng_t>
view_enforce_random_access(rng_t &&) -> view_enforce_random_access<std::views::all_t<rng_t>>;
//!\}

// ============================================================================
//  pseudo_random_access_fn (adaptor definition)
// ============================================================================

//!\brief View adaptor definition for seqan3::views::enforce_random_access.
struct pseudo_random_access_fn : public adaptor_base<pseudo_random_access_fn>
{
private:
    //!\brief Type of the CRTP-base.
    using base_t = adaptor_base<pseudo_random_access_fn>;

public:
    //!\brief Inherit the base class's Constructors.
    using base_t::base_t;

private:
    //!\brief Befriend the base class so it can call impl().
    friend base_t;

    /*!\brief Call the view's constructor with the underlying view as argument.
     * \returns An instance of the adapted range.
     */
    template <std::ranges::viewable_range urng_t>
    static constexpr auto impl(urng_t && urange)
    {
        static_assert(std::ranges::random_access_range<urng_t> || pseudo_random_access_range<urng_t>,
                      "The adapted range must either model std::ranges::random_access_range or must be "
                      "a specific SeqAn range type that supports pseudo random access.");
        static_assert(std::ranges::forward_range<urng_t>,
                      "The underlying range must model std::ranges::forward_range.");

        if constexpr (std::ranges::random_access_range<urng_t>)
        { // Nothing to do, just return as ref_view or original view.
            return std::views::all(std::forward<urng_t>(urange));
        }
        else
        { // Get a subrange using the random access iterators of the container.
            return view_enforce_random_access{std::forward<urng_t>(urange)};
        }
    }
};

} // namespace seqan3::detail

namespace seqan3::views
{
/*!\brief            A view adaptor that converts a pseudo random access range to a std::ranges::random_access_range.
 * \tparam urng_t    The type of the range being processed. See below for requirements. [template parameter is
 *                   omitted in pipe notation]
 * \param[in] urange The range being processed. [parameter is omitted in pipe notation]
 * \returns          A std::ranges::random_access_range over the given range.
 * \ingroup utility_views
 *
 * \details
 *
 * \header_file{seqan3/utility/views/enforce_random_access.hpp}
 *
 * A pseudo random access range is a range whose iterator typically defines all the interfaces necessary to allow
 * random access, but cannot guarantee accessing an arbitrary element in constant time.
 * Thus, the highest category it can support by default is std::ranges::bidirectional_range. However, for many of these
 * pseudo random access ranges better algorithms and data structures with sub-linear runtime complexities can be used
 * (for example logarithmic time complexity). To enforce the faster behaviour of the range in a generic
 * range-based context you can use this range adaptor, which will return a range that models
 * std::ranges::random_access_range. Note, that this does not mean that the complexity of accessing an arbitrary element
 * of the adapted range improves to constant time, but merely all syntactical requirements are fulfilled including the
 * iterator tag.
 *
 * ### View properties
 *
 * | range concepts and reference_t   | `urng_t` (underlying range type)  | `rrng_t` (returned range type)             |
 * |----------------------------------|:---------------------------------:|:------------------------------------------:|
 * | std::ranges::input_range         | *required*                        | *guaranteed*                               |
 * | std::ranges::forward_range       | *required*                        | *guaranteed*                               |
 * | std::ranges::bidirectional_range |                                   | *guaranteed*                               |
 * | std::ranges::random_access_range |                                   | *guaranteed*                               |
 * | std::ranges::contiguous_range    |                                   | *preserved*                                |
 * |                                  |                                   |                                            |
 * | std::ranges::viewable_range      | *required*                        | *guaranteed*                               |
 * | std::ranges::view                |                                   | *guaranteed*                               |
 * | std::ranges::sized_range         |                                   | *preserved*                                |
 * | std::ranges::common_range        |                                   | *preserved*                                |
 * | std::ranges::output_range        |                                   | *preserved*                                |
 * | seqan3::const_iterable_range     |                                   | *preserved*                                |
 * |                                  |                                   |                                            |
 * | std::ranges::range_reference_t   |                                   | std::ranges::range_reference_t<urng_t>     |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * This adaptor requires that the underlying range models either std::ranges::random_access_range or
 * seqan3::pseudo_random_access_range.
 *
 * ### Return type
 *
 * | `urng_t` (underlying range type)       | `rrng_t` (returned range type)                       |
 * |:--------------------------------------:|:----------------------------------------------------:|
 * | `std::ranges::random_access_range`     | `std::ranges::ref_view<urng_t>`                      |
 * | `seqan3::pseudo_random_access_range`   | `seqan3::detail::view_enforce_random_access`         |
 *
 * The adaptor returns exactly the type specified above. In the second case a view is returned whose iterator wraps
 * the iterator of the underlying range and adapts all of its functionality but overwrites the
 * iterator category to be std::random_access_iterator_tag.
 *
 * ### Example
 *
 * \include test/snippet/utility/views/enforce_random_access.cpp
 *
 * ### Complexity
 *
 * Construction of the returned view is in \f$ O(1) \f$.
 *
 * \hideinitializer
 *
 * \experimentalapi{Experimental since version 3.1.}
 */
inline constexpr auto enforce_random_access = detail::pseudo_random_access_fn{};

} // namespace seqan3::views
