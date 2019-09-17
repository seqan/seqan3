// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::views::pseudo_random_access.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/range/views/detail.hpp>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

/*!\interface seqan3::detail::pseudo_random_access_range <>
 * \extends   std::ranges::range
 * \brief     This concept checks if a type models a random access range via `begin_ra` and `end_ra` member
 *            functions, albeit constant access time might not be guaranteed.
 * \ingroup   range
 *
 * \details
 *
 * This concept describes the requirements a range must fulfil in order to be used as a random access range.
 *
 * ### Concepts and doxygen
 *
 * The requirements for this concept are given as related functions and type traits.
 * Types that model this concept are shown as "implementing this interface".
 */
/*!\name Requirements for seqan3::detail::pseudo_random_access_range
 * \brief You can expect these functions on all types that model seqan3::pseudo_random_access_range.
 * \relates seqan3::detail::pseudo_random_access_range
 * \{
 */
/*!\fn pseudo_random_access_range::begin_ra()
 * \brief Returns a pseudo random access iterator.
 *
 * \details
 *
 * The returned iterator must model the std::random_access_iterator, although accessing elements through this iterator
 * might not be constant.
 *
 * \attention This is a concept requirement, not an actual function (however types
 *            modelling this concept will provide an implementation).
 */
/*!\fn pseudo_random_access_range::end_ra()
 * \brief Returns a pseudo random access sentinel.
 *
 * \details
 *
 * The returned sentinel must model std::sentinel_for with the corresponding iterator obtained by begin_ra().
 *
 * \attention This is a concept requirement, not an actual function (however types
 *            modelling this concept will provide an implementation).
 */
//!\cond
template <typename rng_t>
SEQAN3_CONCEPT pseudo_random_access_range = std::ranges::range<rng_t> && requires (rng_t rng)
{
    rng.begin_ra();
    rng.end_ra();

    requires std::random_access_iterator<decltype(rng.begin_ra())>;
    requires std::sentinel_for<decltype(rng.begin_ra()), decltype(rng.end_ra())>;
};
//!\endcond

// ============================================================================
//  pseudo_random_access_fn (adaptor definition)
// ============================================================================

//!\brief View adaptor definition for seqan3::views::pseudo_random_access.
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

        if constexpr (std::ranges::random_access_range<urng_t>)
        { // Nothing to do, just return as ref_view or original view.
            return std::views::all(std::forward<urng_t>(urange));
        }
        else
        { // Get a subrange using the random access iterators of the container.
            using iterator_ra_t = decltype(std::declval<urng_t>().begin_ra());
            using sentinel_ra_t = decltype(std::declval<urng_t>().end_ra());

            return std::ranges::subrange<iterator_ra_t, sentinel_ra_t>{urange.begin_ra(), urange.end_ra()};
        }
    }
};

}  // namespace seqan3::detail

namespace seqan3::views
{

/*!\name General purpose views
 * \{
 */

/*!\brief            A view adaptor that converts a pseudo random access range to a std::random_access_range.
 * \tparam urng_t    The type of the range being processed. See below for requirements. [template parameter is
 *                   omitted in pipe notation]
 * \param[in] urange The range being processed. [parameter is omitted in pipe notation]
 * \returns          A std::ranges::random_access_range over the given range.
 * \ingroup views
 *
 * \details
 *
 * **Header**
 * ```cpp
 *      #include <seqan3/range/views/pseudo_random_access.hpp>
 * ```
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
 * | std::ranges::forward_range       |                                   | *guaranteed*                               |
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
 * | seqan3::reference_t              |                                   | seqan3::reference_t<urng_t>                |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * This adaptor requires that the underlying range models either std::ranges::random_access_range or is one of the
 * following range types:
 *
 *  * seqan3::gap_decorator.
 *
 * ### Return type
 *
 * | `urng_t` (underlying range type)       | `rrng_t` (returned range type)                       |
 * |:--------------------------------------:|:----------------------------------------------------:|
 * | `std::ranges::random_access_range`     | `std::ranges::ref_view<urng_t>`                      |
 * | `seqan3::gap_decorator`                | `std::ranges::subrange<ra_iterator_t, ra_iterator_t>`|
 *
 * The adaptor returns exactly the type specified above. In the second case `ra_iterator_t` refers to
 * seqan3::gap_decorator::gap_decorator_iterator_ra.
 *
 * ### Example
 *
 * \include test/snippet/range/view/pseudo_random_access.cpp
 *
 * ### Complexity
 *
 * Construction of the returned view is in \f$ O(1) \f$.
 *
 * \hideinitializer
 */
inline constexpr auto pseudo_random_access = detail::pseudo_random_access_fn{};
//!\}
} // namespace seqan3::views
