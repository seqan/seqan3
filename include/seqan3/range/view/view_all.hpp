// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::view::all.
 */

#pragma once

#include <string_view>

#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/view/detail.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>
#include <seqan3/std/span>

namespace seqan3::detail
{

// ============================================================================
//  all_fn (adaptor definition)
// ============================================================================

/*!\brief View adaptor definition for view::all.
 */
class all_fn : public adaptor_base<all_fn>
{
private:
    //!\brief Type of the CRTP-base.
    using base_t = adaptor_base<all_fn>;

public:
    //!\brief Inherit the base class's Constructors.
    using base_t::base_t;

private:
    //!\brief Befriend the base class so it can call impl().
    friend base_t;

    /*!\brief Type erase if possible and delegate to std::view::all otherwise.
     * \returns An instance of std::span, std::basic_string_view, std::ranges::subrange or std::view::all's return type.
     */
    template <std::ranges::Range urng_t>
    static constexpr auto impl(urng_t && urange)
    {
        static_assert(std::ranges::ViewableRange<urng_t>,
                      "The view::all adaptor can only be passed ViewableRanges, i.e. Views or &-to-non-View.");

        // views are always passed as-is
        if constexpr (std::ranges::View<remove_cvref_t<urng_t>>)
        {
            return std::view::all(std::forward<urng_t>(urange));
        }
        // string const &
        else if constexpr (is_type_specialisation_of_v<remove_cvref_t<urng_t>, std::basic_string> &&
                           std::is_const_v<std::remove_reference_t<urng_t>>)
        {
            return std::basic_string_view{std::ranges::data(urange), std::ranges::size(urange)};
        }
        // contiguous
        else if constexpr (ForwardingRange<urng_t> &&
                           std::ranges::ContiguousRange<urng_t> &&
                           std::ranges::SizedRange<urng_t>)
        {
            return std::span{std::ranges::data(urange), std::ranges::size(urange)};
        }
        // random_access
        else if constexpr (ForwardingRange<urng_t> &&
                           std::ranges::RandomAccessRange<urng_t> &&
                           std::ranges::SizedRange<urng_t>)
        {
            return std::ranges::subrange<std::ranges::iterator_t<urng_t>, std::ranges::iterator_t<urng_t>>
            {
                std::ranges::begin(urange),
                std::ranges::begin(urange) + std::ranges::size(urange),
                std::ranges::size(urange)
            };
        }
        // pass to std::view::all (will return ref-view)
        else
        {
            return std::view::all(std::forward<urng_t>(urange));
        }
    }
};

} // namespace seqan3::detail

// ============================================================================
//  view::all (adaptor instance definition)
// ============================================================================

namespace seqan3::view
{

/*!\name General purpose views
 * \{
 */

/*!\brief               A view adaptor that behaves like std::view:all, but type erases contiguous ranges.
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \returns             Either `std::basic_string_view{urange}`, or `std::span{urange}`, or `std::view::all(urange)`.
 * \ingroup view
 *
 * \details
 *
 * **Header**
 * ```cpp
 *      #include <seqan3/range/view/view_all.hpp>
 * ```
 *
 * ### View properties
 *
 * | range concepts and reference_t  | `urng_t` (underlying range type)  | `rrng_t` (returned range type)  |
 * |---------------------------------|:---------------------------------:|:-------------------------------:|
 * | std::ranges::InputRange         | *required*                        | *preserved*                     |
 * | std::ranges::ForwardRange       |                                   | *preserved*                     |
 * | std::ranges::BidirectionalRange |                                   | *preserved*                     |
 * | std::ranges::RandomAccessRange  |                                   | *preserved*                     |
 * | std::ranges::ContiguousRange    |                                   | *preserved*                     |
 * |                                 |                                   |                                 |
 * | std::ranges::ViewableRange      | *required*                        | *guaranteed*                    |
 * | std::ranges::View               |                                   | *guaranteed*                    |
 * | std::ranges::SizedRange         |                                   | *preserved*                     |
 * | std::ranges::CommonRange        |                                   | *preserved*                     |
 * | std::ranges::OutputRange        |                                   | *preserved*                     |
 * | seqan3::ConstIterableRange      |                                   | *preserved*                     |
 * |                                 |                                   |                                 |
 * | seqan3::reference_t             |                                   | seqan3::reference_t<urng_t>     |
 *
 * See the \link view view submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Return type
 *
 * | `urng_t` (underlying range type)                                                       | `rrng_t` (returned range type)  |
 * |:--------------------------------------------------------------------------------------:|:-------------------------------:|
 * | `std::basic_string const &` *or* `std::basic_string_view`                              | `std::basic_string_view`        |
 * | `seqan3::ForwardingRange && std::ranges::SizedRange && std::ranges::ContiguousRange`   | `std::span`                     |
 * | `seqan3::ForwardingRange && std::ranges::SizedRange && std::ranges::RandomAccessRange` | `std::ranges::subrange`         |
 * | *else*                                                                                 | *implementation defined type*   |
 *
 * This adaptor is different from std::view::take in that it performs type erasure for some underlying ranges.
 * It returns exactly the type specified above.
 *
 * ### Example
 *
* \include test/snippet/range/view/view_all.cpp
 *
 * \hideinitializer
 */
inline constexpr auto all = detail::all_fn{};

//!\}

} // namespace seqan3::view

namespace seqan3
{
//!\brief Deduces the return value of seqan3::view::all.
template <typename t>
using all_view = decltype(view::all(std::declval<t>()));
}
