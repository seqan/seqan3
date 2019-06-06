// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::view::drop.
 */

#pragma once

#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/view/detail.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>
#include <seqan3/std/span>
#include <seqan3/std/type_traits>

namespace seqan3::detail
{

// ============================================================================
//  drop_fn (adaptor definition)
// ============================================================================

/*!\brief View adaptor definition for view::drop and view::drop_or_throw.
 */
struct drop_fn
{
    //!\brief Store the argument and return a range adaptor closure object.
    constexpr auto operator()(size_t drop_size) const noexcept
    {
        return detail::adaptor_from_functor{*this, drop_size};
    }

     /*!\brief Type erase if possible and forward to std::view::drop if not.
     * \returns An instance of std::span, std::basic_string_view, std::ranges::subrange or std::ranges::drop_view.
     */
    template <std::ranges::Range urng_t>
    constexpr auto operator()(urng_t && urange, size_t drop_size) const
    {
        static_assert(std::ranges::ViewableRange<urng_t>,
                      "The view::drop adaptor can only be passed ViewableRanges, i.e. Views or &-to-non-View.");

        [[maybe_unused]] size_t new_size = -1;

        // safeguard against wrong size
        if constexpr (std::ranges::SizedRange<urng_t>)
        {
            drop_size = std::min(drop_size, static_cast<size_t>(std::ranges::size(urange)));
            new_size = std::ranges::size(urange) - drop_size;
        }

        // string_view
        if constexpr (is_type_specialisation_of_v<remove_cvref_t<urng_t>, std::basic_string_view>)
        {
            return urange.substr(drop_size);
        }
        // string const &
        else if constexpr (is_type_specialisation_of_v<remove_cvref_t<urng_t>, std::basic_string> &&
                           std::is_const_v<std::remove_reference_t<urng_t>>)
        {
            return std::basic_string_view{std::ranges::data(urange) + drop_size, new_size};
        }
        // contiguous
        else if constexpr (ForwardingRange<urng_t> &&
                           std::ranges::ContiguousRange<urng_t> &&
                           std::ranges::SizedRange<urng_t>)
        {
            return std::span{std::ranges::data(urange) + drop_size, new_size};
        }
        // random_access
        else if constexpr (ForwardingRange<urng_t> &&
                           std::ranges::RandomAccessRange<urng_t> &&
                           std::ranges::SizedRange<urng_t>)
        {
            return std::ranges::subrange<std::ranges::iterator_t<urng_t>, std::ranges::iterator_t<urng_t>>
            {
                std::ranges::begin(urange) + drop_size,
                std::ranges::begin(urange) + drop_size + new_size,
                new_size
            };
        }
        // std::view::drop
        else
        {
            return std::forward<urng_t>(urange) | std::view::drop(drop_size);
        }
    }
};

} // namespace seqan3::detail

// ============================================================================
//  view::drop (adaptor instance definition)
// ============================================================================

namespace seqan3::view
{

/*!\name General purpose views
 * \{
 */

/*!\brief               A view adaptor that returns all elements after n from the underlying range (or an empty range
 *                      if the underlying range is shorter).
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \param[in] drop_size The number of elements to drop from the beginning.
 * \returns             All elements of the underlying range after the first `drop_size`.
 * \ingroup view
 *
 * \details
 *
 * **Header**
 * ```cpp
 *      #include <seqan3/range/view/drop.hpp>
 * ```
 *
 * ### View properties
 *
 * | range concepts and reference_t  | `urng_t` (underlying range type)   | `rrng_t` (returned range type)  |
 * |---------------------------------|:----------------------------------:|:-------------------------------:|
 * | std::ranges::InputRange         | *required*                         | *preserved*                     |
 * | std::ranges::ForwardRange       |                                    | *preserved*                     |
 * | std::ranges::BidirectionalRange |                                    | *preserved*                     |
 * | std::ranges::RandomAccessRange  |                                    | *preserved*                     |
 * | std::ranges::ContiguousRange    |                                    | *preserved*                     |
 * |                                 |                                    |                                 |
 * | std::ranges::ViewableRange      | *required*                         | *guaranteed*                    |
 * | std::ranges::View               |                                    | *guaranteed*                    |
 * | std::ranges::SizedRange         |                                    | *preserved*                     |
 * | std::ranges::CommonRange        |                                    | *preserved*                     |
 * | std::ranges::OutputRange        |                                    | *preserved*                     |
 * | seqan3::ConstIterableRange      |                                    | *preserved*                     |
 * |                                 |                                    |                                 |
 * | seqan3::reference_t             |                                    | seqan3::reference_t<urng_t>     |
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
 * The adaptor is different from std::view::drop in that it performs type erasure for some underlying ranges.
 * It returns exactly the type specified above.
 *
 * ### Complexity
 *
 * Construction time of the returned view is in \f$ O(1) \f$ if the underlying range models at least
 * std::ranges::RandomAccessRange and std::ranges::SizedRange; otherwise in \f$ O(drop\_size) \f$.
 *
 * ### Example
 *
 * \include test/snippet/range/view/drop.cpp
 *
 * \hideinitializer
 */
inline constexpr auto drop = detail::drop_fn{};

//!\}

} // namespace seqan3::view
