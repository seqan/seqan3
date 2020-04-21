// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::views::drop.
 */

#pragma once

#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/views/detail.hpp>
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

/*!\brief View adaptor definition for views::drop and views::drop_or_throw.
 */
struct drop_fn
{
    //!\brief Store the argument and return a range adaptor closure object.
    constexpr auto operator()(size_t drop_size) const noexcept
    {
        return detail::adaptor_from_functor{*this, drop_size};
    }

     /*!\brief Type erase if possible and forward to std::views::drop if not.
     * \returns An instance of std::span, std::basic_string_view, std::ranges::subrange or std::ranges::drop_view.
     */
    template <std::ranges::range urng_t>
    constexpr auto operator()(urng_t && urange, size_t drop_size) const
    {
        static_assert(std::ranges::viewable_range<urng_t>,
                      "The views::drop adaptor can only be passed viewable_ranges, i.e. Views or &-to-non-View.");

        [[maybe_unused]] size_t new_size = -1;

        // safeguard against wrong size
        if constexpr (std::ranges::sized_range<urng_t>)
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
        else if constexpr (forwarding_range<urng_t> &&
                           std::ranges::contiguous_range<urng_t> &&
                           std::ranges::sized_range<urng_t>)
        {
            return std::span{std::ranges::data(urange) + drop_size, new_size};
        }
        // random_access
        else if constexpr (forwarding_range<urng_t> &&
                           std::ranges::random_access_range<urng_t> &&
                           std::ranges::sized_range<urng_t>)
        {
            return std::ranges::subrange<std::ranges::iterator_t<urng_t>, std::ranges::iterator_t<urng_t>>
            {
                std::ranges::begin(urange) + drop_size,
                std::ranges::begin(urange) + drop_size + new_size,
                new_size
            };
        }
        // std::views::drop
        else
        {
            return std::forward<urng_t>(urange) | std::views::drop(drop_size);
        }
    }
};

} // namespace seqan3::detail

// ============================================================================
//  views::drop (adaptor instance definition)
// ============================================================================

namespace seqan3::views
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
 * \ingroup views
 *
 * \details
 *
 * \header_file{seqan3/range/views/drop.hpp}
 *
 * ### View properties
 *
 * | Concepts and traits              | `urng_t` (underlying range type)   | `rrng_t` (returned range type)         |
 * |----------------------------------|:----------------------------------:|:--------------------------------------:|
 * | std::ranges::input_range         | *required*                         | *preserved*                            |
 * | std::ranges::forward_range       |                                    | *preserved*                            |
 * | std::ranges::bidirectional_range |                                    | *preserved*                            |
 * | std::ranges::random_access_range |                                    | *preserved*                            |
 * | std::ranges::contiguous_range    |                                    | *preserved*                            |
 * |                                  |                                    |                                        |
 * | std::ranges::viewable_range      | *required*                         | *guaranteed*                           |
 * | std::ranges::view                |                                    | *guaranteed*                           |
 * | std::ranges::sized_range         |                                    | *preserved*                            |
 * | std::ranges::common_range        |                                    | *preserved*                            |
 * | std::ranges::output_range        |                                    | *preserved*                            |
 * | seqan3::const_iterable_range     |                                    | *preserved*                            |
 * |                                  |                                    |                                        |
 * | std::ranges::range_reference_t   |                                    | std::ranges::range_reference_t<urng_t> |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ### Return type
 *
 * | `urng_t` (underlying range type)                                                           | `rrng_t` (returned range type)  |
 * |:------------------------------------------------------------------------------------------:|:-------------------------------:|
 * | `std::basic_string const &` *or* `std::basic_string_view`                                  | `std::basic_string_view`        |
 * | `seqan3::forwarding_range && std::ranges::sized_range && std::ranges::contiguous_range`    | `std::span`                     |
 * | `seqan3::forwarding_range && std::ranges::sized_range && std::ranges::random_access_range` | `std::ranges::subrange`         |
 * | *else*                                                                                     | *implementation defined type*   |
 *
 * The adaptor is different from std::views::drop in that it performs type erasure for some underlying ranges.
 * It returns exactly the type specified above.
 *
 * ### Complexity
 *
 * Construction time of the returned view is in \f$ O(1) \f$ if the underlying range models at least
 * std::ranges::random_access_range and std::ranges::sized_range; otherwise in \f$ O(drop\_size) \f$.
 *
 * ### Example
 *
 * \include test/snippet/range/views/drop.cpp
 *
 * \hideinitializer
 */
inline constexpr auto drop = detail::drop_fn{};

//!\}

} // namespace seqan3::views
