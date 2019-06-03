// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::view::slice.
 */

#pragma once

#include <stdexcept>

#include <seqan3/io/exception.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/view/detail.hpp>
#include <seqan3/range/view/drop.hpp>
#include <seqan3/range/view/take.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>
#include <seqan3/std/span>
#include <seqan3/std/type_traits>

namespace seqan3::detail
{

// ============================================================================
//  slice_fn (adaptor definition)
// ============================================================================

//!\brief View adaptor definition for view::slice.
struct slice_fn
{
    //!\brief Store the arguments and return a range adaptor closure object.
    constexpr auto operator()(ptrdiff_t begin_pos, ptrdiff_t end_pos) const noexcept
    {
        return detail::adaptor_from_functor{*this, begin_pos, end_pos};
    }

    /*!\brief       Call the view's constructor with the underlying view as argument.
     * \returns     An instance of seqan3::detail::view_slice.
     */
    template <std::ranges::ViewableRange urng_t>
    constexpr auto operator()(urng_t && urange, ptrdiff_t begin_pos, ptrdiff_t end_pos) const
    {
        if constexpr (std::ranges::SizedRange<urng_t>)
        {
            begin_pos = std::min(begin_pos, static_cast<ptrdiff_t>(std::ranges::size(urange)));
            end_pos = std::min(end_pos, static_cast<ptrdiff_t>(std::ranges::size(urange)));
        }

        if (end_pos < begin_pos)
            throw std::invalid_argument{"end_pos argument to seqan3::view::slice must be >= the begin_pos argument."};

        return std::forward<urng_t>(urange) | view::drop(begin_pos) | view::take(end_pos - begin_pos);
    }

    // does not require special overloads, because view::drop and view::take handle the flattening.
};

} // namespace seqan3::detail

// ============================================================================
//  view::slice (adaptor instance definition)
// ============================================================================

namespace seqan3::view
{

/*!\name General purpose views
 * \{
 */

/*!\brief               A view adaptor that returns a half-open interval on the underlying range.
 * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
 *                      omitted in pipe notation]
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \param[in] begin_pos The beginning of the interval (index of first element returned).
 * \param[in] end_pos   The end of the interval (index **behind** the last element returned).
 * \returns             Up to `end_pos - begin_pos` elements of the underlying range.
 * \throws std::invalid_argument If `end_pos < begin_pos`.
 * \ingroup view
 *
 * \details
 *
 * **Header**
 * ```cpp
 *      #include <seqan3/range/view/slice.hpp>
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
 * This adaptor is a combination of seqan3::view::drop and seqan3::view::take.
 *
 * If `begin_pos` is larger than the size of the underlying range an empty range is returned.
 * If `end_pos` is larger than the size of the underlying range less elements are returned.
 *
 * If `end_pos < begin_pos` an exception of type std::invalid_argument is thrown.
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
 * The adaptor returns exactly the type specified above.
 *
 * ### Complexity
 *
 * Construction of the returned view is in \f$ O(begin\_pos) \f$ for some views, see seqan3::view::drop.
 *
 * ### Example
 *
 * \include test/snippet/range/view/slice.cpp
 *
 * \hideinitializer
 */
inline constexpr auto slice = detail::slice_fn{};

//!\}

} // namespace seqan3::view
