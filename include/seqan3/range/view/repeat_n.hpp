// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Provides seqan3::view::repeat_n.
 */

#pragma once

#include <seqan3/range/view/repeat.hpp>
#include <seqan3/range/view/take_exactly.hpp>
#include <seqan3/std/concepts>

namespace seqan3::detail
{

/*!\brief The underlying type of seqan3::view::repeat_n.
 * \ingroup view
 *
 * Under the hood this delegates to `view::repeat(value) | view::take_exactly(count)`.
 */
struct repeat_n_fn
{
    /*!\brief Creates a range of size `count`, where each element equals `value`.
     * \tparam    value_t The type of value to repeat; must be std::CopyConstructible.
     * \param[in] value   The value to repeat.
     * \param[in] count   The number of times to repeat `value`.
     * \returns A range of size `count`, where each element equals `value`.
     */
    template <typename value_t>
    constexpr auto operator()(value_t && value, size_t const count) const
    {
        static_assert(std::CopyConstructible<value_t>, "The value passed to repeat_n must be copy constructible.");

        return view::repeat(std::forward<value_t>(value)) | view::take_exactly(count);
    }
};

} // namespace seqan3::detail

namespace seqan3::view
{

/*!\name General purpose views
 * \{
 */
/*!\brief A view factory that repeats a given value `n` times.
 * \tparam    value_t The type of value to repeat; must be std::CopyConstructible.
 * \param[in] value   The value to repeat.
 * \param[in] count   The number of times to repeat `value`.
 * \returns A range of size `count`, where each element equals `value`.
 * \ingroup view
 *
 * \details
 *
 * **Header**
 * ```cpp
 *      #include <seqan3/range/view/repeat_n.hpp>
 * ```
 *
 * ### View properties
 *
 * This view is **source-only**, it can only be at the beginning of a pipe of range transformations.
 *
 * | range concepts and reference_t  | `rrng_t` (returned range type)                     |
 * |---------------------------------|:--------------------------------------------------:|
 * | std::ranges::InputRange         | *guaranteed*                                       |
 * | std::ranges::ForwardRange       | *guaranteed*                                       |
 * | std::ranges::BidirectionalRange | *guaranteed*                                       |
 * | std::ranges::RandomAccessRange  | *guaranteed*                                       |
 * | std::ranges::ContiguousRange    |                                                    |
 * |                                 |                                                    |
 * | std::ranges::ViewableRange      | *guaranteed*                                       |
 * | std::ranges::View               | *guaranteed*                                       |
 * | std::ranges::SizedRange         | *guaranteed*                                       |
 * | std::ranges::CommonRange        |                                                    |
 * | std::ranges::OutputRange        | *guaranteed*                                       |
 * | seqan3::ConstIterableRange      | *guaranteed*                                       |
 * |                                 |                                                    |
 * | seqan3::reference_t             | std::remove_reference_t<value_t> &                 |
 *
 * See the \link view view submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * \attention The given value to repeat is always **copied** into the range.
 *
 * ### Example
 *
 * \include test/snippet/range/view/repeat_n.cpp
 *
 * \hideinitializer
 */
constexpr inline auto repeat_n = detail::repeat_n_fn{};
//!\}

} // namespace seqan3::view
