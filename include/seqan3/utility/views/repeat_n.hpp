// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Provides seqan3::views::repeat_n.
 */

#pragma once

#include <concepts>

#include <seqan3/io/views/detail/take_exactly_view.hpp>
#include <seqan3/utility/views/repeat.hpp>

namespace seqan3::detail
{

/*!\brief The underlying type of seqan3::views::repeat_n.
 * \ingroup utility_views
 *
 * Under the hood this delegates to `views::repeat(value) | detail::take_exactly(count)`.
 */
struct repeat_n_fn
{
    /*!\brief Creates a range of size `count`, where each element equals `value`.
     * \tparam    value_t The type of value to repeat; must be std::copy_constructible.
     * \param[in] value   The value to repeat.
     * \param[in] count   The number of times to repeat `value`.
     * \returns A range of size `count`, where each element equals `value`.
     */
    template <typename value_t>
    constexpr auto operator()(value_t && value, size_t const count) const
    {
        static_assert(std::copy_constructible<value_t>, "The value passed to repeat_n must be copy constructible.");

        return views::repeat(std::forward<value_t>(value)) | detail::take_exactly(count);
    }
};

} // namespace seqan3::detail

namespace seqan3::views
{
/*!\brief A view factory that repeats a given value `n` times.
 * \tparam    value_t The type of value to repeat; must be std::copy_constructible.
 * \param[in] value   The value to repeat.
 * \param[in] count   The number of times to repeat `value`.
 * \returns A range of size `count`, where each element equals `value`.
 * \ingroup utility_views
 *
 * \details
 *
 * \header_file{seqan3/utility/views/repeat_n.hpp}
 *
 * ### View properties
 *
 * This view is **source-only**, it can only be at the beginning of a pipe of range transformations.
 *
 * | Concepts and traits              | `rrng_t` (returned range type)                     |
 * |----------------------------------|:--------------------------------------------------:|
 * | std::ranges::input_range         | *guaranteed*                                       |
 * | std::ranges::forward_range       | *guaranteed*                                       |
 * | std::ranges::bidirectional_range | *guaranteed*                                       |
 * | std::ranges::random_access_range | *guaranteed*                                       |
 * | std::ranges::contiguous_range    |                                                    |
 * |                                  |                                                    |
 * | std::ranges::viewable_range      | *guaranteed*                                       |
 * | std::ranges::view                | *guaranteed*                                       |
 * | std::ranges::sized_range         | *guaranteed*                                       |
 * | std::ranges::common_range        |                                                    |
 * | std::ranges::output_range        | *guaranteed*                                       |
 * | seqan3::const_iterable_range     | *guaranteed*                                       |
 * |                                  |                                                    |
 * | std::ranges::range_reference_t   | std::remove_reference_t<value_t> &                 |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * \attention The given value to repeat is always **copied** into the range.
 *
 * ### Example
 *
 * \include test/snippet/utility/views/repeat_n.cpp
 *
 * \hideinitializer
 */
inline constexpr auto repeat_n = detail::repeat_n_fn{};

} // namespace seqan3::views
