// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::views::trim_quality.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <ranges>

#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/io/views/detail/take_until_view.hpp>
#include <seqan3/utility/views/deep.hpp>

namespace seqan3::detail
{

/*!\brief The underlying type of seqan3::views::trim_quality.
 * \ingroup alphabet_views
 *
 * Under the hood, this delegates to std::views::take_while.
 */
struct trim_fn
{
    //!\brief Store the argument and return a range adaptor closure object.
    template <typename threshold_t>
    constexpr auto operator()(threshold_t const threshold) const
    {
        static_assert(quality_alphabet<threshold_t> || std::integral<threshold_t>,
                      "The threshold must either be a quality alphabet or an integral type "
                      "in which case it is compared with the underlying Phred score type.");

        return adaptor_from_functor{*this, threshold};
    }

    /*!\brief Trim based on minimum Phred score.
     * \tparam irng_t The type of the range being processed. See seqan3::views::trim_quality for requirements.
     * \param irange The range being processed.
     * \param threshold The minimum quality as a Phred score [integral type].
     */
    template <std::ranges::input_range irng_t, typename threshold_t>
    constexpr auto operator()(irng_t && irange, threshold_t const threshold) const
    {
        static_assert(quality_alphabet<std::remove_reference_t<std::ranges::range_reference_t<irng_t>>>,
                      "views::trim_quality can only operate on ranges over seqan3::quality_alphabet.");
        static_assert(
            std::same_as<std::remove_cvref_t<threshold_t>, std::remove_cvref_t<std::ranges::range_reference_t<irng_t>>>
                || std::integral<std::remove_cvref_t<threshold_t>>,
            "The threshold must either be a letter of the underlying alphabet or an integral type "
            "in which case it is compared with the underlying Phred score type.");

        return detail::take_until(
            std::forward<irng_t>(irange),
            [threshold](auto const value)
            {
                if constexpr (std::same_as<std::remove_cvref_t<threshold_t>,
                                           std::remove_cvref_t<std::ranges::range_reference_t<irng_t>>>)
                {
                    return to_phred(value) < to_phred(threshold);
                }
                else
                {
                    using c_t = std::common_type_t<decltype(to_phred(value)), threshold_t>;
                    return static_cast<c_t>(to_phred(value)) < static_cast<c_t>(threshold);
                }
            });
    }
};

} // namespace seqan3::detail

namespace seqan3::views
{
/*!\brief               A view that does quality-threshold trimming on a range of seqan3::quality_alphabet.
 * \tparam urng_t       The type of the range being processed. See below for requirements.
 * \tparam threshold_t  Either std::ranges::range_value_t<urng_t> or
 *                      seqan3::alphabet_phred_t<std::ranges::range_value_t<urng_t>>.
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \param[in] threshold The minimum quality.
 * \returns             A trimmed range. See below for the properties of the returned range.
 * \ingroup alphabet_views
 *
 * \details
 *
 * \header_file{seqan3/alphabet/views/trim_quality.hpp}
 *
 * This view can be used to trim sequences based on quality. Only bases at the end of the sequence not meeting the
 * specified threshold are discarded.
 *
 * ### View properties
 *
 * This view is a **deep view** Given a range-of-range as input (as opposed to just a range), it will apply
 * the transformation on the innermost range (instead of the outermost range).
 *
 * | Concepts and traits              | `urng_t` (underlying range type)      | `rrng_t` (returned range type)         |
 * |----------------------------------|:-------------------------------------:|:--------------------------------------:|
 * | std::ranges::input_range         | *required*                            | *preserved*                            |
 * | std::ranges::forward_range       |                                       | *preserved*                            |
 * | std::ranges::bidirectional_range |                                       | *preserved*                            |
 * | std::ranges::random_access_range |                                       | *preserved*                            |
 * |                                  |                                       |                                        |
 * | std::ranges::view                |                                       | *guaranteed*                           |
 * | std::ranges::sized_range         |                                       | *lost*                                 |
 * | std::ranges::common_range        |                                       | *lost*                                 |
 * | std::ranges::output_range        |                                       | *preserved*                            |
 * | seqan3::const_iterable_range     |                                       | *preserved*                            |
 * |                                  |                                       |                                        |
 * | std::ranges::range_reference_t   | seqan3::quality_alphabet              | std::ranges::range_reference_t<urng_t> |
 *
 * See the \link views views submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ###Example
 *
 * Operating on a range of seqan3::phred42:
 * \include test/snippet/alphabet/views/trim_quality_phred42.cpp
 *
 * Or operating on a range of seqan3::dna5q:
 * \include test/snippet/alphabet/views/trim_quality_dna5q.cpp
 * \hideinitializer
 *
 * \experimentalapi{Experimental since version 3.1.}
 */
inline constexpr auto trim_quality = deep{seqan3::detail::trim_fn{}};

} // namespace seqan3::views
