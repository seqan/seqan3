// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::view::trim.
 */

#pragma once

#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/core/type_traits/all.hpp>
#include <seqan3/range/view/deep.hpp>
#include <seqan3/range/view/take_until.hpp>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

/*!\brief The underlying type of seqan3::view::trim.
 * \ingroup view
 *
 * Under the hood this delegates to view::take_until.
 */
struct trim_fn
{
    //!\brief Store the argument and return a range adaptor closure object.
    template <typename threshold_t>
    constexpr auto operator()(threshold_t const threshold) const
    {
        static_assert(QualityAlphabet<threshold_t> || std::Integral<threshold_t>,
                      "The threshold must either be a quality alphabet or an integral type "
                      "in which case it is compared with the underlying phred type.");

        return adaptor_from_functor{*this, threshold};
    }

    /*!\brief Trim based on minimum phred score.
     * \tparam irng_t The type of the range being processed. See seqan3::view::trim for requirements.
     * \param irange The range being processed.
     * \param threshold The minimum quality as a phred score [integral type].
     */
    template <std::ranges::InputRange irng_t, typename threshold_t>
    constexpr auto operator()(irng_t && irange, threshold_t const threshold) const
    {
        static_assert(QualityAlphabet<std::remove_reference_t<reference_t<irng_t>>>,
                      "view:trim can only operate on ranges over seqan3::QualityAlphabet.");
        static_assert(std::Same<remove_cvref_t<threshold_t>, remove_cvref_t<reference_t<irng_t>>> ||
                      std::Integral<remove_cvref_t<threshold_t>>,
                      "The threshold must either be a letter of the underlying alphabet or an integral type "
                      "in which case it is compared with the underlying phred type.");

        return view::take_until(std::forward<irng_t>(irange), [threshold] (auto const value)
        {
            if constexpr (std::Same<remove_cvref_t<threshold_t>, remove_cvref_t<reference_t<irng_t>>>)
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

namespace seqan3::view
{

/*!\name Alphabet related views
 * \{
 */

/*!\brief               A view that does quality-threshold trimming on a range of seqan3::QualityAlphabet.
 * \tparam urng_t       The type of the range being processed. See below for requirements.
 * \tparam threshold_t  Either seqan3::value_type_t<urng_t> or
 *                      seqan3::alphabet_phred_t<seqan3::value_type_t<urng_t>>.
 * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
 * \param[in] threshold The minimum quality.
 * \returns             A trimmed range. See below for the properties of the returned range.
 * \ingroup view
 *
 * \details
 *
 * This view can be used to do easy quality based trimming of sequences.
 *
 * **Header**
 * ```cpp
 *      #include <seqan3/range/view/trim.hpp>
 * ```
 *
 * ### View properties
 *
 * This view is a **deep view:** Given a range-of-range as input (as opposed to just a range), it will apply
 * the transformation on the innermost range (instead of the outermost range).
 *
 * | range concepts and reference_t  | `urng_t` (underlying range type)      | `rrng_t` (returned range type)  |
 * |---------------------------------|:-------------------------------------:|:-------------------------------:|
 * | std::ranges::InputRange         | *required*                            | *preserved*                     |
 * | std::ranges::ForwardRange       |                                       | *preserved*                     |
 * | std::ranges::BidirectionalRange |                                       | *preserved*                     |
 * | std::ranges::RandomAccessRange  |                                       | *preserved*                     |
 * |                                 |                                       |                                 |
 * | std::ranges::View               |                                       | *guaranteed*                    |
 * | std::ranges::SizedRange         |                                       | *lost*                          |
 * | std::ranges::CommonRange        |                                       | *lost*                          |
 * | std::ranges::OutputRange        |                                       | *preserved*                     |
 * | seqan3::ConstIterableRange      |                                       | *preserved*                     |
 * |                                 |                                       |                                 |
 * | seqan3::reference_t             | seqan3::QualityAlphabet               | seqan3::reference_t<urng_t>     |
 *
 * See the \link view view submodule documentation \endlink for detailed descriptions of the view properties.
 *
 * ###Example
 *
 * Operating on a range of seqan3::phred42:
 * \snippet test/snippet/range/view/trim.cpp phred42
 *
 * Or operating on a range of seqan3::dna5q:
 * \snippet test/snippet/range/view/trim.cpp dna5q
 * \hideinitializer
 */

inline constexpr auto trim = deep{seqan3::detail::trim_fn{}};

//!\}

} // namespace seqan3::view
