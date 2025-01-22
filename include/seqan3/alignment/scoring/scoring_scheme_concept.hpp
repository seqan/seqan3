// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::scoring_scheme_for.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>

namespace seqan3
{

/*!\interface seqan3::scoring_scheme_for <>
 * \brief A concept that requires that type be able to score two letters.
 * \tparam t            The type the concept check is performed on (the putative scoring scheme).
 * \tparam alphabet_t   The type of the first letter that you wish to score; must model seqan3::alphabet.
 * \tparam alphabet2_t  The type of the second letter that you wish to score; must model seqan3::alphabet;
 *                      defaults to `alphabet_t`.
 * \ingroup alignment_scoring
 *
 * \details
 *
 * This concept makes no assumptions about configurability or assignability of the scoring scheme, only the
 * ability to score the two letters is required.
 *
 */
/*!\name Requirements for seqan3::scoring_scheme_for
 * \brief You can expect these members on all types that implement seqan3::scoring_scheme_for.
 * \memberof seqan3::scoring_scheme_for
 * \{
 */

/*!\typedef     typedef IMPLEMENTATION_DEFINED score_type;
 * \brief       The type returned by seqan3::scoring_scheme_for::score(), usually a seqan3::arithmetic.
 *
 * \details
 * \attention This is a concept requirement, not an actual typedef (however types satisfying this concept
 * will provide an implementation).
 */

/*!\fn          score_type score(alph1_t const alph1, alph2_t const alph2);
 * \brief       Compute the score of two letters
 * \param alph1 First letter.
 * \param alph2 Second letter.
 *
 * \details
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */
//!\}

//!\cond
template <typename t, typename alphabet_t, typename alphabet2_t = alphabet_t>
concept scoring_scheme_for = requires (t scheme, alphabet_t const alph1, alphabet2_t const alph2) {
    requires alphabet<alphabet_t>;
    requires alphabet<alphabet2_t>;

    { scheme.score(alph1, alph2) };
    requires std::common_reference_with<decltype(scheme.score(alph1, alph2)),
                                        typename std::remove_reference_t<t>::score_type>;

    { scheme.score(alphabet_t{}, alphabet2_t{}) };
    requires std::common_reference_with<decltype(scheme.score(alphabet_t{}, alphabet2_t{})),
                                        typename std::remove_reference_t<t>::score_type>;
};
//!\endcond

} // namespace seqan3
