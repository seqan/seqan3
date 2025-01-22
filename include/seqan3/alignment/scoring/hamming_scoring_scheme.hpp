// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::hamming_scoring_scheme.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>

#include <seqan3/core/platform.hpp>

namespace seqan3
{

/*!\brief A scoring scheme that assigns a score of `0` to matching letters and `-1` to mismatching letters.
 * \ingroup alignment_scoring
 *
 * \details
 *
 * This stateless scoring scheme is equivalent to the Hamming distance and assigns a score of `0` to matching letters
 * and `-1` to mismatching letters.
 * This scheme is independent of the alphabet type and can be used whenever the two compared alphabets model the
 * std::equality_comparable_with concept.
 * Use this scheme if you want to use the use the more efficient bitparallel alignment algorithm often used in the
 * context of computing the edit distance.
 * See the documentation for \ref seqan3::align_cfg::edit_scheme for more details.
 */
class hamming_scoring_scheme
{
public:
    //!\privatesection
    //!\brief The alphabet type of the scoring scheme. This type is only used for the alignment configuration machinery.
    using alphabet_type = char;
    //!\publicsection

    //!\brief The underlying score type.
    using score_type = int32_t;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    hamming_scoring_scheme() noexcept = default;                                      //!< Defaulted.
    hamming_scoring_scheme(hamming_scoring_scheme const &) = default;                 //!< Defaulted.
    hamming_scoring_scheme(hamming_scoring_scheme &&) noexcept = default;             //!< Defaulted.
    hamming_scoring_scheme & operator=(hamming_scoring_scheme const &) = default;     //!< Defaulted.
    hamming_scoring_scheme & operator=(hamming_scoring_scheme &&) noexcept = default; //!< Defaulted.
    ~hamming_scoring_scheme() = default;                                              //!< Defaulted.
    //!\}

    /*!\name Accessors
     * \{
     */

    /*!\brief Returns the score of two letters.
     * \tparam alph1_t The type of the first letter.
     * \tparam alph2_t The type of the second letter.
     * \param[in] alph1 The first letter to compare.
     * \param[in] alph2 The second letter to compare.
     *
     * \details
     *
     * This function returns `0` if the two letters are equal and `-1` otherwise. Note that both alphabet types of the
     * compared letters must model the std::equality_comparable_with concept or no overload of this function
     * is available.
     *
     * \returns The score of the two letters. `0` if the letters are equal, `-1` otherwise.
     */
    template <typename alph1_t, typename alph2_t>
        requires std::equality_comparable_with<alph1_t, alph2_t>
    constexpr score_type score(alph1_t const alph1, alph2_t const alph2) const noexcept
    {
        return alph1 == alph2 ? 0 : -1;
    }
    //!\}

    //!\name Comparison operators
    //!\{

    //!\brief Always true.
    constexpr bool operator==(hamming_scoring_scheme const &) const noexcept
    {
        return true;
    }

    //!\brief Always false.
    constexpr bool operator!=(hamming_scoring_scheme const &) const noexcept
    {
        return false;
    }
    //!\}
};

} // namespace seqan3
