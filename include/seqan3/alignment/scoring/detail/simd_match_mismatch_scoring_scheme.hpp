// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::simd_match_mismatch_scoring_scheme.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/align_config_mode.hpp>
#include <seqan3/alignment/scoring/scoring_scheme_concept.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/core/simd/concept.hpp>
#include <seqan3/core/simd/simd_algorithm.hpp>
#include <seqan3/std/concepts>

namespace seqan3::detail
{

/*!\brief A vectorised scoring scheme handling matches and mismatches only.
 * \ingroup scoring
 * \tparam simd_score_t The type of the simd vector; must model seqan3::detail::simd_concept.
 * \tparam alphabet_t The type of the alphabet over which to define the scoring scheme; must model seqan3::semialphabet
 *                    and must have an alphabet size of at least 2.
 * \tparam alignment_t The type of the alignment to compute; must be either seqan3::detail::global_alignment_type or
 *                     seqan3::detail::local_alignment_type.
 *
 * \details
 *
 * Wraps a regular scoring scheme by extracting the scores for a match and a mismatch and converts them into
 * seqan3::detail::simd vectors. Only symmetric scoring schemes are preserved, i.e. in the vectorised scoring
 * scheme elements with the same rank are assigned the match score and elements with a different rank are assigned the
 * mismatch score.
 * Note during the conversion to the simd vectors the alphabet type information is lost and
 * only the ranks of the alphabet are used.
 *
 * ### Handling special padding symbols
 *
 * During the vectorised alignment multiple sequences are packed into one simd vector.
 * To handle sequences with different sizes in the vectorised alignment algorithm the smaller sequences are filled up
 * with special padding symbols. These padding symbols are chosen in a way that allows the computation of the alignments
 * without the need of masking the results for invalid positions within the matrix because a specific position might
 * have exceeded the original sequence size.
 * To do so, the global alignment uses the same padding symbol for the first sequence pack and the second sequence
 * pack. This padding symbol is distinct to any symbol in the underlying alphabet of the sequences.
 * The score function is adapted in a way that a comparison with a padding symbol always yields a match.
 * Thus, after the end of a sequence within the pack is reached the score can only grow.
 * The respective score can then be inferred from the projected position of the last row or column of the
 * vectorised matrix depending on the the corresponding alignment configuration.
 *
 * In case of the local alignment the second sequence pack are padded with a symbol that is distinct to any symbol
 * of the corresponding alphabet and to the padding symbol of the first sequence pack.
 * Comparing any symbol with the padding symbols will yield a mismatch, such that the score can only get smaller after
 * the end of a sequence has reached. This way the specific optimum of one sequence pair in the pack is not affected
 * during the computation of the vectorised alignment.
 */
template <simd_concept simd_score_t, semialphabet alphabet_t, typename alignment_t>
//!\cond
    requires seqan3::alphabet_size<alphabet_t> > 1 &&
             (std::same_as<alignment_t, detail::local_alignment_type> ||
              std::same_as<alignment_t, detail::global_alignment_type>)
//!\endcond
class simd_match_mismatch_scoring_scheme
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Defaulted.
    constexpr simd_match_mismatch_scoring_scheme() = default;
    //!\brief Defaulted.
    constexpr simd_match_mismatch_scoring_scheme(simd_match_mismatch_scoring_scheme const &) = default;
    //!\brief Defaulted.
    constexpr simd_match_mismatch_scoring_scheme(simd_match_mismatch_scoring_scheme &&) = default;
    //!\brief Defaulted.
    constexpr simd_match_mismatch_scoring_scheme & operator=(simd_match_mismatch_scoring_scheme const &) = default;
    //!\brief Defaulted.
    constexpr simd_match_mismatch_scoring_scheme & operator=(simd_match_mismatch_scoring_scheme &&) = default;
    //!\brief Defaulted.
    ~simd_match_mismatch_scoring_scheme() = default;

    //!\copydoc seqan3::detail::simd_match_mismatch_scoring_scheme::initialise_from_scalar_scoring_scheme
    template <typename scoring_scheme_t>
    //!\cond
        requires scoring_scheme<scoring_scheme_t, alphabet_t>
    //!\endcond
    constexpr explicit simd_match_mismatch_scoring_scheme(scoring_scheme_t const & scoring_scheme)
    {
        initialise_from_scalar_scoring_scheme(scoring_scheme);
    }

    //!\copydoc seqan3::detail::simd_match_mismatch_scoring_scheme::initialise_from_scalar_scoring_scheme
    template <typename scoring_scheme_t>
    //!\cond
        requires scoring_scheme<scoring_scheme_t, alphabet_t>
    //!\endcond
    constexpr simd_match_mismatch_scoring_scheme & operator=(scoring_scheme_t const & scoring_scheme)
    {
        initialise_from_scalar_scoring_scheme(scoring_scheme);
        return *this;
    }
    //!\}

    /*!\name Score computation
     * \{
     */
    /*!\brief Computes the score for two simd vectors.
     * \param[in] lhs The left operand to compare.
     * \param[in] rhs The right operand to compare.
     *
     * \returns The simd score with match and mismatch scores after comparing both input operands.
     *
     * \details
     *
     * This function compares packed elements in both simd vectors and returns a new simd vector filled with match and
     * mismatch scores depending on the result of the comparison. For global alignments the comparison yields a match
     * if any of the elements is a padding symbol. The padding symbol must have the signed bit set.
     * For local alignments two different padding symbols are assumed, which will always yield a mismatch when compared
     * with any other symbol.
     *
     * ### Exception
     *
     * No-throw guarantee.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Thread safety
     *
     * Thread-safe.
     */
    constexpr simd_score_t score(simd_score_t const & lhs, simd_score_t const & rhs) const noexcept
    {
        typename simd_traits<simd_score_t>::mask_type mask;
        // For global and local alignment there are slightly different formulas because
        // in global alignment padded characters always match
        if constexpr (std::same_as<alignment_t, detail::global_alignment_type>)
            mask = (lhs ^ rhs) <= simd::fill<simd_score_t>(0);
        else // and in local alignment type padded characters always mismatch.
            mask = (lhs ^ rhs) == simd::fill<simd_score_t>(0);

        return mask ? match_score : mismatch_score;
    }
    //!\}

    //!\brief Returns the match score used for padded symbols.
    constexpr auto padding_match_score() noexcept
    {
        return match_score[0];
    }

private:
    /*!\brief Initialises the simd vector match score and mismatch score from the given scoring scheme.
     * \tparam scoring_scheme_t The type of the underlying scoring scheme; must model seqan3::scoring_scheme for
     *                          `alphabet_t`.
     *
     * \param[in] scoring_scheme The scoring scheme to initialise the vectorised match and mismatch score from.
     *
     * \throws std::invalid_argument if the value of the match or mismatch score exceed the value range of the
     *         scalar type of the used simd vector type `simd_score_t`.
     *
     * \details
     *
     * Obtains the score for a match and a mismatch respectively and fills the corresponding simd vectors with
     * these scores. In addition, some safety checks are performed in order to avoid that a score is used which
     * cannot be represented by the scalar type of the used simd vector.
     */
    template <typename scoring_scheme_t>
    constexpr void initialise_from_scalar_scoring_scheme(scoring_scheme_t const & scoring_scheme)
    {
        using score_t = decltype(std::declval<scoring_scheme_t const &>().score(alphabet_t{}, alphabet_t{}));
        using simd_scalar_t = typename simd_traits<simd_score_t>::scalar_type;

        score_t scalar_match_score = scoring_scheme.score(seqan3::assign_rank_to(0, alphabet_t{}),
                                                          seqan3::assign_rank_to(0, alphabet_t{}));
        score_t scalar_mismatch_score = scoring_scheme.score(seqan3::assign_rank_to(0, alphabet_t{}),
                                                             seqan3::assign_rank_to(1, alphabet_t{}));

        // Check if the scoring scheme match and mismatch scores do not overflow with the respective scalar type.
        if constexpr (sizeof(simd_scalar_t) < sizeof(score_t))
        {
            if (scalar_match_score > static_cast<score_t>(std::numeric_limits<simd_scalar_t>::max()) ||
                scalar_mismatch_score < static_cast<score_t>(std::numeric_limits<simd_scalar_t>::lowest()))
            {
                throw std::invalid_argument{"The selected scoring scheme score overflows "
                                            "for the selected scalar type of the simd type."};
            }
        }

        match_score = simd::fill<simd_score_t>(static_cast<simd_scalar_t>(scalar_match_score));
        mismatch_score = simd::fill<simd_score_t>(static_cast<simd_scalar_t>(scalar_mismatch_score));
    }

    simd_score_t match_score; //!< The simd vector for a match score.
    simd_score_t mismatch_score; //!< The simd vector for a mismatch score.
};

} // namespace seqan3::detail
