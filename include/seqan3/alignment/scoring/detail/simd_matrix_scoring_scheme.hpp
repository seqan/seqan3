// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::simd_matrix_scoring_scheme.
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
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

/*!\brief A vectorised scoring scheme to handle scoring matrices.
 * \ingroup scoring
 * \tparam simd_score_t The type of the simd vector; must model seqan3::simd::simd_concept.
 * \tparam alphabet_t The type of the alphabet over which to define the scoring scheme; must model seqan3::semialphabet.
 * \tparam alignment_t The type of the alignment to compute; must be either seqan3::detail::global_alignment_type or
 *                     seqan3::detail::local_alignment_type.
 * \tparam scoring_scheme_t The type of the scoring scheme which will be used. Must model seqan3::scoring_scheme with
 *                          the given `alphabet_t`.
 *
 * \details
 *
 * When scoring two seqan3::detail::simd vectors, this performs element-wise lookups of the compared simd vectors
 * using the underlying scoring matrix and returns the result in another seqan3::detail::simd vector.
 *
 * \note Note that the alphabet type information is lost during the conversion to the simd vectors and
 * only the ranks of the alphabet are used.
 */
template <simd_concept simd_score_t, semialphabet alphabet_t, typename alignment_t, typename scoring_scheme_t>
//!\cond
    requires scoring_scheme<scoring_scheme_t, alphabet_t> &&
             (std::same_as<alignment_t, detail::local_alignment_type> ||
              std::same_as<alignment_t, detail::global_alignment_type>)
//!\endcond
class simd_matrix_scoring_scheme
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr simd_matrix_scoring_scheme() = default; //!< Defaulted.
    constexpr simd_matrix_scoring_scheme(simd_matrix_scoring_scheme const &) = default; //!< Defaulted.
    constexpr simd_matrix_scoring_scheme(simd_matrix_scoring_scheme &&) = default; //!< Defaulted.
    constexpr simd_matrix_scoring_scheme & operator=(simd_matrix_scoring_scheme const &) = default; //!< Defaulted.
    constexpr simd_matrix_scoring_scheme & operator=(simd_matrix_scoring_scheme &&) = default; //!< Defaulted.
    ~simd_matrix_scoring_scheme() = default; //!< Defaulted.

    //!\copydoc seqan3::detail::simd_matrix_scoring_scheme::initialise_from_scalar_scoring_scheme
    constexpr explicit simd_matrix_scoring_scheme(scoring_scheme_t const & scoring_scheme)
    {
        initialise_from_scalar_scoring_scheme(scoring_scheme);
    }

    //!\copydoc seqan3::detail::simd_matrix_scoring_scheme::initialise_from_scalar_scoring_scheme
    constexpr simd_matrix_scoring_scheme & operator=(scoring_scheme_t const & scoring_scheme)
    {
        initialise_from_scalar_scoring_scheme(scoring_scheme);
        return *this;
    }
    //!\}

    /*!\name Score computation
     *\{
     */
    /*!\brief Given two simd vectors, compute an element-wise score and return another simd vector.
     * \param[in] lhs The left operand to compare.
     * \param[in] rhs The right operand to compare.
     * \returns A simd vector of the same type as the input, with scores computed between both inputs.
     *
     * ### Exception
     *
     * No-throw guarantee.
     *
     * ### Complexity
     *
     * Linear in the length of one input vector (`lhs` and `rhs` are equally sized).
     *
     * ### Thread safety
     *
     * Thread-safe.
     */
    constexpr simd_score_t score(simd_score_t const & lhs, simd_score_t const & rhs) const noexcept
    {
        simd_score_t result{};

        for (size_t i = 0; i < simd_traits<simd_score_t>::length; ++i)
        {
            if (is_padded(lhs[i]) || is_padded(rhs[i]))
            {
                if constexpr (std::same_as<alignment_t, detail::global_alignment_type>)
                    result[i] = 1;
                else
                    result[i] = -1;
            }
            else
            {
                result[i] = internal_scoring_scheme.score(assign_rank_to(lhs[i], alphabet_t{}),
                                                          assign_rank_to(rhs[i], alphabet_t{}));
            }
        }

        return result;
    }
    //!\}

    //!\brief Returns the match score used for padded symbols.
    constexpr typename scoring_scheme_t::score_type padding_match_score() noexcept
    {
        return 1;
    }

private:
    //!\brief Internally stores the given scalar scoring scheme matrix.
    scoring_scheme_t internal_scoring_scheme;

    /*!\brief Store the given scoring scheme matrix into a private member variable.
     * \param[in] scoring_scheme The scoring scheme to initialise the vectorised match and mismatch score with.
     */
    constexpr void initialise_from_scalar_scoring_scheme(scoring_scheme_t const & scoring_scheme)
    {
        using score_t = decltype(std::declval<scoring_scheme_t const &>().score(alphabet_t{}, alphabet_t{}));
        using simd_scalar_t = typename simd_traits<simd_score_t>::scalar_type;

        // Check if the scoring scheme match and mismatch scores do not overflow with the respective scalar type.
        if constexpr (sizeof(simd_scalar_t) < sizeof(score_t))
        {
            if (min_or_max_exceeded<score_t, simd_scalar_t>(scoring_scheme))
                throw std::invalid_argument{"The selected scoring scheme score overflows "
                                            "for the selected scalar type of the simd type."};
        }

        internal_scoring_scheme = scoring_scheme;
    }

    /*!\brief Check if any score in the scoring scheme matrix exceeds the min or max value allowed by the simd vector.
     * \param[in] scoring_scheme The scoring scheme to check for min or max violations.
     * \tparam score_t The type of the scores in the scoring scheme matrix.
     * \tparam simd_scalar_t The type of the simd vector.
     * \returns Returns `true` if it encounters a value larger than the maximum or lower than the minimum allowed by
     *          the simd vector, and `false` if otherwise.
     */
    template <typename score_t, typename simd_scalar_t>
    constexpr bool min_or_max_exceeded(scoring_scheme_t const & scoring_scheme)
    {
        score_t max_val = static_cast<score_t>(std::numeric_limits<simd_scalar_t>::max());
        score_t min_val = static_cast<score_t>(std::numeric_limits<simd_scalar_t>::lowest());

        for (uint8_t i = 0; i < alphabet_size<alphabet_t>; ++i)
        {
            for (uint8_t j = 0; j < alphabet_size<alphabet_t>; ++j)
            {
                score_t tmp_score = scoring_scheme.score(assign_rank_to(i, alphabet_t{}),
                                                         assign_rank_to(j, alphabet_t{}));

                if (tmp_score > max_val || tmp_score < min_val)
                    return true;
            }
        }

        return false;
    }

    /*!\brief Checks whether the given value is a padded value or not.
     * \param[in] value The value to check if it is padded.
     * \returns Returns `true` if the value is padded and `false` otherwise.
     *
     * \details This check is done by casting the scalar type of the value into an unsigned value, and then checking
     *          to see if this is larger than the alphabet size, which would imply that it is a padded value.
     */
    constexpr bool is_padded(typename simd_traits<simd_score_t>::scalar_type const & value) const
    {
        using unsigned_scalar_t = typename std::make_unsigned<typename simd_traits<simd_score_t>::scalar_type>::type;

        return static_cast<unsigned_scalar_t>(value) >= alphabet_size<alphabet_t>;
    }
};

} // namespace seqan3::detail
