// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::detail::simd_matrix_scoring_scheme.
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 */

#pragma once

#include <concepts>

#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/alignment/scoring/scoring_scheme_concept.hpp>
#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/utility/simd/algorithm.hpp>
#include <seqan3/utility/simd/concept.hpp>

namespace seqan3::detail
{

/*!\brief A vectorised scoring scheme to handle scoring matrices using gather strategy.
 * \ingroup alignment_scoring
 * \tparam simd_score_t The type of the simd vector; must model seqan3::simd::simd_concept.
 * \tparam alphabet_t The type of the alphabet over which to define the scoring scheme; must model seqan3::semialphabet.
 * \tparam alignment_t The type of the alignment to compute; must be either seqan3::align_cfg::method_global or
 *                     seqan3::align_cfg::method_local.
 *
 * \details
 *
 * When scoring two seqan3::detail::simd vectors, this performs element-wise lookups of the compared simd vectors
 * using a [gather operation](https://en.wikipedia.org/wiki/Gather-scatter_(vector_addressing)) on the scoring scheme.
 * The given scoring scheme is first transferred into linear memory such that a simple index gather can be used to
 * retrieve the actual score.
 * The index for the column vector of the alignment matrix must be precomputed using the
 * seqan3::detail::simd_matrix_scoring_scheme::make_score_profile member function.
 * This function computes the starting index of the respective matrix entry within the linearised
 * scoring scheme. To improve the performance this is only done once per column inside of the alignment algorithm.
 *
 * This simd scoring scheme matrix only needs one padding symbol, whose rank is initialised with the size of the
 * alphabet. Accordingly, the internal alphabet size increases by one.
 * Depending on the selected algorithm method the corresponding score values are either set to `1` for the global
 * alignment or `-1` for the local alignment. In the global alignment the score with a padding symbol will always
 * increase, and the projected maximum will be rescaled with the number of matches added.
 * In the local alignment the score with a padding symbol will always decrease, and the optimum can only be found inside
 * of the valid score matrix area.
 *
 * \note Note that the alphabet type information is lost during the conversion to the simd vectors and
 * only the ranks of the alphabet are used.
 */
template <simd_concept simd_score_t, semialphabet alphabet_t, typename alignment_t>
    requires (std::same_as<alignment_t, align_cfg::method_local> || std::same_as<alignment_t, align_cfg::method_global>)
class simd_matrix_scoring_scheme
{
private:
    //!\brief The underlying scalar type of the simd vector.
    using scalar_type = typename simd_traits<simd_score_t>::scalar_type;
    //!\brief The score profile type used for this scoring scheme, which is the same as the simd score type.
    using simd_score_profile_type = simd_score_t;
    //!\brief The type of the simd vector representing the alphabet ranks of one sequence batch.
    using simd_alphabet_ranks_type = simd_score_t;

    static_assert(seqan3::alphabet_size<alphabet_t> <= std::numeric_limits<scalar_type>::max(),
                  "The selected simd scalar type is not large enough to represent the given alphabet including an "
                  "additional padding symbol!");

    static_assert(seqan3::alphabet_size<alphabet_t> < std::numeric_limits<scalar_type>::max(),
                  "The selected simd scalar type is not large enough to represent the given alphabet including an "
                  "additional padding symbol!");

    //!\brief A flag that indicates wether the alignment mode is global.
    static constexpr bool is_global = std::same_as<alignment_t, align_cfg::method_global>;
    //!\brief The offset used to jump to the correct index in the linearised scoring scheme data.
    static constexpr size_t index_offset = seqan3::alphabet_size<alphabet_t> + 1; // scheme is extended by one.
    //!\brief The score used for the padding symbol (global -> increases score; local -> decreases score).
    static constexpr scalar_type score_for_padding_symbol = (is_global) ? 1 : -1;

    //!\brief The scoring scheme stored as a linear array.
    std::vector<scalar_type> scoring_scheme_data{};

public:
    //!\brief The padding symbol used to fill up smaller sequences in a simd batch.
    static constexpr scalar_type padding_symbol = static_cast<scalar_type>(seqan3::alphabet_size<alphabet_t>);

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr simd_matrix_scoring_scheme() = default;                                               //!< Defaulted.
    constexpr simd_matrix_scoring_scheme(simd_matrix_scoring_scheme const &) = default;             //!< Defaulted.
    constexpr simd_matrix_scoring_scheme(simd_matrix_scoring_scheme &&) = default;                  //!< Defaulted.
    constexpr simd_matrix_scoring_scheme & operator=(simd_matrix_scoring_scheme const &) = default; //!< Defaulted.
    constexpr simd_matrix_scoring_scheme & operator=(simd_matrix_scoring_scheme &&) = default;      //!< Defaulted.
    ~simd_matrix_scoring_scheme() = default;                                                        //!< Defaulted.

    //!\copydoc seqan3::detail::simd_matrix_scoring_scheme::initialise_from_scalar_scoring_scheme
    template <typename scoring_scheme_t>
    constexpr explicit simd_matrix_scoring_scheme(scoring_scheme_t const & scoring_scheme)
    {
        initialise_from_scalar_scoring_scheme(scoring_scheme);
    }

    //!\copydoc seqan3::detail::simd_matrix_scoring_scheme::initialise_from_scalar_scoring_scheme
    template <typename scoring_scheme_t>
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
     * \param[in] score_profile The precomputed score profile.
     * \param[in] ranks A simd vector over alphabet ranks.
     * \returns A score simd vector computed based on the given score profile and alphabet ranks.
     *
     * ### Exception
     *
     * No-throw guarantee.
     *
     * ### Complexity
     *
     * Linear in the length of one input vector (`score_profile` and `ranks` are equally sized).
     *
     * ### Thread safety
     *
     * Thread-safe.
     *
     * \attention You need to call seqan3::detail::simd_matrix_scoring_scheme::make_score_profile before invoking
     *            the score interface to convert the column alphabet ranks into a column profile. Not doing this cannot
     *            be detected by the program and can lead to wrong results.
     */
    constexpr simd_score_t score(simd_score_profile_type const & score_profile,
                                 simd_alphabet_ranks_type const & ranks) const noexcept
    {
        simd_score_t const matrix_index = score_profile + ranks; // Compute the matrix indices for the lookup.
        simd_score_t result{};

        for (size_t idx = 0; idx < simd_traits<simd_score_t>::length; ++idx)
            result[idx] = scoring_scheme_data.data()[matrix_index[idx]];

        return result;
    }
    //!\}

    //!\brief Returns the score used when aligning a padding symbol.
    constexpr scalar_type padding_match_score() const noexcept
    {
        return score_for_padding_symbol;
    }

    /*!\brief Converts the simd alphabet ranks into a score profile used for scoring it later with the alphabet ranks
     *        of another sequence batch.
     *
     * \details
     *
     * In this implementation of the simd matrix gather operations are used to select the score values from the
     * underlying scoring scheme. Since the scoring scheme matrix is represented as a linear vector, the corresponding
     * indices are computed by the alphabet rank of one sequence batch times the alphabet size plus the alphabet rank
     * of another sequence batch.
     */
    constexpr simd_score_profile_type make_score_profile(simd_alphabet_ranks_type const & ranks) const noexcept
    {
        return ranks * simd::fill<simd_score_t>(index_offset);
    }

private:
    /*!\brief Store the given scoring scheme matrix into a private member variable.
     * \tparam scoring_scheme_t The type of the scoring scheme; must model seqan3::scoring_scheme_for the given
     *                          alphabet type.
     * \param[in] scoring_scheme The scoring scheme to initialise the vectorised match and mismatch score with.
     *
     * \throws std::invalid_argument if the score of the given scoring scheme exceeds the score range covered by the
     *         selected simd vector type.
     */
    template <typename scoring_scheme_t>
        requires scoring_scheme_for<scoring_scheme_t, alphabet_t>
    constexpr void initialise_from_scalar_scoring_scheme(scoring_scheme_t const & scoring_scheme)
    {
        using score_t = decltype(std::declval<scoring_scheme_t const &>().score(alphabet_t{}, alphabet_t{}));

        // Helper function to check if the matrix score is in the value range representable by the selected simd vector.
        [[maybe_unused]] auto check_score_range = [&]([[maybe_unused]] score_t score)
        {
            // Note only if the size of the scalar type of the simd vector is smaller than the size of the score type
            // of the original scoring scheme, the score might exceed the valid value range of the scalar type. In this
            // case an exception will be thrown.
            if constexpr (sizeof(scalar_type) < sizeof(score_t))
            {
                constexpr score_t max_score_value = static_cast<score_t>(std::numeric_limits<scalar_type>::max());
                constexpr score_t min_score_value = static_cast<score_t>(std::numeric_limits<scalar_type>::lowest());

                if (score > max_score_value || score < min_score_value)
                    throw std::invalid_argument{"The selected scoring scheme score overflows "
                                                "for the selected scalar type of the simd type."};
            }
        };

        // For the global alignment we extend the alphabet by one symbol to handle sequences with different size.
        scoring_scheme_data.resize(index_offset * index_offset, score_for_padding_symbol);

        // Convert the scoring matrix into a linear vector to allow gather operations later on.
        using alphabet_size_t = std::remove_const_t<decltype(seqan3::alphabet_size<alphabet_t>)>;
        auto data_it = scoring_scheme_data.begin();
        for (alphabet_size_t lhs_rank = 0; lhs_rank < seqan3::alphabet_size<alphabet_t>; ++lhs_rank)
        {
            for (alphabet_size_t rhs_rank = 0; rhs_rank < seqan3::alphabet_size<alphabet_t>; ++rhs_rank, ++data_it)
            {
                score_t tmp_score = scoring_scheme.score(seqan3::assign_rank_to(lhs_rank, alphabet_t{}),
                                                         seqan3::assign_rank_to(rhs_rank, alphabet_t{}));

                check_score_range(tmp_score);
                *data_it = tmp_score;
            }
            ++data_it; // skip one for the padded symbol.
        }
    }
};

} // namespace seqan3::detail
