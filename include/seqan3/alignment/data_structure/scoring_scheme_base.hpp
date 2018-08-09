// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

/*!\file
 * \brief Provides seqan3::scoring_scheme_base.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <range/v3/algorithm/copy.hpp>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/concept/cereal.hpp>
#include <seqan3/core/detail/strong_type.hpp>

#if SEQAN3_WITH_CEREAL
#include <cereal/types/array.hpp>
#endif // SEQAN3_WITH_CEREAL

namespace seqan3
{

// ------------------------------------------------------------------
// seqan3::match_score
// ------------------------------------------------------------------

/*!\brief A strong type of underlying type `score_type` that represents the score of two matching characters.
 * \tparam score_type The underlying type.
 * \ingroup data_structure
 * \see scoring_scheme_base::set_simple_scheme
 */
template <arithmetic_concept score_type>
struct match_score : detail::strong_type<score_type, match_score<score_type>, detail::strong_type_skill::convert>
{
     using detail::strong_type<score_type, match_score<score_type>, detail::strong_type_skill::convert>::strong_type;
};

/*!\name Template argument type deduction guides
 * \relates seqan3::match_score
 * \{
 */
template <arithmetic_concept score_type>
match_score(score_type &&) -> match_score<score_type>;
//!\}

// ------------------------------------------------------------------
// seqan3::mismatch_score
// ------------------------------------------------------------------

/*!\brief A strong type of underlying type `score_type` that represents the score two different characters.
 * \tparam score_type The underlying type.
 * \ingroup data_structure
 * \see scoring_scheme_base::set_simple_scheme
 */
template <arithmetic_concept score_type>
struct mismatch_score : detail::strong_type<score_type, mismatch_score<score_type>, detail::strong_type_skill::convert>
{
     using detail::strong_type<score_type, mismatch_score<score_type>, detail::strong_type_skill::convert>::strong_type;
};

/*!\name Template argument type deduction guides
 * \relates seqan3::mismatch_score
 * \{
 */
template <arithmetic_concept score_type>
mismatch_score(score_type &&) -> mismatch_score<score_type>;
//!\}

// ------------------------------------------------------------------
// seqan3::scoring_scheme_base
// ------------------------------------------------------------------

/*!\brief A CRTP base class for scoring schemes.
 * \tparam derived_t  The derived type.
 * \tparam alphabet_t Type of the largest target alphabet.
 * \tparam score_type Type of the score values in the internal matrix.
 * \ingroup data_structure
 *
 * \details
 *
 * This type is never used directly, instead use seqan3::nucleotide_scoring_scheme or seqan3::aminoacid_scoring_scheme.
 */
template <typename derived_t, alphabet_concept alphabet_t, arithmetic_concept score_type>
class scoring_scheme_base
{
public:
    /*!\name Member types
     * \{
     */
    //!\brief Type of the score values.
    using score_t = score_type;
    //!\brief Size type that can hold the dimension of the matrix (i.e. size of the alphabet).
    using matrix_size_t = uint8_t;
    //!\}

    //!\brief Size of the matrix dimensions (i.e. size of the alphabet).
    static constexpr matrix_size_t matrix_size = alphabet_size_v<alphabet_t>;

    /*!\name Member types
     * \{
     */
    //!\brief Type of the internal matrix (a two-dimensional array).
    using matrix_t = std::array<std::array<score_t, matrix_size>, matrix_size>;
    //!\}

private:
    //!\brief Befriend derived_t so it can instantiate.
    friend derived_t;

    //!\publicsection
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr scoring_scheme_base(scoring_scheme_base const &) noexcept = default;
    constexpr scoring_scheme_base(scoring_scheme_base &&) noexcept = default;
    constexpr scoring_scheme_base & operator=(scoring_scheme_base const &) noexcept = default;
    constexpr scoring_scheme_base & operator=(scoring_scheme_base &&) noexcept = default;
    ~scoring_scheme_base() noexcept = default;

    //!\brief The default constructor (delegates to set_hamming_distance()).
    constexpr scoring_scheme_base() noexcept
    {
        set_hamming_distance();
    }

    /*!\brief Constructor for the simple scheme (delegates to set_simple_scheme()).
     * \copydetails set_simple_scheme()
     */
    template <arithmetic_concept score_arg_t>
    constexpr scoring_scheme_base(match_score<score_arg_t> const ms, mismatch_score<score_arg_t> const mms)
    {
        set_simple_scheme(ms, mms);
    }

    /*!\brief Constructor for a custom scheme (delegates to set_custom_matrix()).
     * \copydetails set_custom_matrix()
     */
    constexpr scoring_scheme_base(matrix_t const & _matrix) noexcept
    {
        set_custom_matrix(_matrix);
    }
    //!\}

public:
    /*!\name Scheme selection
     * \{
     */
    //!\brief Set the hamming scheme, a variant of the simple scheme where match is scored `0` and mismatch `-1`.
    constexpr void set_hamming_distance() noexcept
    {
        set_simple_scheme(match_score<score_t>{0}, mismatch_score<score_t>{-1});
    }

    /*!\brief Set the simple scheme (everything is either match or mismatch).
     * \tparam score_arg_t The underlying type of the arguments.
     * \param[in] ms  Matches shall be given this value (of type seqan3::match_score).
     * \param[in] mms Mismatches shall be given this value (of type seqan3::mismatch_score).
     * \throws std::invalid_argument Thrown if you pass a value that is to large/low to be represented by `score_t`.
     */
    template <arithmetic_concept score_arg_t>
    constexpr void set_simple_scheme(match_score<score_arg_t> const ms, mismatch_score<score_arg_t> const mms)
    {
        std::conditional_t<integral_concept<score_t>, int64_t, double> i_ms = static_cast<score_arg_t>(ms);
        std::conditional_t<integral_concept<score_t>, int64_t, double> i_mms = static_cast<score_arg_t>(mms);
        if ((i_ms  < std::numeric_limits<score_t>::lowest() || i_ms  > std::numeric_limits<score_t>::max()) ||
            (i_mms < std::numeric_limits<score_t>::lowest() || i_mms > std::numeric_limits<score_t>::max()))
        {
            throw std::invalid_argument{"You passed a score value to set_simple_scheme that is out of range of the "
                                        "scoring scheme's underlying type. Define your scoring scheme with a larger "
                                        "template parameter or down-cast you score value beforehand to prevent "
                                        "this exception."};
        }

        for (matrix_size_t i = 0; i < matrix_size; ++i)
            for (matrix_size_t j = 0; j < matrix_size; ++j)
                matrix[i][j] = (i == j) ? static_cast<score_t>(i_ms) : static_cast<score_t>(i_mms);
    }

    /*!\brief Set a custom scheme by passing a full matrix with arbitrary content.
     * \param[in] _matrix A full matrix that is copied into the scheme.
     */
    constexpr void set_custom_matrix(matrix_t const & _matrix) noexcept
    {
        ranges::copy(_matrix, ranges::begin(matrix));
    }
    //!\}

    /*!\name Accessors
     * \{
     */
    /*!\brief Score two letters (either two nucleotids or two amino acids).
     * \tparam    alph1_t Type of the first letter.
     * \tparam    alph2_t Type of the second letter (needn't be the same as alph1_t).
     * \param[in] alph1   The first letter to score.
     * \param[in] alph2   The second letter to score.
     * \return The score of the two letters in the current scheme.
     */
    template <explicitly_convertible_to_concept<alphabet_t> alph1_t,
              explicitly_convertible_to_concept<alphabet_t> alph2_t>
    constexpr score_t & get_score(alph1_t const alph1, alph2_t const alph2) noexcept
    {
        return matrix[to_rank(static_cast<alphabet_t>(alph1))][to_rank(static_cast<alphabet_t>(alph2))];
    }

    //!\copydoc get_score
    template <explicitly_convertible_to_concept<alphabet_t> alph1_t,
              explicitly_convertible_to_concept<alphabet_t> alph2_t>
    constexpr score_t get_score(alph1_t const alph1, alph2_t const alph2) const noexcept
    {
        return matrix[to_rank(static_cast<alphabet_t>(alph1))][to_rank(static_cast<alphabet_t>(alph2))];
    }
    //!\}

    //!\name Comparison operators
    //!\{
    constexpr bool operator==(derived_t const & rhs) const noexcept
    {
        return matrix == rhs.matrix;
    }

    constexpr bool operator!=(derived_t const & rhs) const noexcept
    {
        return matrix != rhs.matrix;
    }
    //!\}

    /*!\cond DEV
     * \brief Serialisation support function.
     * \tparam archive_t Type of `archive`; must satisfy seqan3::cereal_archive_concept.
     * \param  archive   The archive being serialised from/to.
     *
     * \attention These functions are never called directly, see \ref serialisation for more details.
     */
    template <cereal_archive_concept archive_t>
    void CEREAL_SERIALIZE_FUNCTION_NAME(archive_t & archive)
    {
        archive(matrix);
    }
    //!\endcond

private:
    //!\brief The actual data member.
    matrix_t matrix;
};

} // namespace seqan3
