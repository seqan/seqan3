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
 * \brief Provides seqan3::nucleotide_scoring_scheme.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/data_structure/scoring_scheme_base.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>

namespace seqan3
{

/*!\brief A data structure for managing and computing the score of two nucleotides.
 * \tparam score_type The underlying type.
 * \ingroup data_structure
 *
 * \details
 *
 * You can use an instance of this class to score two nucleotides, the nucleotides need not be of the same type.
 * Different scoring behaviour can be set via the member functions.
 *
 * ### Example
 *
 * Score two letters:
 * ```cpp
 * nucleotide_scoring_scheme scheme; // hamming is default
 * std::cout << scheme.get_score(dna5::A, dna5::G); // == -1
 * std::cout << scheme.get_score(dna5::A, dna5::A); // == 0
 *
 * scheme.set_scheme_simple(match_score{3}, mismatch_score{-2});
 * std::cout << scheme.get_score(dna5::A, rna15::G); // == -2
 * std::cout << scheme.get_score(dna5::A, rna15::A); // == 3
 * // you can score differenct nucleotides  ^
 * ```
 *
 * You can "edit" a given matrix directly:
 * ```cpp
 * scheme.set_scheme_similarity_matrix;
 * std::cout << scheme.get_score(aa27::T, aa27::S); // == 2
 * scheme.get_score(aa27::T, aa27::S) = 3;
 * std::cout << scheme.get_score(aa27::T, aa27::S); // == 3
 * ```
 *
 * Score two sequences:
 * ```cpp
 * std::vector<dna15> one = "AGAATA"_dna5;
 * std::vector<dna15> two = "ATACTA"_dna5;
 *
 * nucleotide_scoring_scheme scheme; // haming distance is default
 * int score = 0;
 * for (auto pair : ranges::view::zip(one, two))
 *     score += scheme.get_score(std::get<0>(pair), std::get<1>(pair));
 * std::cout << score; // == 1 - 1 + 1 - 1 + 1 + 1 == 2
 *
 * ```
 */
template <arithmetic_concept score_type = int8_t>
class nucleotide_scoring_scheme : public scoring_scheme_base<nucleotide_scoring_scheme<score_type>, dna15, score_type>
{
private:
    //!\brief Type of the CRTP-base.
    using base_t = scoring_scheme_base<nucleotide_scoring_scheme, dna15, score_type>;
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    //NOTE(h-2): unfortunately these are not inherited in the documentation
    //!\copydoc scoring_scheme_base::scoring_scheme_base()
    SEQAN3_DOXYGEN_ONLY(( constexpr nucleotide_scoring_scheme() noexcept {} ))
    //!\copydoc scoring_scheme_base::scoring_scheme_base(match_score<score_arg_t> const ms, mismatch_score<score_arg_t> const mms)
    SEQAN3_DOXYGEN_ONLY((
      template <arithmetic_concept score_arg_t>
      constexpr nucleotide_scoring_scheme(match_score<score_arg_t> const ms, mismatch_score<score_arg_t> const mms) {}
    ))
    //!\copydoc scoring_scheme_base::scoring_scheme_base(matrix_t const & _matrix)
    SEQAN3_DOXYGEN_ONLY(( constexpr nucleotide_scoring_scheme(matrix_t const & _matrix) noexcept {} ))
    //!\}

    //!\privatesection
    //!\brief Inherit the base class's constructors.
    using base_t::base_t;
};

/*!\name Type deduction guides
 * \relates seqan3::nucleotide_scoring_scheme
 * \{
 */
nucleotide_scoring_scheme() -> nucleotide_scoring_scheme<int8_t>;

/*!\brief Attention: This guide does not actually deduce from the underlying type, but always defaults to `int8_t`.
 * To use a larger type, specify the template argument manually.
 */
template <arithmetic_concept score_arg_type>
nucleotide_scoring_scheme(match_score<score_arg_type>,
                          mismatch_score<score_arg_type>) -> nucleotide_scoring_scheme<int8_t>;

template <arithmetic_concept score_arg_type>
nucleotide_scoring_scheme(std::array<std::array<score_arg_type, 15>, 15>) -> nucleotide_scoring_scheme<score_arg_type>;
//!\}

} // namespace seqan3
