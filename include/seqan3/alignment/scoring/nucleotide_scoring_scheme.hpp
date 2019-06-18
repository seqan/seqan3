// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::nucleotide_scoring_scheme.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/scoring/scoring_scheme_base.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>

namespace seqan3
{

/*!\brief A data structure for managing and computing the score of two nucleotides.
 * \tparam score_type The underlying type.
 * \ingroup scoring
 * \implements seqan3::ScoringScheme
 *
 * \details
 *
 * You can use an instance of this class to score two nucleotides, the nucleotides need not be of the same type.
 * Different scoring behaviour can be set via the member functions.
 *
 * ### Example
 *
 * Score two letters:
 * \snippet test/snippet/alignment/scoring/nucleotide_scoring_scheme.cpp two letters
 *
 * You can "edit" a given matrix directly:
 * \snippet test/snippet/alignment/scoring/nucleotide_scoring_scheme.cpp edit matrix
 *
 * Score two sequences:
 * \snippet test/snippet/alignment/scoring/nucleotide_scoring_scheme.cpp score sequences
 */
template <Arithmetic score_type = int8_t>
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
      template <Arithmetic score_arg_t>
      constexpr nucleotide_scoring_scheme(match_score<score_arg_t> const ms, mismatch_score<score_arg_t> const mms) {}
    ))
    //!\copydoc scoring_scheme_base::scoring_scheme_base(matrix_type const & _matrix)
    SEQAN3_DOXYGEN_ONLY(( constexpr nucleotide_scoring_scheme(matrix_type const & _matrix) noexcept {} ))
    //!\}

    //!\privatesection
    //!\brief Inherit the base class's constructors.
    using base_t::base_t;
};

/*!\name Type deduction guides
 * \relates seqan3::nucleotide_scoring_scheme
 * \{
 */

//!\brief Default constructed objects deduce to `int8_t`.
nucleotide_scoring_scheme() -> nucleotide_scoring_scheme<int8_t>;

/*!\brief Attention: This guide does not actually deduce from the underlying type, but always defaults to `int8_t`.
 * To use a larger type, specify the template argument manually.
 */
template <Arithmetic score_arg_type>
nucleotide_scoring_scheme(match_score<score_arg_type>,
                          mismatch_score<score_arg_type>) -> nucleotide_scoring_scheme<int8_t>;

//!\brief Deduce the score type from the provided matrix.
template <Arithmetic score_arg_type>
nucleotide_scoring_scheme(std::array<std::array<score_arg_type, 15>, 15>) -> nucleotide_scoring_scheme<score_arg_type>;
//!\}

} // namespace seqan3
