// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

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
 * \ingroup alignment_scoring
 * \implements seqan3::scoring_scheme
 *
 * \details
 *
 * You can use an instance of this class to score two nucleotides, the nucleotides need not be of the same type.
 * Different scoring behaviour can be set via the member functions.
 *
 * ### Example
 *
 * \include test/snippet/alignment/scoring/nucleotide_scoring_scheme.cpp
 */
template <arithmetic score_type = int8_t>
class nucleotide_scoring_scheme : public scoring_scheme_base<nucleotide_scoring_scheme<score_type>, dna15, score_type>
{
private:
    //!\brief Type of the CRTP-base.
    using base_t = scoring_scheme_base<nucleotide_scoring_scheme, dna15, score_type>;

    //!\brief Befriend base_t so it can access itself through this derived type.
    friend base_t;

public:
    //!\privatesection
    //!\copydoc scoring_scheme_base::matrix_type
    using typename base_t::matrix_type;
    //!\publicsection

    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\copydoc scoring_scheme_base::scoring_scheme_base()
    constexpr nucleotide_scoring_scheme() noexcept = default;
    //!\copydoc scoring_scheme_base::scoring_scheme_base(match_score<score_arg_t> const ms, mismatch_score<score_arg_t> const mms)
    template <arithmetic score_arg_t>
    constexpr nucleotide_scoring_scheme(match_score<score_arg_t> const ms, mismatch_score<score_arg_t> const mms) :
        base_t{ms, mms}
    {}
    //!\copydoc scoring_scheme_base::scoring_scheme_base(matrix_type const & matrix)
    constexpr nucleotide_scoring_scheme(matrix_type const & matrix) noexcept : base_t{matrix}
    {}
    //!\}
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
template <arithmetic score_arg_type>
nucleotide_scoring_scheme(match_score<score_arg_type>,
                          mismatch_score<score_arg_type>) -> nucleotide_scoring_scheme<int8_t>;

//!\brief Deduce the score type from the provided matrix.
template <arithmetic score_arg_type>
nucleotide_scoring_scheme(std::array<std::array<score_arg_type, 15>, 15>) -> nucleotide_scoring_scheme<score_arg_type>;
//!\}

} // namespace seqan3
