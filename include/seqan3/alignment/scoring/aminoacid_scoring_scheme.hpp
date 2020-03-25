// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::aminoacid_scoring_scheme.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <range/v3/algorithm/copy.hpp>

#include <seqan3/alignment/scoring/scoring_scheme_base.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/std/algorithm>

namespace seqan3
{

/*!\brief Identifiers for amino acid similarity matrixes.
 * \ingroup scoring
 *
 * \details
 *
 * This enum provides IDs for amino acid similarity matrixes of the
 * [BLOSUM](https://en.wikipedia.org/wiki/BLOSUM) and [PAM](https://en.wikipedia.org/wiki/Point_accepted_mutation)
 * families.
 *
 * \see seqan3::aminoacid_scoring_scheme::set_similarity_matrix
 */
enum class aminoacid_similarity_matrix
{
    //ATTENTION: when you change this, also update set_similarity_matrix() below
    /*!\brief The BLOSUM30 matrix for very distantly related proteins.
     * \details
     * \snippet include/seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp matrix_data blosum30
     */
    BLOSUM30,
    /*!\brief The BLOSUM45 matrix for distantly related proteins.
     * \details
     * \snippet include/seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp matrix_data blosum45
     */
    BLOSUM45,
    /*!\brief The BLOSUM62 matrix recommended for most use-cases.
     * \details
     * \snippet include/seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp matrix_data blosum62
     */
    BLOSUM62,
    /*!\brief The BLOSUM80 matrix for closely related proteins.
     * \details
     * \snippet include/seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp matrix_data blosum80
     */
    BLOSUM80
};

/*!\brief A data structure for managing and computing the score of two amino acids.
 * \ingroup scoring
 * \implements seqan3::scoring_scheme
 *
 * \details
 *
 * You can use an instance of this class to score two amino acids. The amino acids need not be of the same type.
 * Different scoring behaviour can be set via the member functions.
 *
 * ### Example
 *
 * \include test/snippet/alignment/scoring/aminoacid_scoring_scheme.cpp
 */
template <arithmetic score_type = int8_t>
class aminoacid_scoring_scheme : public scoring_scheme_base<aminoacid_scoring_scheme<score_type>, aa27, score_type>
{
private:
    //!\brief Type of the CRTP-base.
    using base_t = scoring_scheme_base<aminoacid_scoring_scheme<score_type>, aa27, score_type>;

    //!\brief Make the base's private member visible.
    using base_t::matrix;

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
    constexpr aminoacid_scoring_scheme() noexcept = default;
    //!\copydoc scoring_scheme_base::scoring_scheme_base(match_score<score_arg_t> const ms, mismatch_score<score_arg_t> const mms)
    template <arithmetic score_arg_t>
    constexpr aminoacid_scoring_scheme(match_score<score_arg_t> const ms, mismatch_score<score_arg_t> const mms)
        : base_t{ms, mms}
    {}
    //!\copydoc scoring_scheme_base::scoring_scheme_base(matrix_type const & matrix)
    constexpr aminoacid_scoring_scheme(matrix_type const & matrix) noexcept
        : base_t{matrix}
    {}

    //!\brief Construct for seqan3::aminoacid_similarity_matrix.
    //!\copydetails set_similarity_matrix()
    constexpr aminoacid_scoring_scheme(aminoacid_similarity_matrix const matrix_id)
    {
        set_similarity_matrix(matrix_id);
    }
    //!\}

    /*!\name Scheme selection
     * \{
     */
    /*!\brief Set the similarity matrix scheme (e.g. BLOSUM62).
     * \param[in] matrix_id The enum ID of the matrix, see seqan3::aminoacid_similarity_matrix.
     * \throws std::invalid_argument If there is no matrix data for the given ID (usually a BUG).
     */
    constexpr void set_similarity_matrix(aminoacid_similarity_matrix const matrix_id)
    {
        switch (matrix_id)
        {
            case aminoacid_similarity_matrix::BLOSUM30: std::ranges::copy(blosum30, begin(matrix)); break;
            case aminoacid_similarity_matrix::BLOSUM45: std::ranges::copy(blosum45, begin(matrix)); break;
            case aminoacid_similarity_matrix::BLOSUM62: std::ranges::copy(blosum62, begin(matrix)); break;
            case aminoacid_similarity_matrix::BLOSUM80: std::ranges::copy(blosum80, begin(matrix)); break;
            default:
                throw std::invalid_argument{"ERROR in set_similarity_matrix(), matrix_id has no matrix."};
        }
    }
    //!\}

private:
    //!\brief The matrix data corresponding to seqan3::aminoacid_similarity_matrix::BLOSUM30.
    static constexpr matrix_type blosum30
    {{
        //! [matrix_data blosum30]
        //A   B   C   D   E   F   G   H   I   J   K   L   M   N   O   P   Q   R   S   T   U   V   W   X,  Y   Z   *
        { 4,  0, -3,  0,  0, -2,  0, -2,  0, -1,  0, -1,  1,  0,  0, -1,  1, -1,  1,  1,  0,  1, -5,  0, -4,  0, -7},//A
        { 0,  5, -2,  5,  0, -3,  0, -2, -2, -2,  0, -1, -2,  4, -1, -2, -1, -2,  0,  0, -1, -2, -5, -1, -3,  0, -7},//B
        {-3, -2, 17, -3,  1, -3, -4, -5, -2, -1, -3,  0, -2, -1, -2, -3, -2, -2, -2, -2, -2, -2, -2, -2, -6,  0, -7},//C
        { 0,  5, -3,  9,  1, -5, -1, -2, -4, -3,  0, -1, -3,  1, -1, -1, -1, -1,  0, -1, -1, -2, -4, -1, -1,  0, -7},//D
        { 0,  0,  1,  1,  6, -4, -2,  0, -3, -2,  2, -1, -1, -1, -1,  1,  2, -1,  0, -2, -1, -3, -1, -1, -2,  5, -7},//E
        {-2, -3, -3, -5, -4, 10, -3, -3,  0,  1, -1,  2, -2, -1, -1, -4, -3, -1, -1, -2, -1,  1,  1, -1,  3, -4, -7},//F
        { 0,  0, -4, -1, -2, -3,  8, -3, -1, -2, -1, -2, -2,  0, -1, -1, -2, -2,  0, -2, -1, -3,  1, -1, -3, -2, -7},//G
        {-2, -2, -5, -2,  0, -3, -3, 14, -2, -2, -2, -1,  2, -1, -1,  1,  0, -1, -1, -2, -1, -3, -5, -1,  0,  0, -7},//H
        { 0, -2, -2, -4, -3,  0, -1, -2,  6,  4, -2,  2,  1,  0,  0, -3, -2, -3, -1,  0,  0,  4, -3,  0, -1, -3, -7},//I
        {-1, -2, -1, -3, -2,  1, -2, -2,  4,  4, -2,  3,  2, -1,  0, -3, -2, -3, -2,  0,  0,  3, -3,  0,  1, -2, -7},//J
        { 0,  0, -3,  0,  2, -1, -1, -2, -2, -2,  4, -2,  2,  0,  0,  1,  0,  1,  0, -1,  0, -2, -2,  0, -1,  1, -7},//K
        {-1, -1,  0, -1, -1,  2, -2, -1,  2,  3, -2,  4,  2, -2,  0, -3, -2, -2, -2,  0,  0,  1, -2,  0,  3, -1, -7},//L
        { 1, -2, -2, -3, -1, -2, -2,  2,  1,  2,  2,  2,  6,  0,  0, -4, -1,  0, -2,  0,  0,  0, -3,  0, -1, -1, -7},//M
        { 0,  4, -1,  1, -1, -1,  0, -1,  0, -1,  0, -2,  0,  8,  0, -3, -1, -2,  0,  1,  0, -2, -7,  0, -4, -1, -7},//N
        { 0, -1, -2, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0, -1, -1,  0, -1,  0,  0, -1,  0, -2, -1, -1,  0, -7},//O
        {-1, -2, -3, -1,  1, -4, -1,  1, -3, -3,  1, -3, -4, -3, -1, 11,  0, -1, -1,  0, -1, -4, -3, -1, -2,  0, -7},//P
        { 1, -1, -2, -1,  2, -3, -2,  0, -2, -2,  0, -2, -1, -1,  0,  0,  8,  3, -1,  0,  0, -3, -1,  0, -1,  4, -7},//Q
        {-1, -2, -2, -1, -1, -1, -2, -1, -3, -3,  1, -2,  0, -2, -1, -1,  3,  8, -1, -3, -1, -1,  0, -1,  0,  0, -7},//R
        { 1,  0, -2,  0,  0, -1,  0, -1, -1, -2,  0, -2, -2,  0,  0, -1, -1, -1,  4,  2,  0, -1, -3,  0, -2, -1, -7},//S
        { 1,  0, -2, -1, -2, -2, -2, -2,  0,  0, -1,  0,  0,  1,  0,  0,  0, -3,  2,  5,  0,  1, -5,  0, -1, -1, -7},//T
        { 0, -1, -2, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0, -1, -1,  0, -1,  0,  0, -1,  0, -2, -1, -1,  0, -7},//U
        { 1, -2, -2, -2, -3,  1, -3, -3,  4,  3, -2,  1,  0, -2,  0, -4, -3, -1, -1,  1,  0,  5, -3,  0,  1, -3, -7},//V
        {-5, -5, -2, -4, -1,  1,  1, -5, -3, -3, -2, -2, -3, -7, -2, -3, -1,  0, -3, -5, -2, -3, 20, -2,  5, -1, -7},//W
        { 0, -1, -2, -1, -1, -1, -1, -1,  0,  0,  0,  0,  0,  0, -1, -1,  0, -1,  0,  0, -1,  0, -2, -1, -1,  0, -7},//X
        {-4, -3, -6, -1, -2,  3, -3,  0, -1,  1, -1,  3, -1, -4, -1, -2, -1,  0, -2, -1, -1,  1,  5, -1,  9, -2, -7},//Y
        { 0,  0,  0,  0,  5, -4, -2,  0, -3, -2,  1, -1, -1, -1,  0,  0,  4,  0, -1, -1,  0, -3, -1,  0, -2,  4, -7},//Z
        {-7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7,  1} //*
        //! [matrix_data blosum30]
    }};

    //!\brief The matrix data corresponding to seqan3::aminoacid_similarity_matrix::BLOSUM45.
    static constexpr matrix_type blosum45
    {{
        //! [matrix_data blosum45]
        //A   B   C   D   E   F   G   H   I   J   K   L   M   N   O   P   Q   R   S   T   U   V   W   X,  Y   Z   *
        { 5, -1, -1, -2, -1, -2,  0, -2, -1, -1, -1, -1, -1, -1,  0, -1, -1, -2,  1,  0,  0,  0, -2,  0, -2, -1, -5},//A
        {-1,  4, -2,  5,  1, -3, -1,  0, -3, -3,  0, -3, -2,  4, -1, -2,  0, -1,  0,  0, -1, -3, -4, -1, -2,  2, -5},//B
        {-1, -2, 12, -3, -3, -2, -3, -3, -3, -3, -3, -2, -2, -2, -2, -4, -3, -3, -1, -1, -2, -1, -5, -2, -3, -3, -5},//C
        {-2,  5, -3,  7,  2, -4, -1,  0, -4, -4,  0, -3, -3,  2, -1, -1,  0, -1,  0, -1, -1, -3, -4, -1, -2,  1, -5},//D
        {-1,  1, -3,  2,  6, -3, -2,  0, -3, -3,  1, -2, -2,  0, -1,  0,  2,  0,  0, -1, -1, -3, -3, -1, -2,  4, -5},//E
        {-2, -3, -2, -4, -3,  8, -3, -2,  0,  1, -3,  1,  0, -2, -1, -3, -4, -2, -2, -1, -1,  0,  1, -1,  3, -3, -5},//F
        { 0, -1, -3, -1, -2, -3,  7, -2, -4, -4, -2, -3, -2,  0, -1, -2, -2, -2,  0, -2, -1, -3, -2, -1, -3, -2, -5},//G
        {-2,  0, -3,  0,  0, -2, -2, 10, -3, -3, -1, -2,  0,  1, -1, -2,  1,  0, -1, -2, -1, -3, -3, -1,  2,  0, -5},//H
        {-1, -3, -3, -4, -3,  0, -4, -3,  5,  4, -3,  2,  2, -2, -1, -2, -2, -3, -2, -1, -1,  3, -2, -1,  0, -3, -5},//I
        {-1, -3, -3, -4, -3,  1, -4, -3,  4,  4, -3,  4,  2, -3, -1, -3, -2, -3, -3, -1, -1,  2, -2, -1,  0, -3, -5},//J
        {-1,  0, -3,  0,  1, -3, -2, -1, -3, -3,  5, -3, -1,  0, -1, -1,  1,  3, -1, -1, -1, -2, -2, -1, -1,  1, -5},//K
        {-1, -3, -2, -3, -2,  1, -3, -2,  2,  4, -3,  5,  2, -3, -1, -3, -2, -2, -3, -1, -1,  1, -2, -1,  0, -2, -5},//L
        {-1, -2, -2, -3, -2,  0, -2,  0,  2,  2, -1,  2,  6, -2, -1, -2,  0, -1, -2, -1, -1,  1, -2, -1,  0, -1, -5},//M
        {-1,  4, -2,  2,  0, -2,  0,  1, -2, -3,  0, -3, -2,  6, -1, -2,  0,  0,  1,  0, -1, -3, -4, -1, -2,  0, -5},//N
        { 0, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0, -1, -1, -2, -1, -1, -1, -5},//O
        {-1, -2, -4, -1,  0, -3, -2, -2, -2, -3, -1, -3, -2, -2, -1,  9, -1, -2, -1, -1, -1, -3, -3, -1, -3, -1, -5},//P
        {-1,  0, -3,  0,  2, -4, -2,  1, -2, -2,  1, -2,  0,  0, -1, -1,  6,  1,  0, -1, -1, -3, -2, -1, -1,  4, -5},//Q
        {-2, -1, -3, -1,  0, -2, -2,  0, -3, -3,  3, -2, -1,  0, -1, -2,  1,  7, -1, -1, -1, -2, -2, -1, -1,  0, -5},//R
        { 1,  0, -1,  0,  0, -2,  0, -1, -2, -3, -1, -3, -2,  1,  0, -1,  0, -1,  4,  2,  0, -1, -4,  0, -2,  0, -5},//S
        { 0,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -1,  0,  0, -1, -1, -1,  2,  5,  0,  0, -3,  0, -1, -1, -5},//T
        { 0, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0, -1, -1, -2, -1, -1, -1, -5},//U
        { 0, -3, -1, -3, -3,  0, -3, -3,  3,  2, -2,  1,  1, -3, -1, -3, -3, -2, -1,  0, -1,  5, -3, -1, -1, -3, -5},//V
        {-2, -4, -5, -4, -3,  1, -2, -3, -2, -2, -2, -2, -2, -4, -2, -3, -2, -2, -4, -3, -2, -3, 15, -2,  3, -2, -5},//W
        { 0, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0, -1, -1, -2, -1, -1, -1, -5},//X
        {-2, -2, -3, -2, -2,  3, -3,  2,  0,  0, -1,  0,  0, -2, -1, -3, -1, -1, -2, -1, -1, -1,  3, -1,  8, -2, -5},//Y
        {-1,  2, -3,  1,  4, -3, -2,  0, -3, -3,  1, -2, -1,  0, -1, -1,  4,  0,  0, -1, -1, -3, -2, -1, -2,  4, -5},//Z
        {-5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  1} //*
        //! [matrix_data blosum45]
    }};

    //!\brief The matrix data corresponding to seqan3::aminoacid_similarity_matrix::BLOSUM62.
    static constexpr matrix_type blosum62
    {{
        //! [matrix_data blosum62]
        //A   B   C   D   E   F   G   H   I   J   K   L   M   N   O   P   Q   R   S   T   U   V   W   X   Y   Z   *
        { 4, -2,  0, -2, -1, -2,  0, -2, -1, -1, -1, -1, -1, -2,  0, -1, -1, -1,  1,  0,  0,  0, -3,  0, -2, -1, -4},//A
        {-2,  4, -3,  4,  1, -3, -1,  0, -3, -4,  0, -4, -3,  3, -1, -2,  0, -1,  0, -1, -1, -3, -4, -1, -3,  1, -4},//B
        { 0, -3,  9, -3, -4, -2, -3, -3, -1, -1, -3, -1, -1, -3, -2, -3, -3, -3, -1, -1, -2, -1, -2, -2, -2, -3, -4},//C
        {-2,  4, -3,  6,  2, -3, -1, -1, -3, -4, -1, -4, -3,  1, -1, -1,  0, -2,  0, -1, -1, -3, -4, -1, -3,  1, -4},//D
        {-1,  1, -4,  2,  5, -3, -2,  0, -3, -3,  1, -3, -2,  0, -1, -1,  2,  0,  0, -1, -1, -2, -3, -1, -2,  4, -4},//E
        {-2, -3, -2, -3, -3,  6, -3, -1,  0,  0, -3,  0,  0, -3, -1, -4, -3, -3, -2, -2, -1, -1,  1, -1,  3, -3, -4},//F
        { 0, -1, -3, -1, -2, -3,  6, -2, -4, -4, -2, -4, -3,  0, -1, -2, -2, -2,  0, -2, -1, -3, -2, -1, -3, -2, -4},//G
        {-2,  0, -3, -1,  0, -1, -2,  8, -3, -3, -1, -3, -2,  1, -1, -2,  0,  0, -1, -2, -1, -3, -2, -1,  2,  0, -4},//H
        {-1, -3, -1, -3, -3,  0, -4, -3,  4,  3, -3,  2,  1, -3, -1, -3, -3, -3, -2, -1, -1,  3, -3, -1, -1, -3, -4},//I
        {-1, -4, -1, -4, -3,  0, -4, -3,  3,  3, -3,  3,  2, -3, -1, -3, -3, -3, -2, -1, -1,  2, -3, -1, -1, -3, -4},//J
        {-1,  0, -3, -1,  1, -3, -2, -1, -3, -3,  5, -2, -1,  0, -1, -1,  1,  2,  0, -1, -1, -2, -3, -1, -2,  1, -4},//K
        {-1, -4, -1, -4, -3,  0, -4, -3,  2,  3, -2,  4,  2, -3, -1, -3, -2, -2, -2, -1, -1,  1, -2, -1, -1, -3, -4},//L
        {-1, -3, -1, -3, -2,  0, -3, -2,  1,  2, -1,  2,  5, -2, -1, -2,  0, -1, -1, -1, -1,  1, -1, -1, -1, -1, -4},//M
        {-2,  3, -3,  1,  0, -3,  0,  1, -3, -3,  0, -3, -2,  6, -1, -2,  0,  0,  1,  0, -1, -3, -4, -1, -2,  0, -4},//N
        { 0, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, -1, -1,  0,  0, -1, -1, -2, -1, -1, -1, -4},//O
        {-1, -2, -3, -1, -1, -4, -2, -2, -3, -3, -1, -3, -2, -2, -2,  7, -1, -2, -1, -1, -2, -2, -4, -2, -3, -1, -4},//P
        {-1,  0, -3,  0,  2, -3, -2,  0, -3, -3,  1, -2,  0,  0, -1, -1,  5,  1,  0, -1, -1, -2, -2, -1, -1,  3, -4},//Q
        {-1, -1, -3, -2,  0, -3, -2,  0, -3, -3,  2, -2, -1,  0, -1, -2,  1,  5, -1, -1, -1, -3, -3, -1, -2,  0, -4},//R
        { 1,  0, -1,  0,  0, -2,  0, -1, -2, -2,  0, -2, -1,  1,  0, -1,  0, -1,  4,  1,  0, -2, -3,  0, -2,  0, -4},//S
        { 0, -1, -1, -1, -1, -2, -2, -2, -1, -1, -1, -1, -1,  0,  0, -1, -1, -1,  1,  5,  0,  0, -2,  0, -2, -1, -4},//T
        { 0, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, -1, -1,  0,  0, -1, -1, -2, -1, -1, -1, -4},//U
        { 0, -3, -1, -3, -2, -1, -3, -3,  3,  2, -2,  1,  1, -3, -1, -2, -2, -3, -2,  0, -1,  4, -3, -1, -1, -2, -4},//V
        {-3, -4, -2, -4, -3,  1, -2, -2, -3, -3, -3, -2, -1, -4, -2, -4, -2, -3, -3, -2, -2, -3, 11, -2,  2, -3, -4},//W
        { 0, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2, -1, -1,  0,  0, -1, -1, -2, -1, -1, -1, -4},//X
        {-2, -3, -2, -3, -2,  3, -3,  2, -1, -1, -2, -1, -1, -2, -1, -3, -1, -2, -2, -2, -1, -1,  2, -1,  7, -2, -4},//Y
        {-1,  1, -3,  1,  4, -3, -2,  0, -3, -3,  1, -3, -1,  0, -1, -1,  3,  0,  0, -1, -1, -2, -3, -1, -2,  4, -4},//Z
        {-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1} //*
        //! [matrix_data blosum62]
    }};

    //!\brief The matrix data corresponding to seqan3::aminoacid_similarity_matrix::BLOSUM80.
    static constexpr matrix_type blosum80
    {{
        //! [matrix_data blosum80]
        //A   B   C   D   E   F   G   H   I   J   K   L   M   N   O   P   Q   R   S   T   U   V   W   X,  Y   Z   *
        { 7, -3, -1, -3, -2, -4,  0, -3, -3, -3, -1, -3, -2, -3, -1, -1, -2, -3,  2,  0, -1, -1, -5, -1, -4, -2, -8},//A
        {-3,  6, -6,  6,  1, -6, -2, -1, -6, -7, -1, -7, -5,  5, -3, -4, -1, -2,  0, -1, -3, -6, -8, -3, -5,  0, -8},//B
        {-1, -6, 13, -7, -7, -4, -6, -7, -2, -3, -6, -3, -3, -5, -4, -6, -5, -6, -2, -2, -4, -2, -5, -4, -5, -7, -8},//C
        {-3,  6, -7, 10,  2, -6, -3, -2, -7, -7, -2, -7, -6,  2, -3, -3, -1, -3, -1, -2, -3, -6, -8, -3, -6,  1, -8},//D
        {-2,  1, -7,  2,  8, -6, -4,  0, -6, -6,  1, -6, -4, -1, -2, -2,  3, -1, -1, -2, -2, -4, -6, -2, -5,  6, -8},//E
        {-4, -6, -4, -6, -6, 10, -6, -2, -1, -1, -5,  0,  0, -6, -3, -6, -5, -5, -4, -4, -3, -2,  0, -3,  4, -6, -8},//F
        { 0, -2, -6, -3, -4, -6,  9, -4, -7, -7, -3, -7, -5, -1, -3, -5, -4, -4, -1, -3, -3, -6, -6, -3, -6, -4, -8},//G
        {-3, -1, -7, -2,  0, -2, -4, 12, -6, -6, -1, -5, -4,  1, -2, -4,  1,  0, -2, -3, -2, -5, -4, -2,  3,  0, -8},//H
        {-3, -6, -2, -7, -6, -1, -7, -6,  7,  5, -5,  2,  2, -6, -2, -5, -5, -5, -4, -2, -2,  4, -5, -2, -3, -6, -8},//I
        {-3, -7, -3, -7, -6, -1, -7, -6,  5,  5, -5,  4,  3, -6, -2, -5, -5, -5, -4, -3, -2,  3, -5, -2, -3, -6, -8},//J
        {-1, -1, -6, -2,  1, -5, -3, -1, -5, -5,  8, -4, -3,  0, -2, -2,  2,  3, -1, -1, -2, -4, -6, -2, -4,  1, -8},//K
        {-3, -7, -3, -7, -6,  0, -7, -5,  2,  4, -4,  6,  3, -6, -2, -5, -4, -4, -4, -3, -2,  1, -4, -2, -2, -5, -8},//L
        {-2, -5, -3, -6, -4,  0, -5, -4,  2,  3, -3,  3,  9, -4, -2, -4, -1, -3, -3, -1, -2,  1, -3, -2, -3, -3, -8},//M
        {-3,  5, -5,  2, -1, -6, -1,  1, -6, -6,  0, -6, -4,  9, -2, -4,  0, -1,  1,  0, -2, -5, -7, -2, -4, -1, -8},//N
        {-1, -3, -4, -3, -2, -3, -3, -2, -2, -2, -2, -2, -2, -2, -2, -3, -2, -2, -1, -1, -2, -2, -5, -2, -3, -1, -8},//O
        {-1, -4, -6, -3, -2, -6, -5, -4, -5, -5, -2, -5, -4, -4, -3, 12, -3, -3, -2, -3, -3, -4, -7, -3, -6, -2, -8},//P
        {-2, -1, -5, -1,  3, -5, -4,  1, -5, -5,  2, -4, -1,  0, -2, -3,  9,  1, -1, -1, -2, -4, -4, -2, -3,  5, -8},//Q
        {-3, -2, -6, -3, -1, -5, -4,  0, -5, -5,  3, -4, -3, -1, -2, -3,  1,  9, -2, -2, -2, -4, -5, -2, -4,  0, -8},//R
        { 2,  0, -2, -1, -1, -4, -1, -2, -4, -4, -1, -4, -3,  1, -1, -2, -1, -2,  7,  2, -1, -3, -6, -1, -3, -1, -8},//S
        { 0, -1, -2, -2, -2, -4, -3, -3, -2, -3, -1, -3, -1,  0, -1, -3, -1, -2,  2,  8, -1,  0, -5, -1, -3, -2, -8},//T
        {-1, -3, -4, -3, -2, -3, -3, -2, -2, -2, -2, -2, -2, -2, -2, -3, -2, -2, -1, -1, -2, -2, -5, -2, -3, -1, -8},//U
        {-1, -6, -2, -6, -4, -2, -6, -5,  4,  3, -4,  1,  1, -5, -2, -4, -4, -4, -3,  0, -2,  7, -5, -2, -3, -4, -8},//V
        {-5, -8, -5, -8, -6,  0, -6, -4, -5, -5, -6, -4, -3, -7, -5, -7, -4, -5, -6, -5, -5, -5, 16, -5,  3, -5, -8},//W
        {-1, -3, -4, -3, -2, -3, -3, -2, -2, -2, -2, -2, -2, -2, -2, -3, -2, -2, -1, -1, -2, -2, -5, -2, -3, -1, -8},//X
        {-4, -5, -5, -6, -5,  4, -6,  3, -3, -3, -4, -2, -3, -4, -3, -6, -3, -4, -3, -3, -3, -3,  3, -3, 11, -4, -8},//Y
        {-2,  0, -7,  1,  6, -6, -4,  0, -6, -6,  1, -5, -3, -1, -1, -2,  5,  0, -1, -2, -1, -4, -5, -1, -4,  6, -8},//Z
        {-8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8,  1} //*
        //! [matrix_data blosum80]
    }};

    //TODO(h-2): add more matrixes
};

/*!\name Type deduction guides
 * \relates seqan3::aminoacid_scoring_scheme
 * \{
 */

//!\brief Default constructed objects deduce to `int8_t`.
aminoacid_scoring_scheme() -> aminoacid_scoring_scheme<int8_t>;

/*!\brief Attention: This guide does not actually deduce from the underlying type, but always defaults to `int8_t`.
 * To use a larger type, specify the template argument manually.
 */
template <arithmetic score_arg_type>
aminoacid_scoring_scheme(match_score<score_arg_type>,
                         mismatch_score<score_arg_type>) -> aminoacid_scoring_scheme<int8_t>;

//!\brief Deduce the score type from the provided matrix.
template <arithmetic score_arg_type>
aminoacid_scoring_scheme(std::array<std::array<score_arg_type, 27>, 27>) -> aminoacid_scoring_scheme<score_arg_type>;

/*!\brief Attention: This guide does not actually deduce from the underlying type, but always defaults to `int8_t`.
 * To use a larger type, specify the template argument manually.
 */
aminoacid_scoring_scheme(aminoacid_similarity_matrix) -> aminoacid_scoring_scheme<int8_t>;
//!\}

} // namespace seqan3
