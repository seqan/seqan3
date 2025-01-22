// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <iostream>

#include <seqan3/alignment/matrix/detail/debug_matrix.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using seqan3::detail::debug_matrix;
    using namespace seqan3::literals;

    std::vector<seqan3::dna4> database = "AACCGGTT"_dna4;
    std::vector<seqan3::dna4> query = "ACGT"_dna4;

    seqan3::detail::row_wise_matrix<int> score_matrix{
        seqan3::detail::number_rows{5u},
        seqan3::detail::number_cols{9u},
        std::vector{-0, -1, -2, -3, -4, -5, -6, -7, -8, -1, -0, -1, -2, -3, -4, -5, -6, -7, -2, -1, -1, -1, -2,
                    -3, -4, -5, -6, -3, -2, -2, -2, -2, -2, -3, -4, -5, -4, -3, -3, -3, -3, -3, -3, -3, -4}};

    seqan3::debug_stream << "database:\t" << database << '\n';
    seqan3::debug_stream << "query:\t\t" << query << '\n';
    seqan3::debug_stream << '\n';

    seqan3::debug_stream << "score_matrix: " << score_matrix.cols() << " columns and " << score_matrix.rows()
                         << " rows\n";

    // Prints out the matrix in a convenient way
    seqan3::debug_stream << score_matrix << '\n';                                                   // without sequences
    seqan3::debug_stream << debug_matrix{score_matrix, database, query} << '\n';                    // with sequences
    seqan3::debug_stream << seqan3::fmtflags2::utf8 << debug_matrix{score_matrix, database, query}; // as utf8

    return 0;
}
