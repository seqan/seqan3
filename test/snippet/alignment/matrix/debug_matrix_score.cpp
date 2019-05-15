#include <iostream>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alignment/matrix/debug_matrix.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/to_char.hpp>

int main()
{
    using namespace seqan3;
    using seqan3::detail::debug_matrix;

    std::vector<dna4> database = "AACCGGTT"_dna4;
    std::vector<dna4> query = "ACGT"_dna4;

    detail::row_wise_matrix score_matrix
    {
        std::vector
        {
            -0, -1, -2, -3, -4, -5, -6, -7, -8,
            -1, -0, -1, -2, -3, -4, -5, -6, -7,
            -2, -1, -1, -1, -2, -3, -4, -5, -6,
            -3, -2, -2, -2, -2, -2, -3, -4, -5,
            -4, -3, -3, -3, -3, -3, -3, -3, -4
        },
        5u,
        9u
    };

    debug_stream << "database:\t" << database << std::endl;
    debug_stream << "query:\t\t" << query << std::endl;
    debug_stream << std::endl;

    debug_stream << "score_matrix: " << score_matrix.cols() << " columns and "
              << score_matrix.rows() << " rows" << std::endl;

    // Prints out the matrix in a convenient way
    debug_stream << score_matrix << std::endl; // without sequences
    debug_stream << debug_matrix{score_matrix, database, query} << std::endl; // with sequences
    debug_stream << fmtflags2::utf8 << debug_matrix{score_matrix, database, query}; // as utf8

    return 0;
}
