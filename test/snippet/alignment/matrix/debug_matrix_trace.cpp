#include <iostream>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alignment/matrix/debug_matrix.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/views/to_char.hpp>

int main()
{
    using seqan3::detail::debug_matrix;
    using seqan3::operator""_dna4;
    using seqan3::operator|;

    std::vector<seqan3::dna4> database = "AACCGGTT"_dna4;
    std::vector<seqan3::dna4> query = "ACGT"_dna4;

    auto N = seqan3::detail::trace_directions::none;
    auto D = seqan3::detail::trace_directions::diagonal;
    auto U = seqan3::detail::trace_directions::up;
    auto L = seqan3::detail::trace_directions::left;

    seqan3::detail::row_wise_matrix<seqan3::detail::trace_directions> trace_matrix
    {
        seqan3::detail::number_rows{5u},
        seqan3::detail::number_cols{9u},
        std::vector
        {
            N,L,L  ,L  ,L  ,L  ,L  ,L,L  ,
            U,D,D|L,L  ,L  ,L  ,L  ,L,L  ,
            U,U,D  ,D  ,D|L,L  ,L  ,L,L  ,
            U,U,D|U,D|U,D  ,D  ,D|L,L,L  ,
            U,U,D|U,D|U,D|U,D|U,D  ,D,D|L
        }
    };

    seqan3::debug_stream << "database:\t" << database << '\n';
    seqan3::debug_stream << "query:\t\t" << query << '\n';
    seqan3::debug_stream << '\n';

    seqan3::debug_stream << "trace_matrix: " << trace_matrix.cols() << " columns and "
              << trace_matrix.rows() << " rows\n";

    // Prints out the matrix in a convenient way
    seqan3::debug_stream << trace_matrix << '\n'; // without sequences
    seqan3::debug_stream << debug_matrix{trace_matrix, database, query} << '\n'; // with sequences
    seqan3::debug_stream << seqan3::fmtflags2::utf8 << debug_matrix{trace_matrix, database, query}; // as utf8
    return 0;
}
