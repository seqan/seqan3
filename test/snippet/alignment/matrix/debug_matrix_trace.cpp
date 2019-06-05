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

    auto N = detail::trace_directions::none;
    auto D = detail::trace_directions::diagonal;
    auto U = detail::trace_directions::up;
    auto L = detail::trace_directions::left;

    detail::row_wise_matrix trace_matrix
    {
        std::vector
        {
            N,L,L  ,L  ,L  ,L  ,L  ,L,L  ,
            U,D,D|L,L  ,L  ,L  ,L  ,L,L  ,
            U,U,D  ,D  ,D|L,L  ,L  ,L,L  ,
            U,U,D|U,D|U,D  ,D  ,D|L,L,L  ,
            U,U,D|U,D|U,D|U,D|U,D  ,D,D|L
        },
        5u,
        9u
    };

    debug_stream << "database:\t" << database << std::endl;
    debug_stream << "query:\t\t" << query << std::endl;
    debug_stream << std::endl;

    debug_stream << "trace_matrix: " << trace_matrix.cols() << " columns and "
              << trace_matrix.rows() << " rows" << std::endl;

    // Prints out the matrix in a convenient way
    debug_stream << trace_matrix << std::endl; // without sequences
    debug_stream << debug_matrix{trace_matrix, database, query} << std::endl; // with sequences
    debug_stream << fmtflags2::utf8 << debug_matrix{trace_matrix, database, query}; // as utf8
    return 0;
}
