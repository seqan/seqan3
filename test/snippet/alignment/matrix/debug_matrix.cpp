#include <iostream>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alignment/matrix/debug_matrix.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/to_char.hpp>

int main()
{
    using namespace seqan3;
{
    debug_stream << "//! [score_matrix::out]" << std::endl;
//! [score_matrix]
using namespace seqan3;
using namespace seqan3::detail;

std::vector<dna4> database = "AACCGGTT"_dna4;
std::vector<dna4> query = "ACGT"_dna4;

debug_matrix score_matrix
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
//! [score_matrix]
    debug_stream << "//! [score_matrix::out]" << std::endl;
    debug_stream.unsetf(fmtflags2::utf8);
}
{
    debug_stream << "//! [trace_matrix::out]" << std::endl;
//! [trace_matrix]
using namespace seqan3;
using namespace seqan3::detail;

std::vector<dna4> database = "AACCGGTT"_dna4;
std::vector<dna4> query = "ACGT"_dna4;

auto N = trace_directions::none;
auto D = trace_directions::diagonal;
auto U = trace_directions::up;
auto L = trace_directions::left;

debug_matrix trace_matrix
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
//! [trace_matrix]
    debug_stream << "//! [trace_matrix::out]" << std::endl;
}
    return 0;
}
