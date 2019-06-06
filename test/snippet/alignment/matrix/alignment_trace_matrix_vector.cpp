#include <iostream>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alignment/matrix/alignment_matrix_formatter.hpp>
#include <seqan3/alignment/matrix/alignment_trace_matrix.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/to_char.hpp>

int main()
{

//! [code]
using namespace seqan3;
using namespace seqan3::detail;

std::vector<dna4> database = "AACCGGTT"_dna4;
std::vector<dna4> query = "ACGT"_dna4;

auto N = trace_directions::none;
auto D = trace_directions::diagonal;
auto U = trace_directions::up;
auto L = trace_directions::left;

alignment_trace_matrix trace_matrix
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

debug_stream << "database:\t" << (database | view::to_char) << std::endl;
debug_stream << "query:\t\t" << (query | view::to_char) << std::endl;
debug_stream << std::endl;

debug_stream << "trace_matrix: " << trace_matrix.cols() << " columns and "
          << trace_matrix.rows() << " rows" << std::endl;

// Print out the matrix.
for(unsigned row = 0u; row < trace_matrix.rows(); ++row)
{
    for(unsigned col = 0u; col < trace_matrix.cols(); ++col)
    {
        trace_directions dir = trace_matrix.at(row, col);
        if (dir == N)
            debug_stream << "N";
        if ((dir & D) == D)
            debug_stream << "D";
        if ((dir & U) == U)
            debug_stream << "U";
        if ((dir & L) == L)
            debug_stream << "L";
        debug_stream << ", ";
    }
    debug_stream << std::endl;
}
debug_stream << std::endl;

// Prints out the matrix in a convenient way
alignment_matrix_formatter{trace_matrix}.format(database, query);
//! [code]

    return 0;
}
