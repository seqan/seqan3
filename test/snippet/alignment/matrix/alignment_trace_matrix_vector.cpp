#include <iostream>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alignment/matrix/alignment_matrix_formatter.hpp>
#include <seqan3/alignment/matrix/alignment_trace_matrix.hpp>
#include <seqan3/range/view/to_char.hpp>

int main()
{

//! [code]
using namespace seqan3;
using namespace seqan3::detail;
using namespace seqan3::literal;

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

std::cout << "database:\t" << (database | view::to_char) << std::endl;
std::cout << "query:\t\t" << (query | view::to_char) << std::endl;
std::cout << std::endl;

std::cout << "trace_matrix: " << trace_matrix.cols() << " columns and "
          << trace_matrix.rows() << " rows" << std::endl;

// Print out the matrix.
for(unsigned row = 0u; row < trace_matrix.rows(); ++row)
{
    for(unsigned col = 0u; col < trace_matrix.cols(); ++col)
    {
        trace_directions dir = trace_matrix.at(row, col);
        if (dir == N)
            std::cout << "N";
        if ((dir & D) == D)
            std::cout << "D";
        if ((dir & U) == U)
            std::cout << "U";
        if ((dir & L) == L)
            std::cout << "L";
        std::cout << ", ";
    }
    std::cout << std::endl;
}
std::cout << std::endl;

// Prints out the matrix in a convenient way
alignment_matrix_formatter{trace_matrix}.format(database, query);
//! [code]

    return 0;
}
