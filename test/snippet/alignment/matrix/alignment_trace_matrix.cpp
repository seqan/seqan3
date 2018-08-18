#include <iostream>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alignment/matrix/alignment_matrix_formatter.hpp>
#include <seqan3/alignment/matrix/alignment_score_matrix.hpp>
#include <seqan3/range/view/to_char.hpp>

int main()
{

//! [code]
using namespace seqan3;
using namespace seqan3::literal;

std::vector<dna4> database = "AACCGGTT"_dna4;
std::vector<dna4> query = "ACGT"_dna4;

alignment_score_matrix score_matrix
{
    std::vector
    {
        0, 1, 2, 3, 4, 5, 6, 7, 8,
        1, 0, 1, 2, 3, 4, 5, 6, 7,
        2, 1, 1, 1, 2, 3, 4, 5, 6,
        3, 2, 2, 2, 2, 2, 3, 4, 5,
        4, 3, 3, 3, 3, 3, 3, 3, 4
    },
    database,
    query
};
alignment_trace_matrix trace_matrix{std::move(score_matrix)};

std::cout << "database:\t" << (trace_matrix.database() | view::to_char) << std::endl;
std::cout << "query:\t\t" << (trace_matrix.query() | view::to_char) << std::endl;
std::cout << std::endl;

std::cout << "trace_matrix: " << trace_matrix.cols() << " columns and "
          << trace_matrix.rows() << " rows" << std::endl;

// Print out the matrix.
for(unsigned row = 0u; row < trace_matrix.rows(); ++row)
{
    for(unsigned col = 0u; col < trace_matrix.cols(); ++col)
    {
        trace_matrix_directions dir = trace_matrix.at(row, col);
        if (dir == trace_matrix_directions::none)
            std::cout << "N";
        if ((dir & trace_matrix_directions::diagonal) == trace_matrix_directions::diagonal)
            std::cout << "D";
        if ((dir & trace_matrix_directions::up)       == trace_matrix_directions::up)
            std::cout << "U";
        if ((dir & trace_matrix_directions::left)     == trace_matrix_directions::left)
            std::cout << "L";
        std::cout << ", ";
    }
    std::cout << std::endl;
}
std::cout << std::endl;

// Prints out the matrix in a convenient way
alignment_matrix_formatter{trace_matrix}.format();
//! [code]

    return 0;
}
