#include <iostream>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alignment/matrix/alignment_matrix_formatter.hpp>
#include <seqan3/alignment/matrix/alignment_score_matrix.hpp>
#include <seqan3/range/view/to_char.hpp>

int main()
{

//! [code]
using namespace seqan3;
using namespace seqan3::detail;
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
    5u,
    9u
};

std::cout << "database:\t" << (database | view::to_char) << std::endl;
std::cout << "query:\t\t" << (query | view::to_char) << std::endl;
std::cout << std::endl;

std::cout << "score_matrix: " << score_matrix.cols() << " columns and "
          << score_matrix.rows() << " rows" << std::endl;

// Print out the matrix.
for(unsigned row = 0u; row < score_matrix.rows(); ++row)
{
    for(unsigned col = 0u; col < score_matrix.cols(); ++col)
        std::cout << score_matrix.at(row, col) << ", ";
    std::cout << std::endl;
}
std::cout << std::endl;

// Prints out the matrix in a convenient way
alignment_matrix_formatter{score_matrix}.format(database, query);
//! [code]

    return 0;
}
