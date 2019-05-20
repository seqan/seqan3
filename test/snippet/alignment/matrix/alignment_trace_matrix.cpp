#include <iostream>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alignment/matrix/alignment_matrix_formatter.hpp>
#include <seqan3/alignment/matrix/alignment_score_matrix.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/to_char.hpp>

int main()
{

//! [code]
using namespace seqan3;
using namespace seqan3::detail;

struct no_config
{};

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
alignment_trace_matrix trace_matrix{database, query, no_config{}, std::move(score_matrix)};

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
        if (dir == trace_directions::none)
            debug_stream << "N";
        if ((dir & trace_directions::diagonal) == trace_directions::diagonal)
            debug_stream << "D";
        if ((dir & trace_directions::up)       == trace_directions::up)
            debug_stream << "U";
        if ((dir & trace_directions::left)     == trace_directions::left)
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
