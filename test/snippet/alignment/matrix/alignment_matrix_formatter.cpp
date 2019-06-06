#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alignment/matrix/alignment_matrix_formatter.hpp>
#include <seqan3/alignment/matrix/alignment_trace_matrix.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
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

    alignment_trace_matrix trace_matrix{database, query, no_config{}, score_matrix};

    debug_stream << "//! [alignment_matrix_format::out]" << std::endl;
//! [alignment_matrix_format]
using namespace seqan3;

alignment_matrix_format format
{
    "e", // epsilon
    "|", // col_sep
    "-", // row_sep
    "/", // row_col_sep
    "INF", // infinity
    {" ","D","U","DU","L","DL","UL","DUL"} // trace_dir
};

alignment_matrix_formatter{trace_matrix, format}.format(database, query);
//! [alignment_matrix_format]
    debug_stream << "//! [alignment_matrix_format::out]" << std::endl;

    debug_stream << "//! [alignment_matrix_format::ascii#score]" << std::endl;
    alignment_matrix_formatter{score_matrix, alignment_matrix_format::ascii}.format(database, query);
    debug_stream << "//! [alignment_matrix_format::ascii#score]" << std::endl;

    debug_stream << "//! [alignment_matrix_format::ascii#trace]" << std::endl;
    alignment_matrix_formatter{trace_matrix, alignment_matrix_format::ascii}.format(database, query);
    debug_stream << "//! [alignment_matrix_format::ascii#trace]" << std::endl;

    debug_stream << "//! [alignment_matrix_format::csv#score]" << std::endl;
    alignment_matrix_formatter{score_matrix, alignment_matrix_format::csv}.format(database, query);
    debug_stream << "//! [alignment_matrix_format::csv#score]" << std::endl;

    debug_stream << "//! [alignment_matrix_format::csv#trace]" << std::endl;
    alignment_matrix_formatter{trace_matrix, alignment_matrix_format::csv}.format(database, query);
    debug_stream << "//! [alignment_matrix_format::csv#trace]" << std::endl;

    debug_stream << "//! [alignment_matrix_format::unicode_block#score]" << std::endl;
    alignment_matrix_formatter{score_matrix, alignment_matrix_format::unicode_block}.format(database, query);
    debug_stream << "//! [alignment_matrix_format::unicode_block#score]" << std::endl;

    debug_stream << "//! [alignment_matrix_format::unicode_block#trace]" << std::endl;
    alignment_matrix_formatter{trace_matrix, alignment_matrix_format::unicode_block}.format(database, query);
    debug_stream << "//! [alignment_matrix_format::unicode_block#trace]" << std::endl;

    debug_stream << "//! [alignment_matrix_format::unicode_braille#score]" << std::endl;
    alignment_matrix_formatter{score_matrix, alignment_matrix_format::unicode_braille}.format(database, query);
    debug_stream << "//! [alignment_matrix_format::unicode_braille#score]" << std::endl;

    debug_stream << "//! [alignment_matrix_format::unicode_braille#trace]" << std::endl;
    alignment_matrix_formatter{trace_matrix, alignment_matrix_format::unicode_braille}.format(database, query);
    debug_stream << "//! [alignment_matrix_format::unicode_braille#trace]" << std::endl;

    debug_stream << "//! [alignment_matrix_format::unicode_arrows#score]" << std::endl;
    alignment_matrix_formatter{score_matrix, alignment_matrix_format::unicode_arrows}.format(database, query);
    debug_stream << "//! [alignment_matrix_format::unicode_arrows#score]" << std::endl;

    debug_stream << "//! [alignment_matrix_format::unicode_arrows#trace]" << std::endl;
    alignment_matrix_formatter{trace_matrix, alignment_matrix_format::unicode_arrows}.format(database, query);
    debug_stream << "//! [alignment_matrix_format::unicode_arrows#trace]" << std::endl;

    return 0;
}
