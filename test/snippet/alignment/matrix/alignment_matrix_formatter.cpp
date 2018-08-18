#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alignment/matrix/alignment_matrix_formatter.hpp>
#include <seqan3/alignment/matrix/alignment_trace_matrix.hpp>

int main()
{
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

    alignment_trace_matrix trace_matrix{score_matrix};

    std::cout << "//! [alignment_matrix_format::out]" << std::endl;
//! [alignment_matrix_format]
using namespace seqan3;

alignment_matrix_format format
{
    "e", // epsilon
    "|", // col_sep
    "-", // row_sep
    "/", // row_col_sep
    {" ","D","U","DU","L","DL","UL","DUL"} // trace_dir
};

alignment_matrix_formatter{trace_matrix, format}.format();
//! [alignment_matrix_format]
    std::cout << "//! [alignment_matrix_format::out]" << std::endl;

    std::cout << "//! [alignment_matrix_format::ascii#score]" << std::endl;
    alignment_matrix_formatter{score_matrix, alignment_matrix_format::ascii}.format();
    std::cout << "//! [alignment_matrix_format::ascii#score]" << std::endl;

    std::cout << "//! [alignment_matrix_format::ascii#trace]" << std::endl;
    alignment_matrix_formatter{trace_matrix, alignment_matrix_format::ascii}.format();
    std::cout << "//! [alignment_matrix_format::ascii#trace]" << std::endl;

    std::cout << "//! [alignment_matrix_format::csv#score]" << std::endl;
    alignment_matrix_formatter{score_matrix, alignment_matrix_format::csv}.format();
    std::cout << "//! [alignment_matrix_format::csv#score]" << std::endl;

    std::cout << "//! [alignment_matrix_format::csv#trace]" << std::endl;
    alignment_matrix_formatter{trace_matrix, alignment_matrix_format::csv}.format();
    std::cout << "//! [alignment_matrix_format::csv#trace]" << std::endl;

    std::cout << "//! [alignment_matrix_format::unicode_block#score]" << std::endl;
    alignment_matrix_formatter{score_matrix, alignment_matrix_format::unicode_block}.format();
    std::cout << "//! [alignment_matrix_format::unicode_block#score]" << std::endl;

    std::cout << "//! [alignment_matrix_format::unicode_block#trace]" << std::endl;
    alignment_matrix_formatter{trace_matrix, alignment_matrix_format::unicode_block}.format();
    std::cout << "//! [alignment_matrix_format::unicode_block#trace]" << std::endl;

    std::cout << "//! [alignment_matrix_format::unicode_braille#score]" << std::endl;
    alignment_matrix_formatter{score_matrix, alignment_matrix_format::unicode_braille}.format();
    std::cout << "//! [alignment_matrix_format::unicode_braille#score]" << std::endl;

    std::cout << "//! [alignment_matrix_format::unicode_braille#trace]" << std::endl;
    alignment_matrix_formatter{trace_matrix, alignment_matrix_format::unicode_braille}.format();
    std::cout << "//! [alignment_matrix_format::unicode_braille#trace]" << std::endl;

    std::cout << "//! [alignment_matrix_format::unicode_arrows#score]" << std::endl;
    alignment_matrix_formatter{score_matrix, alignment_matrix_format::unicode_arrows}.format();
    std::cout << "//! [alignment_matrix_format::unicode_arrows#score]" << std::endl;

    std::cout << "//! [alignment_matrix_format::unicode_arrows#trace]" << std::endl;
    alignment_matrix_formatter{trace_matrix, alignment_matrix_format::unicode_arrows}.format();
    std::cout << "//! [alignment_matrix_format::unicode_arrows#trace]" << std::endl;

    return 0;
}
