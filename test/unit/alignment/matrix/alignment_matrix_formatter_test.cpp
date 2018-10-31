#include <gtest/gtest.h>

#include <seqan3/alignment/matrix/alignment_score_matrix.hpp>
#include <seqan3/alignment/matrix/alignment_trace_matrix.hpp>
#include <seqan3/alignment/matrix/alignment_matrix_formatter.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

using namespace seqan3;
using namespace seqan3::literal;

struct no_config
{};

namespace seqan3::detail
{
struct matrix_formatter_test : public ::testing::Test
{
    std::vector<dna4> database = "AACACGTTAACCGGTT"_dna4;
    std::vector<dna4> query = "ACGTACGT"_dna4;
    std::vector<int> scores
    {
        0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16,
        1,  0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15,
        2,  1,  1,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
        3,  2,  2,  2,  2,  3,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13,
        4,  3,  3,  3,  3,  3,  4,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,
        5,  4,  3,  4,  3,  4,  4,  4,  4,  4,  5,  6,  7,  8,  9, 10, 11,
        6,  5,  4,  3,  4,  3,  4,  5,  5,  5,  5,  5,  6,  7,  8,  9, 10,
        7,  6,  5,  4,  4,  4,  3,  4,  5,  6,  6,  6,  6,  6,  7,  8,  9,
        8,  7,  6,  5,  5,  5,  4,  3,  4,  5,  6,  7,  7,  7,  7,  7,  8
    };

    detail::trace_directions N{},
        D{detail::trace_directions::diagonal},
        L{detail::trace_directions::left},
        U{detail::trace_directions::up},
        DL{D|L}, DU{D|U}, UL{U|L}, DUL{D|U|L};

    std::vector<detail::trace_directions> traces
    {
        N,  L,  L,  L,  L,  L,  L,  L,  L,  L,  L,  L,  L,  L,  L,  L,  L,
        U,  D,  DL, L,  DL, L,  L,  L,  L,  DL, DL, L,  L,  L,  L,  L,  L,
        U,  U,  D,  D,  L,  DL, L,  L,  L,  L,  L,  DL, DL, L,  L,  L,  L,
        U,  U,  DU, DU, D,  DL, D,  L,  L,  L,  L,  L,  L,  DL, DL, L,  L,
        U,  U,  DU, DU, DU, D,  DUL,D,  DL, L,  L,  L,  L,  L,  L,  DL, DL,
        U,  DU, D,  DUL,D,  DUL,D,  U,  D,  D,  DL, L,  L,  L,  L,  L,  L,
        U,  U,  U,  D,  UL, D,  L,  DUL,DU, DU, D,  D,  DL, L,  L,  L,  L,
        U,  U,  U,  U,  D,  U,  D,  L,  L,  DUL,DU, DU, D,  D,  DL, L,  L,
        U,  U,  U,  U,  DU, DU, U,  D,  DL, L,  L,  DUL,DU, DU, D,  D,  DL
    };

    auto unicode_str_length(std::string const & str)
    {
        detail::alignment_score_matrix matrix{scores, 9, 17};
        detail::alignment_matrix_formatter formatter{matrix};
        return decltype(formatter)::unicode_str_length(str);
    }
};
} // namespace seqan3::detail

using typename seqan3::detail::matrix_formatter_test;

TEST_F(matrix_formatter_test, unicode_str_length)
{
    EXPECT_EQ(unicode_str_length(" "), 1);
    EXPECT_EQ(unicode_str_length(";"), 1);
    EXPECT_EQ(unicode_str_length(""), 0);
    EXPECT_EQ(unicode_str_length("N"), 1);
    EXPECT_EQ(unicode_str_length("D"), 1);
    EXPECT_EQ(unicode_str_length("U"), 1);
    EXPECT_EQ(unicode_str_length("DU"), 2);
    EXPECT_EQ(unicode_str_length("L"), 1);
    EXPECT_EQ(unicode_str_length("DL"), 2);
    EXPECT_EQ(unicode_str_length("UL"), 2);
    EXPECT_EQ(unicode_str_length("DUL"), 3);
    EXPECT_EQ(unicode_str_length("|"), 1);
    EXPECT_EQ(unicode_str_length("-"), 1);
    EXPECT_EQ(unicode_str_length("/"), 1);
    EXPECT_EQ(unicode_str_length("INF"), 3);

    EXPECT_EQ(unicode_str_length(u8"ε"), 1);
    EXPECT_EQ(unicode_str_length(u8"║"), 1);
    EXPECT_EQ(unicode_str_length(u8"═"), 1);
    EXPECT_EQ(unicode_str_length(u8"╬"), 1);
    EXPECT_EQ(unicode_str_length(u8"∞"), 1);

    EXPECT_EQ(unicode_str_length(u8"█"), 1);
    EXPECT_EQ(unicode_str_length(u8"▘"), 1);
    EXPECT_EQ(unicode_str_length(u8"▝"), 1);
    EXPECT_EQ(unicode_str_length(u8"▀"), 1);
    EXPECT_EQ(unicode_str_length(u8"▖"), 1);
    EXPECT_EQ(unicode_str_length(u8"▌"), 1);
    EXPECT_EQ(unicode_str_length(u8"▞"), 1);
    EXPECT_EQ(unicode_str_length(u8"▛"), 1);
    EXPECT_EQ(unicode_str_length(u8"∞"), 1);

    EXPECT_EQ(unicode_str_length(u8"⠀"), 1);
    EXPECT_EQ(unicode_str_length(u8"⠁"), 1);
    EXPECT_EQ(unicode_str_length(u8"⠈"), 1);
    EXPECT_EQ(unicode_str_length(u8"⠉"), 1);
    EXPECT_EQ(unicode_str_length(u8"⠄"), 1);
    EXPECT_EQ(unicode_str_length(u8"⠅"), 1);
    EXPECT_EQ(unicode_str_length(u8"⠌"), 1);
    EXPECT_EQ(unicode_str_length(u8"⠍"), 1);

    EXPECT_EQ(unicode_str_length(u8"↺"), 1);
    EXPECT_EQ(unicode_str_length(u8"↖"), 1);
    EXPECT_EQ(unicode_str_length(u8"↑"), 1);
    EXPECT_EQ(unicode_str_length(u8"↖↑"), 2);
    EXPECT_EQ(unicode_str_length(u8"←"), 1);
    EXPECT_EQ(unicode_str_length(u8"↖←"), 2);
    EXPECT_EQ(unicode_str_length(u8"↑←"), 2);
    EXPECT_EQ(unicode_str_length(u8"↖↑←"), 3);
}

TEST_F(matrix_formatter_test, score_matrix_ascii)
{
    detail::alignment_score_matrix matrix{scores, 9, 17};
    detail::alignment_matrix_formatter formatter{matrix, detail::alignment_matrix_format::ascii};

    EXPECT_FALSE(formatter.is_traceback_matrix);
    EXPECT_EQ(formatter.auto_width(), 2u);

    std::stringstream stream;
    formatter.format(database, query, stream, 3u);
    EXPECT_EQ(stream.str(),
        " |   |A  |A  |C  |A  |C  |G  |T  |T  |A  |A  |C  |C  |G  |G  |T  |T  |\n"
        " /---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/\n"
        " |0  |1  |2  |3  |4  |5  |6  |7  |8  |9  |10 |11 |12 |13 |14 |15 |16 |\n"
        " /---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/\n"
        "A|1  |0  |1  |2  |3  |4  |5  |6  |7  |8  |9  |10 |11 |12 |13 |14 |15 |\n"
        " /---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/\n"
        "C|2  |1  |1  |1  |2  |3  |4  |5  |6  |7  |8  |9  |10 |11 |12 |13 |14 |\n"
        " /---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/\n"
        "G|3  |2  |2  |2  |2  |3  |3  |4  |5  |6  |7  |8  |9  |10 |11 |12 |13 |\n"
        " /---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/\n"
        "T|4  |3  |3  |3  |3  |3  |4  |3  |4  |5  |6  |7  |8  |9  |10 |11 |12 |\n"
        " /---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/\n"
        "A|5  |4  |3  |4  |3  |4  |4  |4  |4  |4  |5  |6  |7  |8  |9  |10 |11 |\n"
        " /---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/\n"
        "C|6  |5  |4  |3  |4  |3  |4  |5  |5  |5  |5  |5  |6  |7  |8  |9  |10 |\n"
        " /---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/\n"
        "G|7  |6  |5  |4  |4  |4  |3  |4  |5  |6  |6  |6  |6  |6  |7  |8  |9  |\n"
        " /---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/---/\n"
        "T|8  |7  |6  |5  |5  |5  |4  |3  |4  |5  |6  |7  |7  |7  |7  |7  |8  |\n"
    );
}

TEST_F(matrix_formatter_test, score_matrix_unicode)
{
    detail::alignment_score_matrix matrix{scores, 9, 17};
    detail::alignment_matrix_formatter formatter{matrix};

    EXPECT_FALSE(formatter.is_traceback_matrix);
    EXPECT_EQ(formatter.auto_width(), 2u);

    std::stringstream stream;
    formatter.format(database, query, stream, 4u);
    EXPECT_EQ(stream.str(),
        " ║ε   ║A   ║A   ║C   ║A   ║C   ║G   ║T   ║T   ║A   ║A   ║C   ║C   ║G   ║G   ║T   ║T   ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "ε║0   ║1   ║2   ║3   ║4   ║5   ║6   ║7   ║8   ║9   ║10  ║11  ║12  ║13  ║14  ║15  ║16  ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "A║1   ║0   ║1   ║2   ║3   ║4   ║5   ║6   ║7   ║8   ║9   ║10  ║11  ║12  ║13  ║14  ║15  ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "C║2   ║1   ║1   ║1   ║2   ║3   ║4   ║5   ║6   ║7   ║8   ║9   ║10  ║11  ║12  ║13  ║14  ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "G║3   ║2   ║2   ║2   ║2   ║3   ║3   ║4   ║5   ║6   ║7   ║8   ║9   ║10  ║11  ║12  ║13  ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "T║4   ║3   ║3   ║3   ║3   ║3   ║4   ║3   ║4   ║5   ║6   ║7   ║8   ║9   ║10  ║11  ║12  ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "A║5   ║4   ║3   ║4   ║3   ║4   ║4   ║4   ║4   ║4   ║5   ║6   ║7   ║8   ║9   ║10  ║11  ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "C║6   ║5   ║4   ║3   ║4   ║3   ║4   ║5   ║5   ║5   ║5   ║5   ║6   ║7   ║8   ║9   ║10  ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "G║7   ║6   ║5   ║4   ║4   ║4   ║3   ║4   ║5   ║6   ║6   ║6   ║6   ║6   ║7   ║8   ║9   ║\n"
        " ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        "T║8   ║7   ║6   ║5   ║5   ║5   ║4   ║3   ║4   ║5   ║6   ║7   ║7   ║7   ║7   ║7   ║8   ║\n"
    );
}

TEST_F(matrix_formatter_test, trace_matrix_csv)
{
    detail::alignment_trace_matrix matrix{traces, 9, 17};
    detail::alignment_matrix_formatter formatter{matrix, detail::alignment_matrix_format::csv};

    EXPECT_TRUE(formatter.is_traceback_matrix);
    EXPECT_EQ(formatter.auto_width(), 3u);

    std::stringstream stream;
    formatter.format(database, query, stream, 4u);
    EXPECT_EQ(stream.str(),
        " ;    ;A   ;A   ;C   ;A   ;C   ;G   ;T   ;T   ;A   ;A   ;C   ;C   ;G   ;G   ;T   ;T   ;\n"
        " ;N   ;L   ;L   ;L   ;L   ;L   ;L   ;L   ;L   ;L   ;L   ;L   ;L   ;L   ;L   ;L   ;L   ;\n"
        "A;U   ;D   ;DL  ;L   ;DL  ;L   ;L   ;L   ;L   ;DL  ;DL  ;L   ;L   ;L   ;L   ;L   ;L   ;\n"
        "C;U   ;U   ;D   ;D   ;L   ;DL  ;L   ;L   ;L   ;L   ;L   ;DL  ;DL  ;L   ;L   ;L   ;L   ;\n"
        "G;U   ;U   ;DU  ;DU  ;D   ;DL  ;D   ;L   ;L   ;L   ;L   ;L   ;L   ;DL  ;DL  ;L   ;L   ;\n"
        "T;U   ;U   ;DU  ;DU  ;DU  ;D   ;DUL ;D   ;DL  ;L   ;L   ;L   ;L   ;L   ;L   ;DL  ;DL  ;\n"
        "A;U   ;DU  ;D   ;DUL ;D   ;DUL ;D   ;U   ;D   ;D   ;DL  ;L   ;L   ;L   ;L   ;L   ;L   ;\n"
        "C;U   ;U   ;U   ;D   ;UL  ;D   ;L   ;DUL ;DU  ;DU  ;D   ;D   ;DL  ;L   ;L   ;L   ;L   ;\n"
        "G;U   ;U   ;U   ;U   ;D   ;U   ;D   ;L   ;L   ;DUL ;DU  ;DU  ;D   ;D   ;DL  ;L   ;L   ;\n"
        "T;U   ;U   ;U   ;U   ;DU  ;DU  ;U   ;D   ;DL  ;L   ;L   ;DUL ;DU  ;DU  ;D   ;D   ;DL  ;\n"
    );
}

TEST_F(matrix_formatter_test, trace_matrix_unicode)
{
    detail::alignment_trace_matrix matrix{traces, 9, 17};
    detail::alignment_matrix_formatter formatter{matrix, detail::alignment_matrix_format::unicode_arrows};

    EXPECT_TRUE(formatter.is_traceback_matrix);
    EXPECT_EQ(formatter.auto_width(), 3u);

    std::stringstream stream;
    formatter.format(database, query, stream, 4u);
    EXPECT_EQ(stream.str(),
        u8" ║ε   ║A   ║A   ║C   ║A   ║C   ║G   ║T   ║T   ║A   ║A   ║C   ║C   ║G   ║G   ║T   ║T   ║\n"
        u8" ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        u8"ε║↺   ║←   ║←   ║←   ║←   ║←   ║←   ║←   ║←   ║←   ║←   ║←   ║←   ║←   ║←   ║←   ║←   ║\n"
        u8" ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        u8"A║↑   ║↖   ║↖←  ║←   ║↖←  ║←   ║←   ║←   ║←   ║↖←  ║↖←  ║←   ║←   ║←   ║←   ║←   ║←   ║\n"
        u8" ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        u8"C║↑   ║↑   ║↖   ║↖   ║←   ║↖←  ║←   ║←   ║←   ║←   ║←   ║↖←  ║↖←  ║←   ║←   ║←   ║←   ║\n"
        u8" ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        u8"G║↑   ║↑   ║↖↑  ║↖↑  ║↖   ║↖←  ║↖   ║←   ║←   ║←   ║←   ║←   ║←   ║↖←  ║↖←  ║←   ║←   ║\n"
        u8" ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        u8"T║↑   ║↑   ║↖↑  ║↖↑  ║↖↑  ║↖   ║↖↑← ║↖   ║↖←  ║←   ║←   ║←   ║←   ║←   ║←   ║↖←  ║↖←  ║\n"
        u8" ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        u8"A║↑   ║↖↑  ║↖   ║↖↑← ║↖   ║↖↑← ║↖   ║↑   ║↖   ║↖   ║↖←  ║←   ║←   ║←   ║←   ║←   ║←   ║\n"
        u8" ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        u8"C║↑   ║↑   ║↑   ║↖   ║↑←  ║↖   ║←   ║↖↑← ║↖↑  ║↖↑  ║↖   ║↖   ║↖←  ║←   ║←   ║←   ║←   ║\n"
        u8" ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        u8"G║↑   ║↑   ║↑   ║↑   ║↖   ║↑   ║↖   ║←   ║←   ║↖↑← ║↖↑  ║↖↑  ║↖   ║↖   ║↖←  ║←   ║←   ║\n"
        u8" ╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬════╬\n"
        u8"T║↑   ║↑   ║↑   ║↑   ║↖↑  ║↖↑  ║↑   ║↖   ║↖←  ║←   ║←   ║↖↑← ║↖↑  ║↖↑  ║↖   ║↖   ║↖←  ║\n"
    );
}

TEST_F(matrix_formatter_test, trace_matrix_from_score_matrix_unicode)
{
    detail::alignment_trace_matrix matrix{traces, 9, 17};
    detail::alignment_matrix_formatter formatter
    {
        matrix,
        {
            u8"ε", u8"|", u8"═", u8"/",
            {u8"█",u8"▘",u8"↑",u8"⠉",u8"▖",u8"⠅",u8"▞",u8"▛"}
        }
    };

    EXPECT_TRUE(formatter.is_traceback_matrix);
    EXPECT_EQ(formatter.auto_width(), 1u);

    std::stringstream stream;
    formatter.format(database, query, stream, 1u);
    EXPECT_EQ(stream.str(),
        u8" |ε|A|A|C|A|C|G|T|T|A|A|C|C|G|G|T|T|\n"
        u8" /═/═/═/═/═/═/═/═/═/═/═/═/═/═/═/═/═/\n"
        u8"ε|█|▖|▖|▖|▖|▖|▖|▖|▖|▖|▖|▖|▖|▖|▖|▖|▖|\n"
        u8" /═/═/═/═/═/═/═/═/═/═/═/═/═/═/═/═/═/\n"
        u8"A|↑|▘|⠅|▖|⠅|▖|▖|▖|▖|⠅|⠅|▖|▖|▖|▖|▖|▖|\n"
        u8" /═/═/═/═/═/═/═/═/═/═/═/═/═/═/═/═/═/\n"
        u8"C|↑|↑|▘|▘|▖|⠅|▖|▖|▖|▖|▖|⠅|⠅|▖|▖|▖|▖|\n"
        u8" /═/═/═/═/═/═/═/═/═/═/═/═/═/═/═/═/═/\n"
        u8"G|↑|↑|⠉|⠉|▘|⠅|▘|▖|▖|▖|▖|▖|▖|⠅|⠅|▖|▖|\n"
        u8" /═/═/═/═/═/═/═/═/═/═/═/═/═/═/═/═/═/\n"
        u8"T|↑|↑|⠉|⠉|⠉|▘|▛|▘|⠅|▖|▖|▖|▖|▖|▖|⠅|⠅|\n"
        u8" /═/═/═/═/═/═/═/═/═/═/═/═/═/═/═/═/═/\n"
        u8"A|↑|⠉|▘|▛|▘|▛|▘|↑|▘|▘|⠅|▖|▖|▖|▖|▖|▖|\n"
        u8" /═/═/═/═/═/═/═/═/═/═/═/═/═/═/═/═/═/\n"
        u8"C|↑|↑|↑|▘|▞|▘|▖|▛|⠉|⠉|▘|▘|⠅|▖|▖|▖|▖|\n"
        u8" /═/═/═/═/═/═/═/═/═/═/═/═/═/═/═/═/═/\n"
        u8"G|↑|↑|↑|↑|▘|↑|▘|▖|▖|▛|⠉|⠉|▘|▘|⠅|▖|▖|\n"
        u8" /═/═/═/═/═/═/═/═/═/═/═/═/═/═/═/═/═/\n"
        u8"T|↑|↑|↑|↑|⠉|⠉|↑|▘|⠅|▖|▖|▛|⠉|⠉|▘|▘|⠅|\n"
    );
}
