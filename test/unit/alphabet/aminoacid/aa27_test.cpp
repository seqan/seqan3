// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/core/debug_stream/range.hpp>
#include <seqan3/utility/views/zip.hpp>

#include "../alphabet_constexpr_test_template.hpp"
#include "../alphabet_test_template.hpp"
#include "../semi_alphabet_constexpr_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"
#include "aminoacid_test_template.hpp"

using seqan3::operator""_aa27;

INSTANTIATE_TYPED_TEST_SUITE_P(aa27, alphabet, seqan3::aa27, );
INSTANTIATE_TYPED_TEST_SUITE_P(aa27, semi_alphabet_test, seqan3::aa27, );
INSTANTIATE_TYPED_TEST_SUITE_P(aa27, alphabet_constexpr, seqan3::aa27, );
INSTANTIATE_TYPED_TEST_SUITE_P(aa27, semi_alphabet_constexpr, seqan3::aa27, );
INSTANTIATE_TYPED_TEST_SUITE_P(aa27, aminoacid, seqan3::aa27, );

TEST(aa27, assign_char)
{
    std::vector<char> chars{'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'a', 'b', 'c', 'd', 'e',
                            'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W',
                            'X', 'Y', 'Z', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '*', '!'};

    std::vector<seqan3::aa27> alphabets{
        'A'_aa27, 'B'_aa27, 'C'_aa27, 'D'_aa27, 'E'_aa27, 'F'_aa27, 'G'_aa27, 'H'_aa27, 'I'_aa27, 'J'_aa27, 'K'_aa27,
        'L'_aa27, 'M'_aa27, 'A'_aa27, 'B'_aa27, 'C'_aa27, 'D'_aa27, 'E'_aa27, 'F'_aa27, 'G'_aa27, 'H'_aa27, 'I'_aa27,
        'J'_aa27, 'K'_aa27, 'L'_aa27, 'M'_aa27, 'N'_aa27, 'O'_aa27, 'P'_aa27, 'Q'_aa27, 'R'_aa27, 'S'_aa27, 'T'_aa27,
        'U'_aa27, 'V'_aa27, 'W'_aa27, 'X'_aa27, 'Y'_aa27, 'Z'_aa27, 'N'_aa27, 'O'_aa27, 'P'_aa27, 'Q'_aa27, 'R'_aa27,
        'S'_aa27, 'T'_aa27, 'U'_aa27, 'V'_aa27, 'W'_aa27, 'X'_aa27, 'Y'_aa27, 'Z'_aa27, '*'_aa27, 'X'_aa27};

    for (auto [chr, alp] : seqan3::views::zip(chars, alphabets))
        EXPECT_EQ((seqan3::assign_char_to(chr, seqan3::aa27{})), alp);
}

TEST(aa27, to_char)
{
    std::vector<char> chars{'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q',
                            'R', 'S', 'T', 'V', 'W', 'Y', 'B', 'J', 'O', 'U', 'X', 'Z', '*', 'X'};

    std::vector<seqan3::aa27> alphabets{'A'_aa27, 'C'_aa27, 'D'_aa27, 'E'_aa27, 'F'_aa27, 'G'_aa27, 'H'_aa27,
                                        'I'_aa27, 'K'_aa27, 'L'_aa27, 'M'_aa27, 'N'_aa27, 'P'_aa27, 'Q'_aa27,
                                        'R'_aa27, 'S'_aa27, 'T'_aa27, 'V'_aa27, 'W'_aa27, 'Y'_aa27, 'B'_aa27,
                                        'J'_aa27, 'O'_aa27, 'U'_aa27, 'X'_aa27, 'Z'_aa27, '*'_aa27, 'X'_aa27};

    for (auto [alp, chr] : seqan3::views::zip(alphabets, chars))
        EXPECT_EQ(seqan3::to_char(alp), chr);
}

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

TEST(literals, char_literal)
{
    EXPECT_EQ(seqan3::to_char('A'_aa27), 'A');
    EXPECT_EQ(seqan3::to_char('B'_aa27), 'B');
    EXPECT_EQ(seqan3::to_char('C'_aa27), 'C');
    EXPECT_EQ(seqan3::to_char('D'_aa27), 'D');
    EXPECT_EQ(seqan3::to_char('E'_aa27), 'E');
    EXPECT_EQ(seqan3::to_char('F'_aa27), 'F');
    EXPECT_EQ(seqan3::to_char('G'_aa27), 'G');
    EXPECT_EQ(seqan3::to_char('H'_aa27), 'H');
    EXPECT_EQ(seqan3::to_char('I'_aa27), 'I');
    EXPECT_EQ(seqan3::to_char('J'_aa27), 'J');
    EXPECT_EQ(seqan3::to_char('K'_aa27), 'K');
    EXPECT_EQ(seqan3::to_char('L'_aa27), 'L');
    EXPECT_EQ(seqan3::to_char('M'_aa27), 'M');
    EXPECT_EQ(seqan3::to_char('N'_aa27), 'N');
    EXPECT_EQ(seqan3::to_char('O'_aa27), 'O');
    EXPECT_EQ(seqan3::to_char('P'_aa27), 'P');
    EXPECT_EQ(seqan3::to_char('Q'_aa27), 'Q');
    EXPECT_EQ(seqan3::to_char('R'_aa27), 'R');
    EXPECT_EQ(seqan3::to_char('S'_aa27), 'S');
    EXPECT_EQ(seqan3::to_char('T'_aa27), 'T');
    EXPECT_EQ(seqan3::to_char('U'_aa27), 'U');
    EXPECT_EQ(seqan3::to_char('V'_aa27), 'V');
    EXPECT_EQ(seqan3::to_char('W'_aa27), 'W');
    EXPECT_EQ(seqan3::to_char('X'_aa27), 'X');
    EXPECT_EQ(seqan3::to_char('Y'_aa27), 'Y');
    EXPECT_EQ(seqan3::to_char('Z'_aa27), 'Z');
    EXPECT_EQ(seqan3::to_char('*'_aa27), '*');

    EXPECT_EQ(seqan3::to_char('!'_aa27), 'X');
}

TEST(literals, vector)
{
    seqan3::aa27_vector v27;
    v27.resize(5, 'A'_aa27);
    EXPECT_EQ(v27, "AAAAA"_aa27);

    std::vector<seqan3::aa27>
        w27{'A'_aa27, 'Y'_aa27, 'P'_aa27, 'T'_aa27, 'U'_aa27, 'N'_aa27, 'X'_aa27, '!'_aa27, '*'_aa27};
    EXPECT_EQ(w27, "AYPTUNXX*"_aa27);
}

TEST(aa27, char_is_valid)
{
    constexpr auto validator = seqan3::is_alpha || seqan3::is_char<'*'>;
    for (char c : std::views::iota(std::numeric_limits<char>::min(), std::numeric_limits<char>::max()))
        EXPECT_EQ(seqan3::aa27::char_is_valid(c), validator(c));
}
