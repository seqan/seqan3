// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include "../alphabet_constexpr_test_template.hpp"
#include "../alphabet_test_template.hpp"
#include "../semi_alphabet_constexpr_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"
#include "aminoacid_test_template.hpp"

#include <seqan3/alphabet/aminoacid/aa20.hpp>
#include <seqan3/range/views/zip.hpp>

using seqan3::operator""_aa20;

INSTANTIATE_TYPED_TEST_SUITE_P(aa20, alphabet_, seqan3::aa20, );
INSTANTIATE_TYPED_TEST_SUITE_P(aa20, semi_alphabet_test, seqan3::aa20, );
INSTANTIATE_TYPED_TEST_SUITE_P(aa20, alphabet_constexpr, seqan3::aa20, );
INSTANTIATE_TYPED_TEST_SUITE_P(aa20, semi_alphabet_constexpr, seqan3::aa20, );
INSTANTIATE_TYPED_TEST_SUITE_P(aa20, aminoacid, seqan3::aa20, );

TEST(aa20, assign_char)
{
    std::vector<char> chars
    {
        'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
        'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm',
        'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
        'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z',
        '*', '!'
    };

    std::vector<seqan3::aa20> alphabets
    {
        'A'_aa20, 'D'_aa20, 'C'_aa20, 'D'_aa20, 'E'_aa20, 'F'_aa20, 'G'_aa20,
        'H'_aa20, 'I'_aa20, 'L'_aa20, 'K'_aa20, 'L'_aa20, 'M'_aa20,
        'A'_aa20, 'D'_aa20, 'C'_aa20, 'D'_aa20, 'E'_aa20, 'F'_aa20, 'G'_aa20,
        'H'_aa20, 'I'_aa20, 'L'_aa20, 'K'_aa20, 'L'_aa20, 'M'_aa20,
        'N'_aa20, 'L'_aa20, 'P'_aa20, 'Q'_aa20, 'R'_aa20, 'S'_aa20, 'T'_aa20,
        'C'_aa20, 'V'_aa20, 'W'_aa20, 'S'_aa20, 'Y'_aa20, 'E'_aa20,
        'N'_aa20, 'L'_aa20, 'P'_aa20, 'Q'_aa20, 'R'_aa20, 'S'_aa20, 'T'_aa20,
        'C'_aa20, 'V'_aa20, 'W'_aa20, 'S'_aa20, 'Y'_aa20, 'E'_aa20,
        'W'_aa20, 'S'_aa20
    };

    for (auto [ chr, alp ] : seqan3::views::zip(chars, alphabets))
        EXPECT_EQ((seqan3::assign_char_to(chr, seqan3::aa20{})), alp);
}

TEST(aa20, to_char)
{
    std::vector<char> chars
    {
        'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P',
        'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'D', 'L', 'L', 'C', 'S', 'E',
        'W', 'S'
    };

    std::vector<seqan3::aa20> alphabets
    {
        'A'_aa20, 'C'_aa20, 'D'_aa20, 'E'_aa20, 'F'_aa20, 'G'_aa20, 'H'_aa20,
        'I'_aa20, 'K'_aa20, 'L'_aa20, 'M'_aa20, 'N'_aa20, 'P'_aa20,
        'Q'_aa20, 'R'_aa20, 'S'_aa20, 'T'_aa20, 'V'_aa20, 'W'_aa20, 'Y'_aa20,
        'B'_aa20, 'J'_aa20, 'O'_aa20, 'U'_aa20, 'X'_aa20, 'Z'_aa20,
        'W'_aa20, 'S'_aa20
    };

    for (auto [ chr, alp ] : seqan3::views::zip(chars, alphabets))
        EXPECT_EQ(seqan3::to_char(alp), chr);
}

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

TEST(literals, char_literal)
{
    EXPECT_EQ(seqan3::to_char('A'_aa20), 'A');
    EXPECT_EQ(seqan3::to_char('C'_aa20), 'C');
    EXPECT_EQ(seqan3::to_char('D'_aa20), 'D');
    EXPECT_EQ(seqan3::to_char('E'_aa20), 'E');
    EXPECT_EQ(seqan3::to_char('F'_aa20), 'F');
    EXPECT_EQ(seqan3::to_char('G'_aa20), 'G');
    EXPECT_EQ(seqan3::to_char('H'_aa20), 'H');
    EXPECT_EQ(seqan3::to_char('I'_aa20), 'I');
    EXPECT_EQ(seqan3::to_char('K'_aa20), 'K');
    EXPECT_EQ(seqan3::to_char('L'_aa20), 'L');
    EXPECT_EQ(seqan3::to_char('M'_aa20), 'M');
    EXPECT_EQ(seqan3::to_char('N'_aa20), 'N');
    EXPECT_EQ(seqan3::to_char('P'_aa20), 'P');
    EXPECT_EQ(seqan3::to_char('Q'_aa20), 'Q');
    EXPECT_EQ(seqan3::to_char('R'_aa20), 'R');
    EXPECT_EQ(seqan3::to_char('S'_aa20), 'S');
    EXPECT_EQ(seqan3::to_char('T'_aa20), 'T');
    EXPECT_EQ(seqan3::to_char('V'_aa20), 'V');
    EXPECT_EQ(seqan3::to_char('W'_aa20), 'W');
    EXPECT_EQ(seqan3::to_char('Y'_aa20), 'Y');

    EXPECT_EQ(seqan3::to_char('*'_aa20), 'W');
    EXPECT_EQ(seqan3::to_char('!'_aa20), 'S');
}

TEST(literals, vector)
{
    seqan3::aa20_vector v20;
    v20.resize(5, 'B'_aa20);
    EXPECT_EQ(v20, "DDDDD"_aa20);

    std::vector<seqan3::aa20> w20{'A'_aa20, 'B'_aa20, 'J'_aa20, 'O'_aa20, 'U'_aa20, 'X'_aa20, 'Z'_aa20, '!'_aa20,
                                  '*'_aa20, '*'_aa20};
    EXPECT_EQ(w20, "ADLLCSESW*"_aa20);
}

TEST(aa20, char_is_valid)
{
    constexpr auto aa27_validator = (seqan3::is_alpha || seqan3::is_char<'*'>);

    for (char c : std::views::iota(std::numeric_limits<char>::min(), std::numeric_limits<char>::max()))
    {
        bool expect = false;
        switch (c)
        {
            case 'B': case 'J': case 'O': case 'U': case 'X': case 'Z':
            case 'b': case 'j': case 'o': case 'u': case 'x': case 'z':
            case '*':
                break;
            default:
                expect = aa27_validator(c);
                break;
        }

        EXPECT_EQ(seqan3::aa20::char_is_valid(c), expect);
    }
}
