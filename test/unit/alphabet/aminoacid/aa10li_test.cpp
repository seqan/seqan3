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

#include <seqan3/alphabet/aminoacid/aa10li.hpp>
#include <seqan3/range/views/zip.hpp>

using seqan3::operator""_aa10li;

INSTANTIATE_TYPED_TEST_SUITE_P(aa10li, alphabet_, seqan3::aa10li, );
INSTANTIATE_TYPED_TEST_SUITE_P(aa10li, semi_alphabet_test, seqan3::aa10li, );
INSTANTIATE_TYPED_TEST_SUITE_P(aa10li, alphabet_constexpr, seqan3::aa10li, );
INSTANTIATE_TYPED_TEST_SUITE_P(aa10li, semi_alphabet_constexpr, seqan3::aa10li, );
INSTANTIATE_TYPED_TEST_SUITE_P(aa10li, aminoacid, seqan3::aa10li, );

TEST(aa10li, assign_char)
{
    std::vector<char> chars
    {
        'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
        'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm',
        'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
        'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z',
        '*', '!'
    };

    std::vector<seqan3::aa10li> alphabets
    {
        'A'_aa10li, 'B'_aa10li, 'C'_aa10li, 'B'_aa10li, 'B'_aa10li, 'F'_aa10li, 'G'_aa10li,
        'H'_aa10li, 'I'_aa10li, 'J'_aa10li, 'K'_aa10li, 'J'_aa10li, 'J'_aa10li,
        'A'_aa10li, 'B'_aa10li, 'C'_aa10li, 'B'_aa10li, 'B'_aa10li, 'F'_aa10li, 'G'_aa10li,
        'H'_aa10li, 'I'_aa10li, 'J'_aa10li, 'K'_aa10li, 'J'_aa10li, 'J'_aa10li,
        'H'_aa10li, 'K'_aa10li, 'P'_aa10li, 'B'_aa10li, 'K'_aa10li, 'A'_aa10li, 'A'_aa10li,
        'C'_aa10li, 'I'_aa10li, 'F'_aa10li, 'A'_aa10li, 'F'_aa10li, 'B'_aa10li,
        'H'_aa10li, 'K'_aa10li, 'P'_aa10li, 'B'_aa10li, 'K'_aa10li, 'A'_aa10li, 'A'_aa10li,
        'C'_aa10li, 'I'_aa10li, 'F'_aa10li, 'A'_aa10li, 'F'_aa10li, 'B'_aa10li,
        'F'_aa10li, 'A'_aa10li
    };

    for (auto [ chr, alp ] : seqan3::views::zip(chars, alphabets))
        EXPECT_EQ((seqan3::assign_char_to(chr, seqan3::aa10li{})), alp);
}

TEST(aa10li, to_char)
{
    std::vector<char> chars
    {
        'A', 'C', 'B', 'B', 'F', 'G', 'H', 'I', 'K', 'J', 'J', 'H', 'P',
        'B', 'K', 'A', 'A', 'I', 'F', 'F', 'B', 'J', 'K', 'C', 'A', 'B'
    };

    std::vector<seqan3::aa10li> alphabets
    {
        'A'_aa10li, 'C'_aa10li, 'D'_aa10li, 'E'_aa10li, 'F'_aa10li, 'G'_aa10li, 'H'_aa10li,
        'I'_aa10li, 'K'_aa10li, 'L'_aa10li, 'M'_aa10li, 'N'_aa10li, 'P'_aa10li,
        'Q'_aa10li, 'R'_aa10li, 'S'_aa10li, 'T'_aa10li, 'V'_aa10li, 'W'_aa10li, 'Y'_aa10li,
        'B'_aa10li, 'J'_aa10li, 'O'_aa10li, 'U'_aa10li, 'X'_aa10li, 'Z'_aa10li
    };

    for (auto [ chr, alp ] : seqan3::views::zip(chars, alphabets))
        EXPECT_EQ(seqan3::to_char(alp), chr);
}

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

TEST(literals, char_literal)
{
    EXPECT_EQ(seqan3::to_char('A'_aa10li), 'A');
    EXPECT_EQ(seqan3::to_char('B'_aa10li), 'B');
    EXPECT_EQ(seqan3::to_char('C'_aa10li), 'C');
    EXPECT_EQ(seqan3::to_char('F'_aa10li), 'F');
    EXPECT_EQ(seqan3::to_char('G'_aa10li), 'G');
    EXPECT_EQ(seqan3::to_char('H'_aa10li), 'H');
    EXPECT_EQ(seqan3::to_char('I'_aa10li), 'I');
    EXPECT_EQ(seqan3::to_char('J'_aa10li), 'J');
    EXPECT_EQ(seqan3::to_char('K'_aa10li), 'K');
    EXPECT_EQ(seqan3::to_char('P'_aa10li), 'P');

    EXPECT_EQ(seqan3::to_char('*'_aa10li), 'F');
    EXPECT_EQ(seqan3::to_char('!'_aa10li), 'A');
}

TEST(literals, vector)
{
    seqan3::aa10li_vector v20;
    v20.resize(5, 'D'_aa10li);
    EXPECT_EQ(v20, "BBBBB"_aa10li);

    std::vector<seqan3::aa10li> w20{'A'_aa10li, 'D'_aa10li, 'N'_aa10li, 'O'_aa10li, 'U'_aa10li, 'X'_aa10li, 'R'_aa10li,
                                    '!'_aa10li, '*'_aa10li, '*'_aa10li};
    EXPECT_EQ(w20, "ABHKCAKAF*"_aa10li);
}

TEST(aa10li, char_is_valid)
{
    constexpr auto aa27_validator = (seqan3::is_alpha || seqan3::is_char<'*'>);

    for (char c : std::views::iota(std::numeric_limits<char>::min(), std::numeric_limits<char>::max()))
    {
        bool expect = false;
        switch (c)
        {
            case 'D': case 'E': case 'L': case 'M': case 'N': case 'O': case 'Q':
            case 'R': case 'S': case 'T': case 'U': case 'V': case 'W': case 'X':
            case 'Y': case 'Z':
            case 'd': case 'e': case 'l': case 'm': case 'n': case 'o': case 'q':
            case 'r': case 's': case 't': case 'u': case 'v': case 'w': case 'x':
            case 'y': case 'z':
            case '*':
                break;
            default:
                expect = aa27_validator(c);
                break;
        }

        EXPECT_EQ(seqan3::aa10li::char_is_valid(c), expect);
    }
}
