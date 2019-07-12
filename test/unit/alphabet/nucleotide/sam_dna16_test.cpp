// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/alphabet/nucleotide/sam_dna16.hpp>

#include <seqan3/core/char_operations/predicate.hpp>

#include "../alphabet_test_template.hpp"
#include "../alphabet_constexpr_test_template.hpp"

using namespace seqan3;

// ------------------------------------------------------------------
// sam_dna16 alphabet
// ------------------------------------------------------------------

INSTANTIATE_TYPED_TEST_CASE_P(sam_dna16, alphabet, sam_dna16);
INSTANTIATE_TYPED_TEST_CASE_P(sam_dna16, alphabet_constexpr, sam_dna16);

// nucleotide test: (because the complement is not bijective for sam_dna16 we need to test it manually)
TEST(sam_dna16, nucleotide)
{
    EXPECT_TRUE(NucleotideAlphabet<sam_dna16>);
    EXPECT_TRUE(NucleotideAlphabet<sam_dna16 &>);

    EXPECT_EQ(complement('='_sam_dna16), 'N'_sam_dna16);
    EXPECT_EQ(complement('A'_sam_dna16), 'T'_sam_dna16);
    EXPECT_EQ(complement('C'_sam_dna16), 'G'_sam_dna16);
    EXPECT_EQ(complement('M'_sam_dna16), 'K'_sam_dna16);
    EXPECT_EQ(complement('G'_sam_dna16), 'C'_sam_dna16);
    EXPECT_EQ(complement('R'_sam_dna16), 'Y'_sam_dna16);
    EXPECT_EQ(complement('S'_sam_dna16), 'S'_sam_dna16);
    EXPECT_EQ(complement('V'_sam_dna16), 'B'_sam_dna16);
    EXPECT_EQ(complement('T'_sam_dna16), 'A'_sam_dna16);
    EXPECT_EQ(complement('W'_sam_dna16), 'W'_sam_dna16);
    EXPECT_EQ(complement('Y'_sam_dna16), 'R'_sam_dna16);
    EXPECT_EQ(complement('H'_sam_dna16), 'D'_sam_dna16);
    EXPECT_EQ(complement('K'_sam_dna16), 'M'_sam_dna16);
    EXPECT_EQ(complement('D'_sam_dna16), 'H'_sam_dna16);
    EXPECT_EQ(complement('B'_sam_dna16), 'V'_sam_dna16);
    EXPECT_EQ(complement('N'_sam_dna16), 'N'_sam_dna16);
}

TEST(sam_dna16, to_char_assign_char)
{
    for (char chr : sam_dna16::rank_to_char)
        EXPECT_EQ(to_char(sam_dna16{}.assign_char(chr)), chr);

    EXPECT_EQ(to_char(sam_dna16{}.assign_char('a')), 'A');
    EXPECT_EQ(to_char(sam_dna16{}.assign_char('c')), 'C');
    EXPECT_EQ(to_char(sam_dna16{}.assign_char('g')), 'G');
    EXPECT_EQ(to_char(sam_dna16{}.assign_char('t')), 'T');

    EXPECT_EQ(to_char(sam_dna16{}.assign_char('U')), 'T');
    EXPECT_EQ(to_char(sam_dna16{}.assign_char('!')), 'N');
}

TEST(sam_dna16, char_literal)
{
    EXPECT_EQ(to_char('A'_sam_dna16), 'A');
    EXPECT_EQ(to_char('C'_sam_dna16), 'C');
    EXPECT_EQ(to_char('G'_sam_dna16), 'G');

    EXPECT_EQ(to_char('U'_sam_dna16), 'T');
    EXPECT_EQ(to_char('T'_sam_dna16), 'T');

    EXPECT_EQ(to_char('R'_sam_dna16), 'R');
    EXPECT_EQ(to_char('Y'_sam_dna16), 'Y');
    EXPECT_EQ(to_char('S'_sam_dna16), 'S');
    EXPECT_EQ(to_char('W'_sam_dna16), 'W');
    EXPECT_EQ(to_char('K'_sam_dna16), 'K');
    EXPECT_EQ(to_char('M'_sam_dna16), 'M');
    EXPECT_EQ(to_char('B'_sam_dna16), 'B');
    EXPECT_EQ(to_char('D'_sam_dna16), 'D');
    EXPECT_EQ(to_char('H'_sam_dna16), 'H');
    EXPECT_EQ(to_char('V'_sam_dna16), 'V');

    EXPECT_EQ(to_char('='_sam_dna16), '=');

    EXPECT_EQ(to_char('N'_sam_dna16), 'N');
    EXPECT_EQ(to_char('!'_sam_dna16), 'N');
}

TEST(sam_dna16, string_literal)
{
    sam_dna16_vector v;
    v.resize(5, 'A'_sam_dna16);
    EXPECT_EQ(v, "AAAAA"_sam_dna16);

    std::vector<sam_dna16> w{'A'_sam_dna16, '='_sam_dna16, 'G'_sam_dna16, 'T'_sam_dna16, 'U'_sam_dna16, 'N'_sam_dna16};
    EXPECT_EQ(w, "A=GTTN"_sam_dna16);
}

TEST(sam_dna16, char_is_valid)
{
    constexpr auto validator = is_char<'A'> || is_char<'C'> || is_char<'G'> || is_char<'T'> || is_char<'U'> ||
                               is_char<'a'> || is_char<'c'> || is_char<'g'> || is_char<'t'> || is_char<'u'> ||
                               is_char<'N'> || is_char<'n'> ||
                               is_char<'R'> || is_char<'Y'> || is_char<'S'> || is_char<'W'> || is_char<'K'> ||
                               is_char<'M'> || is_char<'B'> || is_char<'D'> || is_char<'H'> || is_char<'V'> ||
                               is_char<'r'> || is_char<'y'> || is_char<'s'> || is_char<'w'> || is_char<'k'> ||
                               is_char<'m'> || is_char<'b'> || is_char<'d'> || is_char<'h'> || is_char<'v'> ||
                               is_char<'='>;

    for (char c : std::view::iota(std::numeric_limits<char>::min(), std::numeric_limits<char>::max()))
        EXPECT_EQ(sam_dna16::char_is_valid(c), validator(c));
}
