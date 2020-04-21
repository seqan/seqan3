// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/alphabet/nucleotide/sam_dna16.hpp>

#include <seqan3/core/char_operations/predicate.hpp>

#include "../alphabet_constexpr_test_template.hpp"
#include "../alphabet_test_template.hpp"
#include "../semi_alphabet_constexpr_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"

using seqan3::operator""_sam_dna16;

// ------------------------------------------------------------------
// sam_dna16 alphabet
// ------------------------------------------------------------------

INSTANTIATE_TYPED_TEST_SUITE_P(sam_dna16, alphabet_, seqan3::sam_dna16, );
INSTANTIATE_TYPED_TEST_SUITE_P(sam_dna16, semi_alphabet_test, seqan3::sam_dna16, );
INSTANTIATE_TYPED_TEST_SUITE_P(sam_dna16, alphabet_constexpr, seqan3::sam_dna16, );
INSTANTIATE_TYPED_TEST_SUITE_P(sam_dna16, semi_alphabet_constexpr, seqan3::sam_dna16, );

// nucleotide test: (because the complement is not bijective for sam_dna16 we need to test it manually)
TEST(sam_dna16, nucleotide)
{
    EXPECT_TRUE(seqan3::nucleotide_alphabet<seqan3::sam_dna16>);
    EXPECT_TRUE(seqan3::nucleotide_alphabet<seqan3::sam_dna16 &>);

    EXPECT_EQ(seqan3::complement('='_sam_dna16), 'N'_sam_dna16);
    EXPECT_EQ(seqan3::complement('A'_sam_dna16), 'T'_sam_dna16);
    EXPECT_EQ(seqan3::complement('C'_sam_dna16), 'G'_sam_dna16);
    EXPECT_EQ(seqan3::complement('M'_sam_dna16), 'K'_sam_dna16);
    EXPECT_EQ(seqan3::complement('G'_sam_dna16), 'C'_sam_dna16);
    EXPECT_EQ(seqan3::complement('R'_sam_dna16), 'Y'_sam_dna16);
    EXPECT_EQ(seqan3::complement('S'_sam_dna16), 'S'_sam_dna16);
    EXPECT_EQ(seqan3::complement('V'_sam_dna16), 'B'_sam_dna16);
    EXPECT_EQ(seqan3::complement('T'_sam_dna16), 'A'_sam_dna16);
    EXPECT_EQ(seqan3::complement('W'_sam_dna16), 'W'_sam_dna16);
    EXPECT_EQ(seqan3::complement('Y'_sam_dna16), 'R'_sam_dna16);
    EXPECT_EQ(seqan3::complement('H'_sam_dna16), 'D'_sam_dna16);
    EXPECT_EQ(seqan3::complement('K'_sam_dna16), 'M'_sam_dna16);
    EXPECT_EQ(seqan3::complement('D'_sam_dna16), 'H'_sam_dna16);
    EXPECT_EQ(seqan3::complement('B'_sam_dna16), 'V'_sam_dna16);
    EXPECT_EQ(seqan3::complement('N'_sam_dna16), 'N'_sam_dna16);
}

TEST(sam_dna16, to_char_assign_char)
{
    using rank_t = seqan3::alphabet_rank_t<seqan3::sam_dna16>;
    for (rank_t rank = 0; rank < seqan3::alphabet_size<seqan3::sam_dna16>; ++rank)
    {
        char chr = seqan3::to_char(seqan3::assign_rank_to(rank, seqan3::sam_dna16{}));
        EXPECT_EQ(seqan3::to_char(seqan3::sam_dna16{}.assign_char(chr)), chr);
    }

    EXPECT_EQ(seqan3::to_char(seqan3::sam_dna16{}.assign_char('a')), 'A');
    EXPECT_EQ(seqan3::to_char(seqan3::sam_dna16{}.assign_char('c')), 'C');
    EXPECT_EQ(seqan3::to_char(seqan3::sam_dna16{}.assign_char('g')), 'G');
    EXPECT_EQ(seqan3::to_char(seqan3::sam_dna16{}.assign_char('t')), 'T');

    EXPECT_EQ(seqan3::to_char(seqan3::sam_dna16{}.assign_char('U')), 'T');
    EXPECT_EQ(seqan3::to_char(seqan3::sam_dna16{}.assign_char('!')), 'N');
}

TEST(sam_dna16, char_literal)
{
    EXPECT_EQ(seqan3::to_char('A'_sam_dna16), 'A');
    EXPECT_EQ(seqan3::to_char('C'_sam_dna16), 'C');
    EXPECT_EQ(seqan3::to_char('G'_sam_dna16), 'G');

    EXPECT_EQ(seqan3::to_char('U'_sam_dna16), 'T');
    EXPECT_EQ(seqan3::to_char('T'_sam_dna16), 'T');

    EXPECT_EQ(seqan3::to_char('R'_sam_dna16), 'R');
    EXPECT_EQ(seqan3::to_char('Y'_sam_dna16), 'Y');
    EXPECT_EQ(seqan3::to_char('S'_sam_dna16), 'S');
    EXPECT_EQ(seqan3::to_char('W'_sam_dna16), 'W');
    EXPECT_EQ(seqan3::to_char('K'_sam_dna16), 'K');
    EXPECT_EQ(seqan3::to_char('M'_sam_dna16), 'M');
    EXPECT_EQ(seqan3::to_char('B'_sam_dna16), 'B');
    EXPECT_EQ(seqan3::to_char('D'_sam_dna16), 'D');
    EXPECT_EQ(seqan3::to_char('H'_sam_dna16), 'H');
    EXPECT_EQ(seqan3::to_char('V'_sam_dna16), 'V');

    EXPECT_EQ(seqan3::to_char('='_sam_dna16), '=');

    EXPECT_EQ(seqan3::to_char('N'_sam_dna16), 'N');
    EXPECT_EQ(seqan3::to_char('!'_sam_dna16), 'N');
}

TEST(sam_dna16, string_literal)
{
    seqan3::sam_dna16_vector v;
    v.resize(5, 'A'_sam_dna16);
    EXPECT_EQ(v, "AAAAA"_sam_dna16);

    std::vector<seqan3::sam_dna16> w{'A'_sam_dna16,
                                     '='_sam_dna16,
                                     'G'_sam_dna16,
                                     'T'_sam_dna16,
                                     'U'_sam_dna16,
                                     'N'_sam_dna16};
    EXPECT_EQ(w, "A=GTTN"_sam_dna16);
}

TEST(sam_dna16, char_is_valid)
{
    constexpr auto validator = seqan3::is_char<'A'> || seqan3::is_char<'C'> || seqan3::is_char<'G'> ||
                               seqan3::is_char<'T'> || seqan3::is_char<'U'> || seqan3::is_char<'a'> ||
                               seqan3::is_char<'c'> || seqan3::is_char<'g'> || seqan3::is_char<'t'> ||
                               seqan3::is_char<'u'> || seqan3::is_char<'N'> || seqan3::is_char<'n'> ||
                               seqan3::is_char<'R'> || seqan3::is_char<'Y'> || seqan3::is_char<'S'> ||
                               seqan3::is_char<'W'> || seqan3::is_char<'K'> || seqan3::is_char<'M'> ||
                               seqan3::is_char<'B'> || seqan3::is_char<'D'> || seqan3::is_char<'H'> ||
                               seqan3::is_char<'V'> || seqan3::is_char<'r'> || seqan3::is_char<'y'> ||
                               seqan3::is_char<'s'> || seqan3::is_char<'w'> || seqan3::is_char<'k'> ||
                               seqan3::is_char<'m'> || seqan3::is_char<'b'> || seqan3::is_char<'d'> ||
                               seqan3::is_char<'h'> || seqan3::is_char<'v'> || seqan3::is_char<'='>;

    for (char c : std::views::iota(std::numeric_limits<char>::min(), std::numeric_limits<char>::max()))
        EXPECT_EQ(seqan3::sam_dna16::char_is_valid(c), validator(c));
}
