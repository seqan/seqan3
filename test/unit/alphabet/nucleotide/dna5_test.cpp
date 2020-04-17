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
#include "nucleotide_test_template.hpp"

using seqan3::operator""_dna5;

INSTANTIATE_TYPED_TEST_SUITE_P(dna5, alphabet_, seqan3::dna5, );
INSTANTIATE_TYPED_TEST_SUITE_P(dna5, semi_alphabet_test, seqan3::dna5, );
INSTANTIATE_TYPED_TEST_SUITE_P(dna5, alphabet_constexpr, seqan3::dna5, );
INSTANTIATE_TYPED_TEST_SUITE_P(dna5, semi_alphabet_constexpr, seqan3::dna5, );
INSTANTIATE_TYPED_TEST_SUITE_P(dna5, nucleotide, seqan3::dna5, );

TEST(dna5, to_char_assign_char)
{
    EXPECT_EQ(seqan3::to_char(seqan3::dna5{}.assign_char('A')), 'A');
    EXPECT_EQ(seqan3::to_char(seqan3::dna5{}.assign_char('C')), 'C');
    EXPECT_EQ(seqan3::to_char(seqan3::dna5{}.assign_char('G')), 'G');

    EXPECT_EQ(seqan3::to_char(seqan3::dna5{}.assign_char('U')), 'T');
    EXPECT_EQ(seqan3::to_char(seqan3::dna5{}.assign_char('T')), 'T');

    EXPECT_EQ(seqan3::to_char(seqan3::dna5{}.assign_char('R')), 'N');
    EXPECT_EQ(seqan3::to_char(seqan3::dna5{}.assign_char('Y')), 'N');
    EXPECT_EQ(seqan3::to_char(seqan3::dna5{}.assign_char('S')), 'N');
    EXPECT_EQ(seqan3::to_char(seqan3::dna5{}.assign_char('W')), 'N');
    EXPECT_EQ(seqan3::to_char(seqan3::dna5{}.assign_char('K')), 'N');
    EXPECT_EQ(seqan3::to_char(seqan3::dna5{}.assign_char('M')), 'N');
    EXPECT_EQ(seqan3::to_char(seqan3::dna5{}.assign_char('B')), 'N');
    EXPECT_EQ(seqan3::to_char(seqan3::dna5{}.assign_char('D')), 'N');
    EXPECT_EQ(seqan3::to_char(seqan3::dna5{}.assign_char('H')), 'N');
    EXPECT_EQ(seqan3::to_char(seqan3::dna5{}.assign_char('V')), 'N');

    EXPECT_EQ(seqan3::to_char(seqan3::dna5{}.assign_char('N')), 'N');
    EXPECT_EQ(seqan3::to_char(seqan3::dna5{}.assign_char('!')), 'N');
}

TEST(dna5, char_literal)
{
    EXPECT_EQ(seqan3::to_char('A'_dna5), 'A');
    EXPECT_EQ(seqan3::to_char('C'_dna5), 'C');
    EXPECT_EQ(seqan3::to_char('G'_dna5), 'G');

    EXPECT_EQ(seqan3::to_char('U'_dna5), 'T');
    EXPECT_EQ(seqan3::to_char('T'_dna5), 'T');

    EXPECT_EQ(seqan3::to_char('R'_dna5), 'N');
    EXPECT_EQ(seqan3::to_char('Y'_dna5), 'N');
    EXPECT_EQ(seqan3::to_char('S'_dna5), 'N');
    EXPECT_EQ(seqan3::to_char('W'_dna5), 'N');
    EXPECT_EQ(seqan3::to_char('K'_dna5), 'N');
    EXPECT_EQ(seqan3::to_char('M'_dna5), 'N');
    EXPECT_EQ(seqan3::to_char('B'_dna5), 'N');
    EXPECT_EQ(seqan3::to_char('D'_dna5), 'N');
    EXPECT_EQ(seqan3::to_char('H'_dna5), 'N');
    EXPECT_EQ(seqan3::to_char('V'_dna5), 'N');

    EXPECT_EQ(seqan3::to_char('N'_dna5), 'N');
    EXPECT_EQ(seqan3::to_char('!'_dna5), 'N');
}

TEST(dna5, string_literal)
{
    seqan3::dna5_vector v;
    v.resize(5, 'A'_dna5);
    EXPECT_EQ(v, "AAAAA"_dna5);

    std::vector<seqan3::dna5> w{'A'_dna5, 'C'_dna5, 'G'_dna5, 'T'_dna5, 'U'_dna5, 'N'_dna5};
    EXPECT_EQ(w, "ACGTTN"_dna5);
}

TEST(dna5, char_is_valid)
{
    constexpr auto validator = seqan3::is_char<'A'> || seqan3::is_char<'C'> || seqan3::is_char<'G'> ||
                               seqan3::is_char<'T'> || seqan3::is_char<'U'> || seqan3::is_char<'a'> ||
                               seqan3::is_char<'c'> || seqan3::is_char<'g'> || seqan3::is_char<'t'> ||
                               seqan3::is_char<'u'> || seqan3::is_char<'N'> || seqan3::is_char<'n'>;
    for (char c : std::views::iota(std::numeric_limits<char>::min(), std::numeric_limits<char>::max()))
        EXPECT_EQ(seqan3::dna5::char_is_valid(c), validator(c));
}
