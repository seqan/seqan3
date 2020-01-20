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

INSTANTIATE_TYPED_TEST_SUITE_P(dna4, alphabet_, dna4, );
INSTANTIATE_TYPED_TEST_SUITE_P(dna4, semi_alphabet_test, dna4, );
INSTANTIATE_TYPED_TEST_SUITE_P(dna4, alphabet_constexpr, dna4, );
INSTANTIATE_TYPED_TEST_SUITE_P(dna4, semi_alphabet_constexpr, dna4, );
INSTANTIATE_TYPED_TEST_SUITE_P(dna4, nucleotide, dna4, );

TEST(dna4, to_char_assign_char)
{
    EXPECT_EQ(to_char(dna4{}.assign_char('A')), 'A');
    EXPECT_EQ(to_char(dna4{}.assign_char('C')), 'C');
    EXPECT_EQ(to_char(dna4{}.assign_char('G')), 'G');

    EXPECT_EQ(to_char(dna4{}.assign_char('U')), 'T');
    EXPECT_EQ(to_char(dna4{}.assign_char('T')), 'T');

    EXPECT_EQ(to_char(dna4{}.assign_char('R')), 'A');
    EXPECT_EQ(to_char(dna4{}.assign_char('Y')), 'C');
    EXPECT_EQ(to_char(dna4{}.assign_char('S')), 'C');
    EXPECT_EQ(to_char(dna4{}.assign_char('W')), 'A');
    EXPECT_EQ(to_char(dna4{}.assign_char('K')), 'G');
    EXPECT_EQ(to_char(dna4{}.assign_char('M')), 'A');
    EXPECT_EQ(to_char(dna4{}.assign_char('B')), 'C');
    EXPECT_EQ(to_char(dna4{}.assign_char('D')), 'A');
    EXPECT_EQ(to_char(dna4{}.assign_char('H')), 'A');
    EXPECT_EQ(to_char(dna4{}.assign_char('V')), 'A');

    EXPECT_EQ(to_char(dna4{}.assign_char('N')), 'A');
    EXPECT_EQ(to_char(dna4{}.assign_char('!')), 'A');
}

TEST(dna4, char_literal)
{
    EXPECT_EQ(to_char('A'_dna4), 'A');
    EXPECT_EQ(to_char('C'_dna4), 'C');
    EXPECT_EQ(to_char('G'_dna4), 'G');

    EXPECT_EQ(to_char('U'_dna4), 'T');
    EXPECT_EQ(to_char('T'_dna4), 'T');

    EXPECT_EQ(to_char('R'_dna4), 'A');
    EXPECT_EQ(to_char('Y'_dna4), 'C');
    EXPECT_EQ(to_char('S'_dna4), 'C');
    EXPECT_EQ(to_char('W'_dna4), 'A');
    EXPECT_EQ(to_char('K'_dna4), 'G');
    EXPECT_EQ(to_char('M'_dna4), 'A');
    EXPECT_EQ(to_char('B'_dna4), 'C');
    EXPECT_EQ(to_char('D'_dna4), 'A');
    EXPECT_EQ(to_char('H'_dna4), 'A');
    EXPECT_EQ(to_char('V'_dna4), 'A');

    EXPECT_EQ(to_char('N'_dna4), 'A');
    EXPECT_EQ(to_char('!'_dna4), 'A');
}

TEST(dna4, string_literal)
{
    dna4_vector v;
    v.resize(5, 'A'_dna4);
    EXPECT_EQ(v, "AAAAA"_dna4);

    std::vector<dna4> w{'A'_dna4, 'C'_dna4, 'G'_dna4, 'T'_dna4, 'U'_dna4, 'N'_dna4};
    EXPECT_EQ(w, "ACGTTA"_dna4);
}

TEST(dna4, char_is_valid)
{
    constexpr auto validator = is_char<'A'> || is_char<'C'> || is_char<'G'> || is_char<'T'> || is_char<'U'>
                            || is_char<'a'> || is_char<'c'> || is_char<'g'> || is_char<'t'> || is_char<'u'>;
    for (char c : std::views::iota(std::numeric_limits<char>::min(), std::numeric_limits<char>::max()))
        EXPECT_EQ(dna4::char_is_valid(c), validator(c));
}
