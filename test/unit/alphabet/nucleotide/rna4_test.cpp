// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include "../alphabet_constexpr_test_template.hpp"
#include "../alphabet_test_template.hpp"
#include "../semi_alphabet_constexpr_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"
#include "nucleotide_test_template.hpp"

INSTANTIATE_TYPED_TEST_SUITE_P(rna4, alphabet_, rna4, );
INSTANTIATE_TYPED_TEST_SUITE_P(rna4, semi_alphabet_test, rna4, );
INSTANTIATE_TYPED_TEST_SUITE_P(rna4, alphabet_constexpr, rna4, );
INSTANTIATE_TYPED_TEST_SUITE_P(rna4, semi_alphabet_constexpr, rna4, );
INSTANTIATE_TYPED_TEST_SUITE_P(rna4, nucleotide, rna4, );

TEST(rna4, to_char_assign_char)
{
    EXPECT_EQ(to_char(rna4{}.assign_char('A')), 'A');
    EXPECT_EQ(to_char(rna4{}.assign_char('C')), 'C');
    EXPECT_EQ(to_char(rna4{}.assign_char('G')), 'G');

    EXPECT_EQ(to_char(rna4{}.assign_char('U')), 'U');
    EXPECT_EQ(to_char(rna4{}.assign_char('T')), 'U');

    EXPECT_EQ(to_char(rna4{}.assign_char('R')), 'A');
    EXPECT_EQ(to_char(rna4{}.assign_char('Y')), 'C');
    EXPECT_EQ(to_char(rna4{}.assign_char('S')), 'C');
    EXPECT_EQ(to_char(rna4{}.assign_char('W')), 'A');
    EXPECT_EQ(to_char(rna4{}.assign_char('K')), 'G');
    EXPECT_EQ(to_char(rna4{}.assign_char('M')), 'A');
    EXPECT_EQ(to_char(rna4{}.assign_char('B')), 'C');
    EXPECT_EQ(to_char(rna4{}.assign_char('D')), 'A');
    EXPECT_EQ(to_char(rna4{}.assign_char('H')), 'A');
    EXPECT_EQ(to_char(rna4{}.assign_char('V')), 'A');

    EXPECT_EQ(to_char(rna4{}.assign_char('N')), 'A');
    EXPECT_EQ(to_char(rna4{}.assign_char('!')), 'A');
}

TEST(rna4, char_literal)
{
    EXPECT_EQ(to_char('A'_rna4), 'A');
    EXPECT_EQ(to_char('C'_rna4), 'C');
    EXPECT_EQ(to_char('G'_rna4), 'G');

    EXPECT_EQ(to_char('U'_rna4), 'U');
    EXPECT_EQ(to_char('T'_rna4), 'U');

    EXPECT_EQ(to_char('R'_rna4), 'A');
    EXPECT_EQ(to_char('Y'_rna4), 'C');
    EXPECT_EQ(to_char('S'_rna4), 'C');
    EXPECT_EQ(to_char('W'_rna4), 'A');
    EXPECT_EQ(to_char('K'_rna4), 'G');
    EXPECT_EQ(to_char('M'_rna4), 'A');
    EXPECT_EQ(to_char('B'_rna4), 'C');
    EXPECT_EQ(to_char('D'_rna4), 'A');
    EXPECT_EQ(to_char('H'_rna4), 'A');
    EXPECT_EQ(to_char('V'_rna4), 'A');

    EXPECT_EQ(to_char('N'_rna4), 'A');
    EXPECT_EQ(to_char('!'_rna4), 'A');
}

TEST(rna4, string_literal)
{
    rna4_vector v;
    v.resize(5, 'A'_rna4);
    EXPECT_EQ(v, "AAAAA"_rna4);

    std::vector<rna4> w{'A'_rna4, 'C'_rna4, 'G'_rna4, 'T'_rna4, 'U'_rna4, 'N'_rna4};
    EXPECT_EQ(w, "ACGUUA"_rna4);
}

TEST(rna4, char_is_valid)
{
    constexpr auto validator = is_char<'A'> || is_char<'C'> || is_char<'G'> || is_char<'T'> || is_char<'U'>
                            || is_char<'a'> || is_char<'c'> || is_char<'g'> || is_char<'t'> || is_char<'u'>;
    for (char c : std::views::iota(std::numeric_limits<char>::min(), std::numeric_limits<char>::max()))
        EXPECT_EQ(rna4::char_is_valid(c), validator(c));
}
