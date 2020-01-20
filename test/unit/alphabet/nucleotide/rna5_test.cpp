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

INSTANTIATE_TYPED_TEST_SUITE_P(rna5, alphabet_, rna5, );
INSTANTIATE_TYPED_TEST_SUITE_P(rna5, semi_alphabet_test, rna5, );
INSTANTIATE_TYPED_TEST_SUITE_P(rna5, alphabet_constexpr, rna5, );
INSTANTIATE_TYPED_TEST_SUITE_P(rna5, semi_alphabet_constexpr, rna5, );
INSTANTIATE_TYPED_TEST_SUITE_P(rna5, nucleotide, rna5, );

TEST(rna5, to_char_assign_char)
{
    EXPECT_EQ(to_char(rna5{}.assign_char('A')), 'A');
    EXPECT_EQ(to_char(rna5{}.assign_char('C')), 'C');
    EXPECT_EQ(to_char(rna5{}.assign_char('G')), 'G');

    EXPECT_EQ(to_char(rna5{}.assign_char('U')), 'U');
    EXPECT_EQ(to_char(rna5{}.assign_char('T')), 'U');

    EXPECT_EQ(to_char(rna5{}.assign_char('R')), 'N');
    EXPECT_EQ(to_char(rna5{}.assign_char('Y')), 'N');
    EXPECT_EQ(to_char(rna5{}.assign_char('S')), 'N');
    EXPECT_EQ(to_char(rna5{}.assign_char('W')), 'N');
    EXPECT_EQ(to_char(rna5{}.assign_char('K')), 'N');
    EXPECT_EQ(to_char(rna5{}.assign_char('M')), 'N');
    EXPECT_EQ(to_char(rna5{}.assign_char('B')), 'N');
    EXPECT_EQ(to_char(rna5{}.assign_char('D')), 'N');
    EXPECT_EQ(to_char(rna5{}.assign_char('H')), 'N');
    EXPECT_EQ(to_char(rna5{}.assign_char('V')), 'N');

    EXPECT_EQ(to_char(rna5{}.assign_char('N')), 'N');
    EXPECT_EQ(to_char(rna5{}.assign_char('!')), 'N');
}

TEST(rna5, char_literal)
{
    EXPECT_EQ(to_char('A'_rna5), 'A');
    EXPECT_EQ(to_char('C'_rna5), 'C');
    EXPECT_EQ(to_char('G'_rna5), 'G');

    EXPECT_EQ(to_char('U'_rna5), 'U');
    EXPECT_EQ(to_char('T'_rna5), 'U');

    EXPECT_EQ(to_char('R'_rna5), 'N');
    EXPECT_EQ(to_char('Y'_rna5), 'N');
    EXPECT_EQ(to_char('S'_rna5), 'N');
    EXPECT_EQ(to_char('W'_rna5), 'N');
    EXPECT_EQ(to_char('K'_rna5), 'N');
    EXPECT_EQ(to_char('M'_rna5), 'N');
    EXPECT_EQ(to_char('B'_rna5), 'N');
    EXPECT_EQ(to_char('D'_rna5), 'N');
    EXPECT_EQ(to_char('H'_rna5), 'N');
    EXPECT_EQ(to_char('V'_rna5), 'N');

    EXPECT_EQ(to_char('N'_rna5), 'N');
    EXPECT_EQ(to_char('!'_rna5), 'N');
}

TEST(rna5, string_literal)
{
    rna5_vector v;
    v.resize(5, 'A'_rna5);
    EXPECT_EQ(v, "AAAAA"_rna5);

    std::vector<rna5> w{'A'_rna5, 'C'_rna5, 'G'_rna5, 'T'_rna5, 'U'_rna5, 'N'_rna5};
    EXPECT_EQ(w, "ACGUUN"_rna5);
}

TEST(rna5, char_is_valid)
{
    constexpr auto validator = is_char<'A'> || is_char<'C'> || is_char<'G'> || is_char<'T'> || is_char<'U'>
                            || is_char<'a'> || is_char<'c'> || is_char<'g'> || is_char<'t'> || is_char<'u'>
                            || is_char<'N'> || is_char<'n'>;
    for (char c : std::views::iota(std::numeric_limits<char>::min(), std::numeric_limits<char>::max()))
        EXPECT_EQ(rna5::char_is_valid(c), validator(c));
}
