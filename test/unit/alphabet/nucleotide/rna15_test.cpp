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

INSTANTIATE_TYPED_TEST_SUITE_P(rna15, alphabet_, rna15, );
INSTANTIATE_TYPED_TEST_SUITE_P(rna15, semi_alphabet_test, rna15, );
INSTANTIATE_TYPED_TEST_SUITE_P(rna15, alphabet_constexpr, rna15, );
INSTANTIATE_TYPED_TEST_SUITE_P(rna15, semi_alphabet_constexpr, rna15, );
INSTANTIATE_TYPED_TEST_SUITE_P(rna15, nucleotide, rna15, );

TEST(rna15, to_char_assign_char)
{
    EXPECT_EQ(to_char(rna15{}.assign_char('A')), 'A');
    EXPECT_EQ(to_char(rna15{}.assign_char('C')), 'C');
    EXPECT_EQ(to_char(rna15{}.assign_char('G')), 'G');

    EXPECT_EQ(to_char(rna15{}.assign_char('U')), 'U');
    EXPECT_EQ(to_char(rna15{}.assign_char('T')), 'U');

    EXPECT_EQ(to_char(rna15{}.assign_char('R')), 'R');
    EXPECT_EQ(to_char(rna15{}.assign_char('Y')), 'Y');
    EXPECT_EQ(to_char(rna15{}.assign_char('S')), 'S');
    EXPECT_EQ(to_char(rna15{}.assign_char('W')), 'W');
    EXPECT_EQ(to_char(rna15{}.assign_char('K')), 'K');
    EXPECT_EQ(to_char(rna15{}.assign_char('M')), 'M');
    EXPECT_EQ(to_char(rna15{}.assign_char('B')), 'B');
    EXPECT_EQ(to_char(rna15{}.assign_char('D')), 'D');
    EXPECT_EQ(to_char(rna15{}.assign_char('H')), 'H');
    EXPECT_EQ(to_char(rna15{}.assign_char('V')), 'V');

    EXPECT_EQ(to_char(rna15{}.assign_char('N')), 'N');
    EXPECT_EQ(to_char(rna15{}.assign_char('!')), 'N');
}

TEST(rna15, char_literal)
{
    EXPECT_EQ(to_char('A'_rna15), 'A');
    EXPECT_EQ(to_char('C'_rna15), 'C');
    EXPECT_EQ(to_char('G'_rna15), 'G');

    EXPECT_EQ(to_char('U'_rna15), 'U');
    EXPECT_EQ(to_char('T'_rna15), 'U');

    EXPECT_EQ(to_char('R'_rna15), 'R');
    EXPECT_EQ(to_char('Y'_rna15), 'Y');
    EXPECT_EQ(to_char('S'_rna15), 'S');
    EXPECT_EQ(to_char('W'_rna15), 'W');
    EXPECT_EQ(to_char('K'_rna15), 'K');
    EXPECT_EQ(to_char('M'_rna15), 'M');
    EXPECT_EQ(to_char('B'_rna15), 'B');
    EXPECT_EQ(to_char('D'_rna15), 'D');
    EXPECT_EQ(to_char('H'_rna15), 'H');
    EXPECT_EQ(to_char('V'_rna15), 'V');

    EXPECT_EQ(to_char('N'_rna15), 'N');
    EXPECT_EQ(to_char('!'_rna15), 'N');
}

TEST(rna15, string_literal)
{
    rna15_vector v;
    v.resize(5, 'A'_rna15);
    EXPECT_EQ(v, "AAAAA"_rna15);

    std::vector<rna15> w{'A'_rna15, 'C'_rna15, 'G'_rna15, 'T'_rna15, 'U'_rna15, 'N'_rna15};
    EXPECT_EQ(w, "ACGTTN"_rna15);
}

TEST(rna15, char_is_valid)
{
    constexpr auto validator = is_char<'A'> || is_char<'C'> || is_char<'G'> || is_char<'T'> || is_char<'U'>
                            || is_char<'a'> || is_char<'c'> || is_char<'g'> || is_char<'t'> || is_char<'u'>
                            || is_char<'N'> || is_char<'n'>
                            || is_char<'R'> || is_char<'Y'> || is_char<'S'> || is_char<'W'> || is_char<'K'>
                            || is_char<'M'> || is_char<'B'> || is_char<'D'> || is_char<'H'> || is_char<'V'>
                            || is_char<'r'> || is_char<'y'> || is_char<'s'> || is_char<'w'> || is_char<'k'>
                            || is_char<'m'> || is_char<'b'> || is_char<'d'> || is_char<'h'> || is_char<'v'>;
    for (char c : std::views::iota(std::numeric_limits<char>::min(), std::numeric_limits<char>::max()))
        EXPECT_EQ(rna15::char_is_valid(c), validator(c));
}
