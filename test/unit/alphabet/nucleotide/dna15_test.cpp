// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <range/v3/view/zip.hpp>

#include "../alphabet_constexpr_test_template.hpp"
#include "../alphabet_test_template.hpp"
#include "../semi_alphabet_constexpr_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"
#include "nucleotide_test_template.hpp"

INSTANTIATE_TYPED_TEST_SUITE_P(dna15, alphabet_, dna15, );
INSTANTIATE_TYPED_TEST_SUITE_P(dna15, semi_alphabet_test, dna15, );
INSTANTIATE_TYPED_TEST_SUITE_P(dna15, alphabet_constexpr, dna15, );
INSTANTIATE_TYPED_TEST_SUITE_P(dna15, semi_alphabet_constexpr, dna15, );
INSTANTIATE_TYPED_TEST_SUITE_P(dna15, nucleotide, dna15, );

TEST(dna15, to_char_assign_char)
{
    EXPECT_EQ(to_char(dna15{}.assign_char('A')), 'A');
    EXPECT_EQ(to_char(dna15{}.assign_char('C')), 'C');
    EXPECT_EQ(to_char(dna15{}.assign_char('G')), 'G');

    EXPECT_EQ(to_char(dna15{}.assign_char('U')), 'T');
    EXPECT_EQ(to_char(dna15{}.assign_char('T')), 'T');

    EXPECT_EQ(to_char(dna15{}.assign_char('R')), 'R');
    EXPECT_EQ(to_char(dna15{}.assign_char('Y')), 'Y');
    EXPECT_EQ(to_char(dna15{}.assign_char('S')), 'S');
    EXPECT_EQ(to_char(dna15{}.assign_char('W')), 'W');
    EXPECT_EQ(to_char(dna15{}.assign_char('K')), 'K');
    EXPECT_EQ(to_char(dna15{}.assign_char('M')), 'M');
    EXPECT_EQ(to_char(dna15{}.assign_char('B')), 'B');
    EXPECT_EQ(to_char(dna15{}.assign_char('D')), 'D');
    EXPECT_EQ(to_char(dna15{}.assign_char('H')), 'H');
    EXPECT_EQ(to_char(dna15{}.assign_char('V')), 'V');

    EXPECT_EQ(to_char(dna15{}.assign_char('N')), 'N');
    EXPECT_EQ(to_char(dna15{}.assign_char('!')), 'N');
}

TEST(dna15, char_literal)
{
    EXPECT_EQ(to_char('A'_dna15), 'A');
    EXPECT_EQ(to_char('C'_dna15), 'C');
    EXPECT_EQ(to_char('G'_dna15), 'G');

    EXPECT_EQ(to_char('U'_dna15), 'T');
    EXPECT_EQ(to_char('T'_dna15), 'T');

    EXPECT_EQ(to_char('R'_dna15), 'R');
    EXPECT_EQ(to_char('Y'_dna15), 'Y');
    EXPECT_EQ(to_char('S'_dna15), 'S');
    EXPECT_EQ(to_char('W'_dna15), 'W');
    EXPECT_EQ(to_char('K'_dna15), 'K');
    EXPECT_EQ(to_char('M'_dna15), 'M');
    EXPECT_EQ(to_char('B'_dna15), 'B');
    EXPECT_EQ(to_char('D'_dna15), 'D');
    EXPECT_EQ(to_char('H'_dna15), 'H');
    EXPECT_EQ(to_char('V'_dna15), 'V');

    EXPECT_EQ(to_char('N'_dna15), 'N');
    EXPECT_EQ(to_char('!'_dna15), 'N');
}

TEST(dna15, string_literal)
{
    dna15_vector v;
    v.resize(5, 'A'_dna15);
    EXPECT_EQ(v, "AAAAA"_dna15);

    std::vector<dna15> w{'A'_dna15, 'C'_dna15, 'G'_dna15, 'T'_dna15, 'U'_dna15, 'N'_dna15};
    EXPECT_EQ(w, "ACGTTN"_dna15);
}

TEST(dna15, char_is_valid)
{
    constexpr auto validator = is_char<'A'> || is_char<'C'> || is_char<'G'> || is_char<'T'> || is_char<'U'>
                            || is_char<'a'> || is_char<'c'> || is_char<'g'> || is_char<'t'> || is_char<'u'>
                            || is_char<'N'> || is_char<'n'>
                            || is_char<'R'> || is_char<'Y'> || is_char<'S'> || is_char<'W'> || is_char<'K'>
                            || is_char<'M'> || is_char<'B'> || is_char<'D'> || is_char<'H'> || is_char<'V'>
                            || is_char<'r'> || is_char<'y'> || is_char<'s'> || is_char<'w'> || is_char<'k'>
                            || is_char<'m'> || is_char<'b'> || is_char<'d'> || is_char<'h'> || is_char<'v'>;
    for (char c : std::views::iota(std::numeric_limits<char>::min(), std::numeric_limits<char>::max()))
        EXPECT_EQ(dna15::char_is_valid(c), validator(c));
}
