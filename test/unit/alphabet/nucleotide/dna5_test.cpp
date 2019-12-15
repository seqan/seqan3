// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <range/v3/view/zip.hpp>

#include "../alphabet_constexpr_test_template.hpp"
#include "../alphabet_test_template.hpp"
#include "../semi_alphabet_constexpr_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"
#include "nucleotide_test_template.hpp"

INSTANTIATE_TYPED_TEST_SUITE_P(dna5, alphabet_, dna5, );
INSTANTIATE_TYPED_TEST_SUITE_P(dna5, semi_alphabet_test, dna5, );
INSTANTIATE_TYPED_TEST_SUITE_P(dna5, alphabet_constexpr, dna5, );
INSTANTIATE_TYPED_TEST_SUITE_P(dna5, semi_alphabet_constexpr, dna5, );
INSTANTIATE_TYPED_TEST_SUITE_P(dna5, nucleotide, dna5, );

TEST(dna5, to_char_assign_char)
{
    EXPECT_EQ(to_char(dna5{}.assign_char('A')), 'A');
    EXPECT_EQ(to_char(dna5{}.assign_char('C')), 'C');
    EXPECT_EQ(to_char(dna5{}.assign_char('G')), 'G');

    EXPECT_EQ(to_char(dna5{}.assign_char('U')), 'T');
    EXPECT_EQ(to_char(dna5{}.assign_char('T')), 'T');

    EXPECT_EQ(to_char(dna5{}.assign_char('R')), 'N');
    EXPECT_EQ(to_char(dna5{}.assign_char('Y')), 'N');
    EXPECT_EQ(to_char(dna5{}.assign_char('S')), 'N');
    EXPECT_EQ(to_char(dna5{}.assign_char('W')), 'N');
    EXPECT_EQ(to_char(dna5{}.assign_char('K')), 'N');
    EXPECT_EQ(to_char(dna5{}.assign_char('M')), 'N');
    EXPECT_EQ(to_char(dna5{}.assign_char('B')), 'N');
    EXPECT_EQ(to_char(dna5{}.assign_char('D')), 'N');
    EXPECT_EQ(to_char(dna5{}.assign_char('H')), 'N');
    EXPECT_EQ(to_char(dna5{}.assign_char('V')), 'N');

    EXPECT_EQ(to_char(dna5{}.assign_char('N')), 'N');
    EXPECT_EQ(to_char(dna5{}.assign_char('!')), 'N');
}

TEST(dna5, char_literal)
{
    EXPECT_EQ(to_char('A'_dna5), 'A');
    EXPECT_EQ(to_char('C'_dna5), 'C');
    EXPECT_EQ(to_char('G'_dna5), 'G');

    EXPECT_EQ(to_char('U'_dna5), 'T');
    EXPECT_EQ(to_char('T'_dna5), 'T');

    EXPECT_EQ(to_char('R'_dna5), 'N');
    EXPECT_EQ(to_char('Y'_dna5), 'N');
    EXPECT_EQ(to_char('S'_dna5), 'N');
    EXPECT_EQ(to_char('W'_dna5), 'N');
    EXPECT_EQ(to_char('K'_dna5), 'N');
    EXPECT_EQ(to_char('M'_dna5), 'N');
    EXPECT_EQ(to_char('B'_dna5), 'N');
    EXPECT_EQ(to_char('D'_dna5), 'N');
    EXPECT_EQ(to_char('H'_dna5), 'N');
    EXPECT_EQ(to_char('V'_dna5), 'N');

    EXPECT_EQ(to_char('N'_dna5), 'N');
    EXPECT_EQ(to_char('!'_dna5), 'N');
}

TEST(dna5, string_literal)
{
    dna5_vector v;
    v.resize(5, 'A'_dna5);
    EXPECT_EQ(v, "AAAAA"_dna5);

    std::vector<dna5> w{'A'_dna5, 'C'_dna5, 'G'_dna5, 'T'_dna5, 'U'_dna5, 'N'_dna5};
    EXPECT_EQ(w, "ACGTTN"_dna5);
}

TEST(dna5, char_is_valid)
{
    constexpr auto validator = is_char<'A'> || is_char<'C'> || is_char<'G'> || is_char<'T'> || is_char<'U'>
                            || is_char<'a'> || is_char<'c'> || is_char<'g'> || is_char<'t'> || is_char<'u'>
                            || is_char<'N'> || is_char<'n'>;
    for (char c : std::views::iota(std::numeric_limits<char>::min(), std::numeric_limits<char>::max()))
        EXPECT_EQ(dna5::char_is_valid(c), validator(c));
}
