// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include "../alphabet_test_template.hpp"
#include "../alphabet_constexpr_test_template.hpp"
#include "nucleotide_test_template.hpp"

INSTANTIATE_TYPED_TEST_CASE_P(rna15, alphabet, rna15);
INSTANTIATE_TYPED_TEST_CASE_P(rna15, alphabet_constexpr, rna15);
INSTANTIATE_TYPED_TEST_CASE_P(rna15, nucleotide, rna15);

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
