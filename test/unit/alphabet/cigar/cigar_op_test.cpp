// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include "../alphabet_test_template.hpp"
#include "../alphabet_constexpr_test_template.hpp"
#include <seqan3/alphabet/cigar/cigar_op.hpp>

INSTANTIATE_TYPED_TEST_CASE_P(cigar_op, alphabet, cigar_op);
INSTANTIATE_TYPED_TEST_CASE_P(cigar_op, alphabet_constexpr, cigar_op);

TEST(cigar_op, to_char_assign_char)
{
    for (char chr : std::string{"MDISHNPX="})
        EXPECT_EQ(to_char(cigar_op{}.assign_char(chr)), chr);
}

TEST(cigar_op, char_literal)
{
    EXPECT_EQ(to_char(cigar_op{'M'_cigar_op}), 'M');
    EXPECT_EQ(to_char(cigar_op{'D'_cigar_op}), 'D');
    EXPECT_EQ(to_char(cigar_op{'I'_cigar_op}), 'I');
    EXPECT_EQ(to_char(cigar_op{'S'_cigar_op}), 'S');
    EXPECT_EQ(to_char(cigar_op{'H'_cigar_op}), 'H');
    EXPECT_EQ(to_char(cigar_op{'N'_cigar_op}), 'N');
    EXPECT_EQ(to_char(cigar_op{'P'_cigar_op}), 'P');
    EXPECT_EQ(to_char(cigar_op{'X'_cigar_op}), 'X');
    EXPECT_EQ(to_char(cigar_op{'='_cigar_op}), '=');
}

TEST(cigar_op, assign_char_strictly_to)
{
    EXPECT_THROW(assign_char_strictly_to('A', cigar_op{}), invalid_char_assignment);
}
