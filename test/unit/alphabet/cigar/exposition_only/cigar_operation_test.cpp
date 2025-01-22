// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alphabet/cigar/exposition_only/cigar_operation.hpp>

#include "../../alphabet_constexpr_test_template.hpp"
#include "../../alphabet_test_template.hpp"
#include "../../semi_alphabet_constexpr_test_template.hpp"
#include "../../semi_alphabet_test_template.hpp"

INSTANTIATE_TYPED_TEST_SUITE_P(cigar_operation, alphabet, seqan3::exposition_only::cigar_operation, );
INSTANTIATE_TYPED_TEST_SUITE_P(cigar_operation, semi_alphabet_test, seqan3::exposition_only::cigar_operation, );
INSTANTIATE_TYPED_TEST_SUITE_P(cigar_operation, alphabet_constexpr, seqan3::exposition_only::cigar_operation, );
INSTANTIATE_TYPED_TEST_SUITE_P(cigar_operation, semi_alphabet_constexpr, seqan3::exposition_only::cigar_operation, );

TEST(cigar_operation, to_char_assign_char)
{
    for (char chr : std::string{"MDISHNPX="})
        EXPECT_EQ(seqan3::to_char(seqan3::exposition_only::cigar_operation{}.assign_char(chr)), chr);
}

TEST(cigar_operation, assign_char_strictly_to)
{
    EXPECT_THROW(seqan3::assign_char_strictly_to('A', seqan3::exposition_only::cigar_operation{}),
                 seqan3::invalid_char_assignment);
}
