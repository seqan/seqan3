// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/cigar/cigar.hpp>

#include "../alphabet_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"
#include "../semi_alphabet_constexpr_test_template.hpp"

INSTANTIATE_TYPED_TEST_SUITE_P(cigar, semi_alphabet_test, seqan3::cigar, );
INSTANTIATE_TYPED_TEST_SUITE_P(cigar, semi_alphabet_constexpr, seqan3::cigar, );

using namespace seqan3;

TEST(cigar, brace_init)
{
    cigar c1{uint32_t{223}, 'M'_cigar_op};
    EXPECT_EQ(c1.to_string(), small_string<11>{"223M"});
}

TEST(cigar, to_string)
{
    cigar c1{};

    assign_rank_to(uint32_t{223}, get<0>(c1));
    assign_char_to('M', get<1>(c1));
    EXPECT_EQ(c1.to_string(), small_string<11>{"223M"});
}

TEST(cigar, assign_string)
{
    cigar c1{};
    c1.assign_string("223M");
    EXPECT_EQ(uint32_t{223}, to_rank(get<0>(c1)));
    EXPECT_EQ('M',           get<1>(c1).to_char());
}
