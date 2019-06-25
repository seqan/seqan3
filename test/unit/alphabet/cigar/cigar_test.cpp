// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/cigar/cigar.hpp>

#include "../alphabet_test_template.hpp"

INSTANTIATE_TYPED_TEST_CASE_P(cigar, alphabet, cigar);
// note: cigar is not usable in constexpr context.

TEST(cigar, brace_init)
{
    cigar c1{uint32_t{223}, 'M'_cigar_op};
    EXPECT_EQ(c1.to_char(), small_string<11>{"223M"});
}

TEST(cigar, to_char)
{
    cigar c1{};

    assign_rank_to(uint32_t{223}, get<0>(c1));
    assign_char_to('M', get<1>(c1));
    EXPECT_EQ(c1.to_char(), small_string<11>{"223M"});
}

TEST(cigar, assign_char)
{
    cigar c1{};
    c1.assign_char("223M");
    EXPECT_EQ(uint32_t{223}, to_rank(get<0>(c1)));
    EXPECT_EQ('M',           to_char(get<1>(c1)));
}

TEST(cigar, assign_char_strictly_to)
{
    EXPECT_THROW(assign_char_strictly_to("FOO", cigar{}), invalid_char_assignment);
    EXPECT_THROW(assign_char_strictly_to("223MZ", cigar{}), invalid_char_assignment);
}

TEST(cigar, debug_streaming)
{
    std::ostringstream o;
    debug_stream_type my_stream{o};

    cigar c1{};
    c1.assign_char("223M");

    my_stream << c1;

    o.flush();
    EXPECT_EQ(o.str().size(), 4u);
    EXPECT_EQ(o.str(), "223M");
}
