// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>
#include <ostream>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/test/pretty_printing.hpp>

using namespace seqan3;

TEST(debug_stream, basic)
{
    testing::internal::CaptureStdout();

    PrintTo('a', &std::cout);

    PrintTo("AGA"_dna4, &std::cout);

    PrintTo(std::make_tuple<int, int>(42, -10), &std::cout);

    EXPECT_EQ((testing::internal::GetCapturedStdout()), "aAGA(42,-10)");
}
