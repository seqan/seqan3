// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alignment/exception.hpp>

TEST(alignment_exception, invalid_alignment_configuration)
{

    auto throw_lambda = []()
    {
        throw seqan3::invalid_alignment_configuration{"test"};
    };

    EXPECT_THROW(throw_lambda(), seqan3::invalid_alignment_configuration);

    try
    {
        throw_lambda();
    }
    catch (seqan3::invalid_alignment_configuration & ex)
    {
        EXPECT_EQ(ex.what(), std::string{"test"});
    }
    catch (...)
    {
        FAIL();
    }
}
