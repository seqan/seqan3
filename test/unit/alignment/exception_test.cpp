// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

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
