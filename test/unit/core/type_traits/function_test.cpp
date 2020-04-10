// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string>

#include <seqan3/core/type_traits/function.hpp>

std::function test_function_object = [] (size_t arg1, std::string & arg2)
{
    assert(arg1 < arg2.size());
    return arg2[arg1];
};

using function_ptr_t = std::string (*) (int, const double &&, bool &);

TEST(function_traits, argument_count)
{
    using function_t = decltype(test_function_object);
    EXPECT_EQ(seqan3::function_traits<function_t>::argument_count, 2u);
    EXPECT_EQ(seqan3::function_traits<function_ptr_t>::argument_count, 3u);
}

TEST(function_traits, result_type)
{
    using function_t = decltype(test_function_object);
    EXPECT_TRUE((std::same_as<seqan3::function_traits<function_t>::result_type, char>));
    EXPECT_TRUE((std::same_as<seqan3::function_traits<function_ptr_t>::result_type, std::string>));
}

TEST(function_traits, argument_type_at)
{
    using function_t = decltype(test_function_object);
    EXPECT_TRUE((std::same_as<seqan3::function_traits<function_t>::argument_type_at<0>, size_t>));
    EXPECT_TRUE((std::same_as<seqan3::function_traits<function_t>::argument_type_at<1>, std::string &>));
    EXPECT_TRUE((std::same_as<seqan3::function_traits<function_ptr_t>::argument_type_at<0>, int>));
    EXPECT_TRUE((std::same_as<seqan3::function_traits<function_ptr_t>::argument_type_at<1>, const double &&>));
    EXPECT_TRUE((std::same_as<seqan3::function_traits<function_ptr_t>::argument_type_at<2>, bool &>));
}
