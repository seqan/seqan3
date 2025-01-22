// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <string>

#include <seqan3/test/expect_same_type.hpp>
#include <seqan3/utility/type_traits/function_traits.hpp>

std::function test_function_object = [](size_t arg1, std::string & arg2)
{
    assert(arg1 < arg2.size());
    return arg2[arg1];
};

using function_ptr_t = std::string (*)(int, double const &&, bool &);

TEST(function_traits, argument_count)
{
    using function_t = decltype(test_function_object);
    EXPECT_EQ(seqan3::function_traits<function_t>::argument_count, 2u);
    EXPECT_EQ(seqan3::function_traits<function_ptr_t>::argument_count, 3u);
}

TEST(function_traits, result_type)
{
    using function_t = decltype(test_function_object);
    EXPECT_SAME_TYPE(seqan3::function_traits<function_t>::result_type, char);
    EXPECT_SAME_TYPE(seqan3::function_traits<function_ptr_t>::result_type, std::string);
}

TEST(function_traits, argument_type_at)
{
    using function_t = decltype(test_function_object);
    EXPECT_SAME_TYPE(seqan3::function_traits<function_t>::argument_type_at<0>, size_t);
    EXPECT_SAME_TYPE(seqan3::function_traits<function_t>::argument_type_at<1>, std::string &);
    EXPECT_SAME_TYPE(seqan3::function_traits<function_ptr_t>::argument_type_at<0>, int);
    EXPECT_SAME_TYPE(seqan3::function_traits<function_ptr_t>::argument_type_at<1>, double const &&);
    EXPECT_SAME_TYPE(seqan3::function_traits<function_ptr_t>::argument_type_at<2>, bool &);
}
