// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string>

#include <seqan3/core/type_traits/basic.hpp>
#include <seqan3/core/type_traits/function.hpp>

constexpr
int    constexpr_nonvoid_free_fun(int i) { return i; }
int nonconstexpr_nonvoid_free_fun(int i) { return i; }

constexpr
int constexpr_nonvoid_free_fun_const_ref(int const & i) { return i; }
int nonconstexpr_nonvoid_free_fun_const_ref(int const & i) { return i; }

constexpr
void    constexpr_void_free_fun(int) { return; }
void nonconstexpr_void_free_fun(int) { return; }

struct constexpr_nonvoid_member_t
{
    int constexpr get_i(int i)
    {
        return i;
    }
};

struct constexpr_void_member_t
{
    void constexpr get_i(int)
    {}
};

struct nonconstexpr_nonvoid_member_t
{
    int get_i(int i)
    {
        return i;
    }
};

struct nonconstexpr_void_member_t
{
    void get_i(int)
    {}
};

TEST(type_trait, is_constexpr_invocable)
{
    int i = 32;
    int constexpr j = 42;

    EXPECT_TRUE((SEQAN3_IS_CONSTEXPR(constexpr_nonvoid_free_fun(3))));
    EXPECT_TRUE((SEQAN3_IS_CONSTEXPR(constexpr_nonvoid_free_fun(j))));
    EXPECT_TRUE((!SEQAN3_IS_CONSTEXPR(constexpr_nonvoid_free_fun(i))));
    EXPECT_TRUE((!SEQAN3_IS_CONSTEXPR(nonconstexpr_nonvoid_free_fun(3))));

    EXPECT_TRUE((SEQAN3_IS_CONSTEXPR(constexpr_nonvoid_free_fun_const_ref(static_cast<int const &>(3)))));
    EXPECT_TRUE((SEQAN3_IS_CONSTEXPR(constexpr_nonvoid_free_fun_const_ref(j))));
    EXPECT_TRUE((!SEQAN3_IS_CONSTEXPR(constexpr_nonvoid_free_fun_const_ref(i))));
    EXPECT_TRUE((!SEQAN3_IS_CONSTEXPR(nonconstexpr_nonvoid_free_fun_const_ref(static_cast<int const &>(3)))));

    EXPECT_TRUE((SEQAN3_IS_CONSTEXPR(constexpr_void_free_fun(3))));
    EXPECT_TRUE((SEQAN3_IS_CONSTEXPR(constexpr_void_free_fun(j))));
    EXPECT_TRUE((!SEQAN3_IS_CONSTEXPR(constexpr_void_free_fun(i))));
    EXPECT_TRUE((!SEQAN3_IS_CONSTEXPR(nonconstexpr_void_free_fun(3))));

    EXPECT_TRUE((SEQAN3_IS_CONSTEXPR(constexpr_nonvoid_member_t{}.get_i(3))));
    EXPECT_TRUE((SEQAN3_IS_CONSTEXPR(constexpr_void_member_t{}.get_i(3))));
    EXPECT_TRUE((!SEQAN3_IS_CONSTEXPR(nonconstexpr_nonvoid_member_t{}.get_i(3))));
    EXPECT_TRUE((!SEQAN3_IS_CONSTEXPR(nonconstexpr_void_member_t{}.get_i(3))));
}

std::function test_function_object = [] (size_t arg1, std::string & arg2)
{
    assert(arg1 < arg2.size());
    return arg2[arg1];
};

TEST(function_traits, argument_count)
{
    using function_t = decltype(test_function_object);
    EXPECT_EQ(seqan3::function_traits<function_t>::argument_count, 2u);
}

TEST(function_traits, result_type)
{
    using function_t = decltype(test_function_object);
    EXPECT_TRUE((std::same_as<seqan3::function_traits<function_t>::result_type, char>));
}

TEST(function_traits, argument_type_at)
{
    using function_t = decltype(test_function_object);
    EXPECT_TRUE((std::same_as<seqan3::function_traits<function_t>::argument_type_at<0>, size_t>));
    EXPECT_TRUE((std::same_as<seqan3::function_traits<function_t>::argument_type_at<1>, std::string &>));
}
