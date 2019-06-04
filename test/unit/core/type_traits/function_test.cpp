// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

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
