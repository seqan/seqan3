// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <list>

#include <seqan3/utility/type_traits/basic.hpp>

//------------------------------------------------------------------------------
// Tests for remove_cvref transformation trait
//------------------------------------------------------------------------------

TEST(type_trait, remove_cvref_t)
{
    EXPECT_TRUE((std::is_same_v<int, std::remove_cvref_t<int>>));
    EXPECT_TRUE((std::is_same_v<int, std::remove_cvref_t<int const>>));
    EXPECT_TRUE((std::is_same_v<int, std::remove_cvref_t<volatile int>>));
    EXPECT_TRUE((std::is_same_v<int, std::remove_cvref_t<int &>>));
    EXPECT_TRUE((std::is_same_v<int, std::remove_cvref_t<int &&>>));
    EXPECT_TRUE((std::is_same_v<int, std::remove_cvref_t<int const &>>));
    EXPECT_TRUE((std::is_same_v<int, std::remove_cvref_t<int const &&>>));
    EXPECT_TRUE((std::is_same_v<int, std::remove_cvref_t<volatile int const &&>>));
    // don't decay pointers and arrays:
    EXPECT_FALSE((std::is_same_v<int, std::remove_cvref_t<int *>>));    // type stays same
    EXPECT_FALSE((std::is_same_v<int, std::remove_cvref_t<int[3]>>));   // type stays same
    EXPECT_FALSE((std::is_same_v<int *, std::remove_cvref_t<int[3]>>)); // type stays same
    // the last example would be true for std::decay_t
}

//------------------------------------------------------------------------------
// Tests for is_constexpr
//------------------------------------------------------------------------------

constexpr int constexpr_nonvoid_free_fun(int i)
{
    return i;
}
int nonconstexpr_nonvoid_free_fun(int i)
{
    return i;
}

constexpr int constexpr_nonvoid_free_fun_const_ref(int const & i)
{
    return i;
}
int nonconstexpr_nonvoid_free_fun_const_ref(int const & i)
{
    return i;
}

constexpr void constexpr_void_free_fun(int)
{
    return;
}
void nonconstexpr_void_free_fun(int)
{
    return;
}

struct constexpr_nonvoid_member_t
{
    constexpr int get_i(int i)
    {
        return i;
    }
};

struct constexpr_void_member_t
{
    constexpr void get_i(int)
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
#if defined(__clang__)
    GTEST_SKIP() << "Skipping on clang: https://github.com/llvm/llvm-project/issues/58078";
#else // Whole test is in else/endif block, because i/j would be unused on clang.
    int i = 32;
    constexpr int j = 42;

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
#endif
}
