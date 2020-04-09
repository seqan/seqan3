// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/core/concept/function.hpp>

TEST(is_function, regular_function)
{
    using fn_void_0_param_t = void();
    using fn_bool_2_param_t = bool(int, double);

    using fn_void_0_param_const_t = void() const ;
    using fn_bool_2_param_const_t = bool(int, double) const;

    using fn_void_0_param_noexcept_t = void() noexcept;
    using fn_bool_2_param_noexcept_t = bool(int, double) noexcept;

    using fn_void_0_param_lvalue_ref_t = void() &;
    using fn_bool_2_param_lvalue_ref_t = bool(int, double) &;

    using fn_void_0_param_rvalue_ref_t = void() &&;
    using fn_bool_2_param_rvalue_ref_t = bool(int, double) &&;

    using fn_void_0_param_complex_t = void() const volatile && noexcept;
    using fn_bool_2_param_complex_t = bool(int, double) const volatile & noexcept;

    EXPECT_TRUE(seqan3::function_like<fn_void_0_param_t>);
    EXPECT_TRUE(seqan3::function_like<fn_bool_2_param_t>);
    EXPECT_TRUE(seqan3::function_like<fn_void_0_param_const_t>);
    EXPECT_TRUE(seqan3::function_like<fn_bool_2_param_const_t>);
    EXPECT_TRUE(seqan3::function_like<fn_void_0_param_noexcept_t>);
    EXPECT_TRUE(seqan3::function_like<fn_bool_2_param_noexcept_t>);
    EXPECT_TRUE(seqan3::function_like<fn_void_0_param_lvalue_ref_t>);
    EXPECT_TRUE(seqan3::function_like<fn_bool_2_param_lvalue_ref_t>);
    EXPECT_TRUE(seqan3::function_like<fn_void_0_param_rvalue_ref_t>);
    EXPECT_TRUE(seqan3::function_like<fn_bool_2_param_rvalue_ref_t>);
    EXPECT_TRUE(seqan3::function_like<fn_void_0_param_complex_t>);
    EXPECT_TRUE(seqan3::function_like<fn_bool_2_param_complex_t>);
}

TEST(is_function, non_member_function_ptr)
{
    using fn_void_0_param_t = void(*)();
    using fn_bool_2_param_t = bool(*)(int, double);

    using fn_bool_2_param_const_t = bool(int const, double const);

    using fn_void_0_param_noexcept_t = void(*)() noexcept;
    using fn_bool_2_param_noexcept_t = bool(*)(int, double) noexcept;

    using fn_bool_2_param_lvalue_ref_t = bool(*)(int &, double &);
    using fn_bool_2_param_rvalue_ref_t = bool(int &&, double &&);
    using fn_bool_2_param_complex_t = bool(int const volatile &, double const volatile &&) noexcept;

    EXPECT_TRUE(seqan3::function_like<fn_void_0_param_t>);
    EXPECT_TRUE(seqan3::function_like<fn_bool_2_param_t>);
    EXPECT_TRUE(seqan3::function_like<fn_bool_2_param_const_t>);
    EXPECT_TRUE(seqan3::function_like<fn_void_0_param_noexcept_t>);
    EXPECT_TRUE(seqan3::function_like<fn_bool_2_param_noexcept_t>);
    EXPECT_TRUE(seqan3::function_like<fn_bool_2_param_lvalue_ref_t>);
    EXPECT_TRUE(seqan3::function_like<fn_bool_2_param_rvalue_ref_t>);
    EXPECT_TRUE(seqan3::function_like<fn_bool_2_param_complex_t>);
}

TEST(is_function, std_function)
{
    using fn_void_0_param_t = std::function<void()>;
    using fn_bool_2_param_t = std::function<bool(int, double)>;

    using fn_void_0_param_const_t = std::function<void()> const;
    using fn_bool_2_param_const_t = std::function<bool(int, double)> const;

    using fn_void_0_param_lvalue_ref_t = std::function<void()> &;
    using fn_bool_2_param_lvalue_ref_t = std::function<bool(int, double)> &;

    using fn_void_0_param_rvalue_ref_t = std::function<void()> &&;
    using fn_bool_2_param_rvalue_ref_t = std::function<bool(int, double)> &&;

    using fn_void_0_param_complex_t = std::function<void()> const &&;
    using fn_bool_2_param_complex_t = std::function<bool(int, double)> const &;

    EXPECT_TRUE(seqan3::function_like<fn_void_0_param_t>);
    EXPECT_TRUE(seqan3::function_like<fn_bool_2_param_t>);
    EXPECT_TRUE(seqan3::function_like<fn_void_0_param_const_t>);
    EXPECT_TRUE(seqan3::function_like<fn_bool_2_param_const_t>);
    EXPECT_TRUE(seqan3::function_like<fn_void_0_param_lvalue_ref_t>);
    EXPECT_TRUE(seqan3::function_like<fn_bool_2_param_lvalue_ref_t>);
    EXPECT_TRUE(seqan3::function_like<fn_void_0_param_rvalue_ref_t>);
    EXPECT_TRUE(seqan3::function_like<fn_bool_2_param_rvalue_ref_t>);
    EXPECT_TRUE(seqan3::function_like<fn_void_0_param_complex_t>);
    EXPECT_TRUE(seqan3::function_like<fn_bool_2_param_complex_t>);
}

TEST(is_function, captureless_lambda)
{
    auto fn_void_0_param = [] () {};
    auto fn_bool_2_param = [] (int, double) { return true; };

    auto fn_void_0_param_const = [] () constexpr {};
    auto fn_bool_2_param_const = [] (int, double) constexpr { return true; };

    auto fn_void_0_param_noexcept = [] () noexcept {};
    auto fn_bool_2_param_noexcept = [] (int, double) noexcept { return true; };

    auto fn_void_0_param_complex = [] () constexpr noexcept {};
    auto fn_bool_2_param_complex = [] (int, double) constexpr noexcept { return true; };

    using fn_void_0_param_t = decltype(fn_void_0_param);
    using fn_bool_2_param_t = decltype(fn_bool_2_param);

    using fn_void_0_param_const_t = decltype(fn_void_0_param_const);
    using fn_bool_2_param_const_t = decltype(fn_bool_2_param_const);

    using fn_void_0_param_noexcept_t = decltype(fn_void_0_param_noexcept);
    using fn_bool_2_param_noexcept_t = decltype(fn_bool_2_param_noexcept);

    using fn_void_0_param_lvalue_ref_t = decltype(fn_void_0_param) &;
    using fn_bool_2_param_lvalue_ref_t = decltype(fn_bool_2_param) &;

    using fn_void_0_param_rvalue_ref_t = decltype(fn_void_0_param) &&;
    using fn_bool_2_param_rvalue_ref_t = decltype(fn_bool_2_param) &&;

    using fn_void_0_param_complex_t = decltype(fn_void_0_param_complex) &;
    using fn_bool_2_param_complex_t = decltype(fn_bool_2_param_complex) &&;

    EXPECT_TRUE(seqan3::function_like<fn_void_0_param_t>);
    EXPECT_TRUE(seqan3::function_like<fn_bool_2_param_t>);
    EXPECT_TRUE(seqan3::function_like<fn_void_0_param_const_t>);
    EXPECT_TRUE(seqan3::function_like<fn_bool_2_param_const_t>);
    EXPECT_TRUE(seqan3::function_like<fn_void_0_param_noexcept_t>);
    EXPECT_TRUE(seqan3::function_like<fn_bool_2_param_noexcept_t>);
    EXPECT_TRUE(seqan3::function_like<fn_void_0_param_lvalue_ref_t>);
    EXPECT_TRUE(seqan3::function_like<fn_bool_2_param_lvalue_ref_t>);
    EXPECT_TRUE(seqan3::function_like<fn_void_0_param_rvalue_ref_t>);
    EXPECT_TRUE(seqan3::function_like<fn_bool_2_param_rvalue_ref_t>);
    EXPECT_TRUE(seqan3::function_like<fn_void_0_param_complex_t>);
    EXPECT_TRUE(seqan3::function_like<fn_bool_2_param_complex_t>);
}

TEST(is_function, capture_lambda)
{
    bool captured_variable = false;
    auto fn_void_0_param = [&] () mutable { captured_variable = true; };
    auto fn_bool_2_param = [=] (int, double) { return true && captured_variable; };

    auto fn_void_0_param_const = [&] () mutable constexpr { captured_variable = true; };
    auto fn_bool_2_param_const = [=] (int, double) constexpr { return true && captured_variable; };

    auto fn_void_0_param_noexcept = [&] () mutable noexcept { captured_variable = true; };
    auto fn_bool_2_param_noexcept = [=] (int, double) noexcept { return true; };

    auto fn_void_0_param_complex = [&] () mutable constexpr noexcept { captured_variable = true; };
    auto fn_bool_2_param_complex = [=] (int, double) constexpr noexcept { return true && captured_variable; };

    using fn_void_0_param_t = decltype(fn_void_0_param);
    using fn_bool_2_param_t = decltype(fn_bool_2_param);

    using fn_void_0_param_const_t = decltype(fn_void_0_param_const);
    using fn_bool_2_param_const_t = decltype(fn_bool_2_param_const);

    using fn_void_0_param_noexcept_t = decltype(fn_void_0_param_noexcept);
    using fn_bool_2_param_noexcept_t = decltype(fn_bool_2_param_noexcept);

    using fn_void_0_param_lvalue_ref_t = decltype(fn_void_0_param) &;
    using fn_bool_2_param_lvalue_ref_t = decltype(fn_bool_2_param) &;

    using fn_void_0_param_rvalue_ref_t = decltype(fn_void_0_param) &&;
    using fn_bool_2_param_rvalue_ref_t = decltype(fn_bool_2_param) &&;

    using fn_void_0_param_complex_t = decltype(fn_void_0_param_complex) &;
    using fn_bool_2_param_complex_t = decltype(fn_bool_2_param_complex) &&;

    EXPECT_TRUE(seqan3::function_like<fn_void_0_param_t>);
    EXPECT_TRUE(seqan3::function_like<fn_bool_2_param_t>);
    EXPECT_TRUE(seqan3::function_like<fn_void_0_param_const_t>);
    EXPECT_TRUE(seqan3::function_like<fn_bool_2_param_const_t>);
    EXPECT_TRUE(seqan3::function_like<fn_void_0_param_noexcept_t>);
    EXPECT_TRUE(seqan3::function_like<fn_bool_2_param_noexcept_t>);
    EXPECT_TRUE(seqan3::function_like<fn_void_0_param_lvalue_ref_t>);
    EXPECT_TRUE(seqan3::function_like<fn_bool_2_param_lvalue_ref_t>);
    EXPECT_TRUE(seqan3::function_like<fn_void_0_param_rvalue_ref_t>);
    EXPECT_TRUE(seqan3::function_like<fn_bool_2_param_rvalue_ref_t>);
    EXPECT_TRUE(seqan3::function_like<fn_void_0_param_complex_t>);
    EXPECT_TRUE(seqan3::function_like<fn_bool_2_param_complex_t>);
}

struct fn_object_simple
{
    auto operator()(int)
    {
        return true;
    }
};

struct fn_object_const_lvalue_ref
{
    auto operator()(int) &
    {
        return true;
    }
};

struct fn_object_const_complex
{
    auto operator()(int) const & noexcept(std::is_nothrow_default_constructible_v<bool>)
    {
        return true;
    }
};

TEST(is_function, function_object)
{
    EXPECT_TRUE(seqan3::function_like<fn_object_simple>);
    EXPECT_TRUE(seqan3::function_like<fn_object_const_lvalue_ref>);
    EXPECT_TRUE(seqan3::function_like<fn_object_const_complex>);
    EXPECT_TRUE(seqan3::function_like<fn_object_simple &>);
    EXPECT_TRUE(seqan3::function_like<fn_object_const_lvalue_ref &>);
    EXPECT_TRUE(seqan3::function_like<fn_object_const_complex &>);
    EXPECT_TRUE(seqan3::function_like<fn_object_simple &&>);
    EXPECT_TRUE(seqan3::function_like<fn_object_const_lvalue_ref &&>);
    EXPECT_TRUE(seqan3::function_like<fn_object_const_complex &&>);
}

struct fn_object_rvalue_ref
{
    auto operator()(int) &&
    {
        return true;
    }
};

TEST(is_function, false_functions)
{
    auto generic_fn = [] (auto x) { return x; };
    EXPECT_FALSE(seqan3::function_like<int>);
    EXPECT_FALSE(seqan3::function_like<int const *>);
    EXPECT_FALSE(seqan3::function_like<fn_object_rvalue_ref>);
    EXPECT_FALSE(seqan3::function_like<decltype(generic_fn)>);
}
