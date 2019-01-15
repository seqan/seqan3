// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin 
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik 
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <iostream>
#include <memory>

#include <gtest/gtest.h>

#include <seqan3/range/view/detail.hpp>

using namespace seqan3;

/*************++ without args ************************/

// not actually a view
template <typename urng_t>
struct test_view
{
    urng_t const & urange;

    test_view(urng_t const & _urange) :
        urange{_urange}
    {}
};

TEST(generic_pipable_view_adaptor_without_args, function_interface)
{
    using test_fn = detail::generic_pipable_view_adaptor<test_view>;

    test_fn f;

    std::vector<int> urange{1, 2, 3};
    auto v = f(urange);

    EXPECT_TRUE((std::is_same_v<decltype(v), test_view<std::vector<int>>>));
    EXPECT_EQ(v.urange[0], urange[0]);
    EXPECT_EQ(v.urange[1], urange[1]);
    EXPECT_EQ(v.urange[2], urange[2]);
}

TEST(generic_pipable_view_adaptor_without_args, pipe_interface)
{
    using test_fn = detail::generic_pipable_view_adaptor<test_view>;

    test_fn f;

    std::vector<int> urange{1, 2, 3};
    auto v = urange | f;

    EXPECT_TRUE((std::is_same_v<decltype(v), test_view<std::vector<int>>>));
    EXPECT_EQ(v.urange[0], urange[0]);
    EXPECT_EQ(v.urange[1], urange[1]);
    EXPECT_EQ(v.urange[2], urange[2]);
}

/******************* with args ************************/

struct test_param
{
    uint64_t number;
};

template <typename urng_t>
struct test_view_args
{
    urng_t const & urange;
    test_param param;

    test_view_args(urng_t const & _urange, test_param const & _param) :
        urange{_urange}, param {_param}
    {}

    test_view_args(urng_t const & _urange, test_param && _param) :
        urange{_urange}, param{std::move(_param)}
    {}
};

TEST(generic_pipable_view_adaptor_with_args, function_interface)
{
    using test_fn = detail::generic_pipable_view_adaptor<test_view_args>;

    test_fn f;

    std::vector<int> urange{1, 2, 3};
    auto v = f(urange, test_param{7});

    EXPECT_TRUE((std::is_same_v<decltype(v), test_view_args<std::vector<int>>>));
    EXPECT_EQ(v.urange[0], urange[0]);
    EXPECT_EQ(v.urange[1], urange[1]);
    EXPECT_EQ(v.urange[2], urange[2]);
    EXPECT_EQ(v.param.number, 7u);
}

TEST(generic_pipable_view_adaptor_with_args, pipe_interface)
{
    using test_fn = detail::generic_pipable_view_adaptor<test_view_args>;

    test_fn f;

    std::vector<int> urange{1, 2, 3};
    auto v = urange | f(test_param{7});

    EXPECT_TRUE((std::is_same_v<decltype(v), test_view_args<std::vector<int>>>));
    EXPECT_EQ(v.urange[0], urange[0]);
    EXPECT_EQ(v.urange[1], urange[1]);
    EXPECT_EQ(v.urange[2], urange[2]);
    EXPECT_EQ(v.param.number, 7u);
}
