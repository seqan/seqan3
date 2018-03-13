// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

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

TEST(declare_view_functor_type_without_args, function_interface)
{
    using test_fn = detail::declare_view_functor_type<test_view>;

    test_fn f;

    std::vector<int> urange{1, 2, 3};
    auto v = f(urange);

    EXPECT_TRUE((std::is_same_v<decltype(v), test_view<std::vector<int>>>));
    EXPECT_EQ(v.urange[0], urange[0]);
    EXPECT_EQ(v.urange[1], urange[1]);
    EXPECT_EQ(v.urange[2], urange[2]);
}

TEST(declare_view_functor_type_without_args, pipe_interface)
{
    using test_fn = detail::declare_view_functor_type<test_view>;

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

TEST(declare_view_functor_type_with_args, function_interface)
{
    using test_fn = detail::declare_view_functor_type<test_view_args>;

    test_fn f;

    std::vector<int> urange{1, 2, 3};
    auto v = f(urange, test_param{7});

    EXPECT_TRUE((std::is_same_v<decltype(v), test_view_args<std::vector<int>>>));
    EXPECT_EQ(v.urange[0], urange[0]);
    EXPECT_EQ(v.urange[1], urange[1]);
    EXPECT_EQ(v.urange[2], urange[2]);
    EXPECT_EQ(v.param.number, 7u);
}

TEST(declare_view_functor_type_with_args, pipe_interface)
{
    using test_fn = detail::declare_view_functor_type<test_view_args>;

    test_fn f;

    std::vector<int> urange{1, 2, 3};
    auto v = urange | f(test_param{7});

    EXPECT_TRUE((std::is_same_v<decltype(v), test_view_args<std::vector<int>>>));
    EXPECT_EQ(v.urange[0], urange[0]);
    EXPECT_EQ(v.urange[1], urange[1]);
    EXPECT_EQ(v.urange[2], urange[2]);
    EXPECT_EQ(v.param.number, 7u);
}
