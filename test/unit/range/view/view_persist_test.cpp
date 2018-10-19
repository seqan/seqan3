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

#include <gtest/gtest.h>

#include <range/v3/algorithm/copy.hpp>
#include <range/v3/view/unique.hpp>

#include <seqan3/io/stream/parse_condition.hpp>
#include <seqan3/range/view/persist.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/view/filter.hpp>
#include <seqan3/std/view/reverse.hpp>
#include <seqan3/std/view/common.hpp>

using namespace seqan3;

// ============================================================================
//  test templates
// ============================================================================

TEST(view_persist, delegate_to_view_all)
{
    std::string vec{"foo"};

    // pipe notation
    auto v = vec | view::persist;
    EXPECT_EQ("foo", std::string{v});

    // function notation
    std::string v2 = view::persist(vec);
    EXPECT_EQ("foo", v2);

    // combinability
    auto v3 = vec | view::persist | ranges::view::unique;
    EXPECT_EQ("fo", std::string{v3});
    std::string v3b = vec | view::reverse | view::persist | ranges::view::unique;
    EXPECT_EQ("of", v3b);
}

TEST(view_persist, wrap_temporary)
{
    // pipe notation
    auto v = std::string{"foo"} | view::persist;
    EXPECT_EQ("foo", std::string(v));

    // function notation
    std::string v2 = view::persist(std::string{"foo"});
    EXPECT_EQ("foo", v2);

    // combinability
    auto v3 = std::string{"foo"} | view::persist | ranges::view::unique;
    EXPECT_EQ("fo", std::string(v3));
    std::string v3b = std::string{"foo"} | view::persist | view::filter(is_char<'o'>) | ranges::view::unique;
    EXPECT_EQ("o", v3b);
}

TEST(view_persist, const_)
{
    // inner const
    using t = std::string const;
    auto v = t{"foo"} | view::persist;
    EXPECT_EQ("foo", std::string{v});

    // outer const
    auto const & v2 = std::string{"foo"} | view::persist;
    EXPECT_EQ("foo", std::string{v2});

    // inner + outer const
    using t = std::string const;
    auto const & v3 = t{"foo"} | view::persist;
    EXPECT_EQ("foo", std::string{v3});
}

TEST(view_persist, concepts)
{
    std::string vec{"foobar"};
    EXPECT_TRUE(std::ranges::InputRange<decltype(std::string{"foo"})>);
    EXPECT_TRUE(std::ranges::ForwardRange<decltype(std::string{"foo"})>);
    EXPECT_TRUE(std::ranges::BidirectionalRange<decltype(std::string{"foo"})>);
    EXPECT_TRUE(std::ranges::RandomAccessRange<decltype(std::string{"foo"})>);
    EXPECT_FALSE(std::ranges::View<decltype(std::string{"foo"})>);
    EXPECT_TRUE(std::ranges::SizedRange<decltype(std::string{"foo"})>);
    EXPECT_TRUE(std::ranges::CommonRange<decltype(std::string{"foo"})>);
    EXPECT_TRUE(const_iterable_concept<decltype(std::string{"foo"})>);
    EXPECT_TRUE((std::ranges::OutputRange<decltype(std::string{"foo"}), char>));

    auto v1 = std::string{"foo"} | view::persist;

    EXPECT_TRUE(std::ranges::InputRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::ForwardRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::BidirectionalRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::RandomAccessRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::View<decltype(v1)>);
    EXPECT_TRUE(std::ranges::SizedRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::CommonRange<decltype(v1)>);
    EXPECT_TRUE(const_iterable_concept<decltype(v1)>);
    EXPECT_TRUE((std::ranges::OutputRange<decltype(v1), char>));
}
