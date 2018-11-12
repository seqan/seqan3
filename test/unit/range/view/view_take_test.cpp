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

#include <seqan3/range/view/take.hpp>
#include <seqan3/range/view/take_exactly.hpp>
#include <seqan3/range/view/single_pass_input.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/view/reverse.hpp>
#include <seqan3/std/view/common.hpp>

using namespace seqan3;

// ============================================================================
//  test templates
// ============================================================================

template <typename adaptor_t>
void do_test(adaptor_t const & adaptor, std::string const & vec)
{
    // pipe notation
    auto v = vec | adaptor(3);
    EXPECT_EQ("foo", std::string(v));

    // function notation
    std::string v2 = adaptor(vec, 3);
    EXPECT_EQ("foo", v2);

    // combinability
    auto v3 = vec | adaptor(3) | adaptor(3) | ranges::view::unique;
    EXPECT_EQ("fo", std::string(v3));
    std::string v3b = vec | view::reverse | adaptor(3) | ranges::view::unique;
    EXPECT_EQ("rab", v3b);
}

template <typename adaptor_t>
void do_concepts(adaptor_t && adaptor, bool const exactly)
{
    std::string vec{"foobar"};
    EXPECT_TRUE(std::ranges::InputRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::ForwardRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::BidirectionalRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::RandomAccessRange<decltype(vec)>);
    EXPECT_FALSE(std::ranges::View<decltype(vec)>);
    EXPECT_TRUE(std::ranges::SizedRange<decltype(vec)>);
    EXPECT_TRUE(std::ranges::CommonRange<decltype(vec)>);
    EXPECT_TRUE(ConstIterableRange<decltype(vec)>);
    EXPECT_TRUE((std::ranges::OutputRange<decltype(vec), char>));

    auto v1 = vec | adaptor;

    EXPECT_TRUE(std::ranges::InputRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::ForwardRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::BidirectionalRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::RandomAccessRange<decltype(v1)>);
    EXPECT_TRUE(std::ranges::View<decltype(v1)>);
    EXPECT_EQ(std::ranges::SizedRange<decltype(v1)>, exactly);
    EXPECT_FALSE(std::ranges::CommonRange<decltype(v1)>);
    EXPECT_TRUE(ConstIterableRange<decltype(v1)>);
    EXPECT_TRUE((std::ranges::OutputRange<decltype(v1), char>));

    auto v2 = vec | view::single_pass_input | adaptor;

    EXPECT_TRUE(std::ranges::InputRange<decltype(v2)>);
    EXPECT_FALSE(std::ranges::ForwardRange<decltype(v2)>);
    EXPECT_FALSE(std::ranges::BidirectionalRange<decltype(v2)>);
    EXPECT_FALSE(std::ranges::RandomAccessRange<decltype(v2)>);
    EXPECT_TRUE(std::ranges::View<decltype(v2)>);
    EXPECT_EQ(std::ranges::SizedRange<decltype(v2)>, exactly);
    EXPECT_FALSE(std::ranges::CommonRange<decltype(v2)>);
    EXPECT_FALSE(ConstIterableRange<decltype(v2)>);
    EXPECT_TRUE((std::ranges::OutputRange<decltype(v2), char>));
}

// ============================================================================
//  view_take
// ============================================================================

TEST(view_take, regular)
{
    do_test(view::take, "foobar");
}

TEST(view_take, concepts)
{
    do_concepts(view::take(3), false);
}

TEST(view_take, underlying_is_shorter)
{
    std::string vec{"foo"};
    EXPECT_NO_THROW(( view::take(vec, 4) )); // no parsing

    std::string v;
    EXPECT_NO_THROW(( v = vec | view::single_pass_input | view::take(4) )); // full parsing on conversion
    EXPECT_EQ("foo", v);
}

// ============================================================================
//  view_take_exactly
// ============================================================================

TEST(view_take_exactly, regular)
{
    do_test(view::take_exactly, "foobar");
}

TEST(view_take_exactly, concepts)
{
    do_concepts(view::take_exactly(3), true);
}

TEST(view_take_exactly, underlying_is_shorter)
{
    std::string vec{"foo"};
    EXPECT_NO_THROW(( view::take_exactly(vec, 4) )); // no parsing

    std::string v;
    EXPECT_NO_THROW(( v = vec | view::single_pass_input | view::take_exactly(4) )); // full parsing on conversion
    EXPECT_EQ("foo", v);

    auto v2 = vec | view::single_pass_input | view::take_exactly(4);
    EXPECT_EQ(size(v2), 4); // here be dragons
}

// ============================================================================
//  view_take_exactly_or_throw
// ============================================================================

TEST(view_take_exactly_or_throw, regular)
{
    do_test(view::take_exactly_or_throw, "foo\nbar");
}

TEST(view_take_exactly_or_throw, concepts)
{
    do_concepts(view::take_exactly_or_throw(3), true);
}

TEST(view_take_exactly_or_throw, underlying_is_shorter)
{
    std::string vec{"foo"};
    EXPECT_THROW(( view::take_exactly_or_throw(vec, 4) ),
                   std::invalid_argument); // no parsing, but throws on construction

    std::string v;
    EXPECT_THROW(( v = vec | view::single_pass_input | view::take_exactly_or_throw(4)),
                   unexpected_end_of_input); // full parsing on conversion, throw on conversion
}

