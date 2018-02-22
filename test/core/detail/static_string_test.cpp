// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2017, Knut Reinert, FU Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
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

#include <gtest/gtest.h>

#include <seqan3/core/detail/static_string.hpp>

#include <string>

using namespace std::literals;
using namespace seqan3;

// standard construction.
TEST(static_string, standard_construction)
{
    EXPECT_FALSE((std::is_default_constructible_v<static_string<4>>));
    EXPECT_TRUE((std::is_copy_constructible_v<static_string<4>>));
    EXPECT_TRUE((std::is_trivially_copy_constructible_v<static_string<4>>));
    EXPECT_TRUE((std::is_nothrow_copy_constructible_v<static_string<4>>));
    EXPECT_TRUE((std::is_move_constructible_v<static_string<4>>));
    EXPECT_TRUE((std::is_trivially_move_constructible_v<static_string<4>>));
    EXPECT_TRUE((std::is_nothrow_move_constructible_v<static_string<4>>));
    EXPECT_TRUE((std::is_copy_assignable_v<static_string<4>>));
    EXPECT_TRUE((std::is_trivially_copy_assignable_v<static_string<4>>));
    EXPECT_TRUE((std::is_nothrow_copy_assignable_v<static_string<4>>));
    EXPECT_TRUE((std::is_move_assignable_v<static_string<4>>));
    EXPECT_TRUE((std::is_trivially_move_assignable_v<static_string<4>>));
    EXPECT_TRUE((std::is_nothrow_move_assignable_v<static_string<4>>));
}

// construction from literal.
TEST(static_string, construct_from_literal)
{
     EXPECT_TRUE((std::is_same_v<static_string<5>, decltype(static_string{"hello"})>));
}

// construction from char.
TEST(static_string, construct_from_char)
{
     EXPECT_TRUE((std::is_same_v<static_string<1>, decltype(static_string{'h'})>));
}

template <size_t S>
struct foo
{
    inline size_t get() const {return S;};
};

TEST(static_string, size)
{
    static_string em{"hello"};
    foo<em.size()> f;

    EXPECT_EQ(em.size(), 5u);
    EXPECT_EQ(f.get(), 5u);
}

TEST(static_string, c_str)
{
    {
        static_string em{"hello"};
        EXPECT_EQ(std::string{em.c_str()}, "hello"s);
    }

    {
        static_string em{'x'};
        EXPECT_EQ(std::string{em.c_str()}, "x"s);
    }
}

TEST(static_string, string)
{
    static_string em{"hello"};
    EXPECT_EQ(em.string(), "hello"s);
}

TEST(static_string, concat)
{
    {
        static_string em = static_string{"hello"} + static_string{' '} + static_string{"world"};
        foo<em.size()> f;
        EXPECT_EQ(f.get(), 11u);
        EXPECT_EQ(em.string(), "hello world"s);
    }

    {
        static constexpr char const a[] = "hello";
        static constexpr char const b[] = " ";
        static constexpr char const c[] = "world";
        auto em = static_string{a} + static_string{b} + static_string{c};
        EXPECT_EQ(em.string(), "hello world"s);
        foo<em.size()> f;
        EXPECT_EQ(f.get(), 11u);
    }
}
