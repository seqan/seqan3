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

#include <seqan3/range/container/constexpr_string.hpp>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/range/concept.hpp>

#include <utility>
#include <string>

using namespace std::literals;
using namespace seqan3;

// standard construction.
TEST(constexpr_string, standard_construction)
{
    EXPECT_TRUE((std::is_default_constructible_v<constexpr_string<4>>));
    EXPECT_TRUE((std::is_trivially_default_constructible_v<constexpr_string<4>>));
    EXPECT_TRUE((std::is_nothrow_default_constructible_v<constexpr_string<4>>));
    EXPECT_TRUE((std::is_copy_constructible_v<constexpr_string<4>>));
    EXPECT_TRUE((std::is_trivially_copy_constructible_v<constexpr_string<4>>));
    EXPECT_TRUE((std::is_nothrow_copy_constructible_v<constexpr_string<4>>));
    EXPECT_TRUE((std::is_move_constructible_v<constexpr_string<4>>));
    EXPECT_TRUE((std::is_trivially_move_constructible_v<constexpr_string<4>>));
    EXPECT_TRUE((std::is_nothrow_move_constructible_v<constexpr_string<4>>));
    EXPECT_TRUE((std::is_copy_assignable_v<constexpr_string<4>>));
    EXPECT_TRUE((std::is_trivially_copy_assignable_v<constexpr_string<4>>));
    EXPECT_TRUE((std::is_nothrow_copy_assignable_v<constexpr_string<4>>));
    EXPECT_TRUE((std::is_move_assignable_v<constexpr_string<4>>));
    EXPECT_TRUE((std::is_trivially_move_assignable_v<constexpr_string<4>>));
    EXPECT_TRUE((std::is_nothrow_move_assignable_v<constexpr_string<4>>));
}

TEST(constexpr_string, container_concept)
{
    EXPECT_TRUE(container_concept<constexpr_string<4>>);
    EXPECT_TRUE(random_access_range_concept<constexpr_string<4>>);
}

// construction from literal.
TEST(constexpr_string, construct_from_literal)
{
     EXPECT_TRUE((std::is_same_v<constexpr_string<5>, decltype(constexpr_string{"hello"})>));
}

// construction from char.
TEST(constexpr_string, construct_from_char)
{
     EXPECT_TRUE((std::is_same_v<constexpr_string<1>, decltype(constexpr_string{'h'})>));
}

// construction from array.
TEST(constexpr_string, construct_from_array)
{
     EXPECT_TRUE((std::is_same_v<constexpr_string<5>,
                                 decltype(constexpr_string{std::array{'h','e','l','l','o'}})>));
}

template <size_t S>
struct foo
{
    inline size_t get() const {return S;};
};

TEST(constexpr_string, size)
{
    constexpr_string em{"hello"};
    foo<em.size()> f;

    EXPECT_EQ(em.size(), 5u);
    EXPECT_EQ(f.get(), 5u);
}

TEST(constexpr_string, max_size)
{
    constexpr_string em{"hello"};
    foo<em.max_size()> f;

    EXPECT_EQ(em.max_size(), 5u);
    EXPECT_EQ(f.get(), 5u);
}

TEST(constexpr_string, c_str)
{
    {
        constexpr_string em{"hello"};
        EXPECT_EQ(std::string{em.c_str()}, "hello"s);
    }

    {
        constexpr_string em{'x'};
        EXPECT_EQ(std::string{em.c_str()}, "x"s);
    }
}

TEST(constexpr_string, string)
{
    constexpr_string em{"hello"};
    EXPECT_EQ(em.string(), "hello"s);  // explicit
    EXPECT_EQ(std::string{em}, "hello"s);  // implicit
}

TEST(constexpr_string, concat)
{
    {
        constexpr_string em = constexpr_string{"hello"} +
                                   constexpr_string{' '} +
                                   constexpr_string{"world"};
        foo<em.size()> f;
        EXPECT_EQ(f.get(), 11u);
        EXPECT_EQ(em.string(), "hello world"s);
    }

    {
        static constexpr char const a[] = "hello";
        static constexpr char const b[] = " ";
        static constexpr char const c[] = "world";
        auto em = constexpr_string{a} + constexpr_string{b} + constexpr_string{c};
        EXPECT_EQ(em.string(), "hello world"s);
        foo<em.size()> f;
        EXPECT_EQ(f.get(), 11u);
    }
}

TEST(constexpr_string, begin)
{
    constexpr_string s{"hello"};
    EXPECT_EQ(*s.begin(), 'h');
    EXPECT_TRUE((std::is_same_v<decltype(s)::iterator, decltype(s.begin())>));

    constexpr_string<5> const cs{s};
    EXPECT_EQ(*cs.begin(), 'h');
    EXPECT_TRUE((std::is_same_v<decltype(cs)::const_iterator, decltype(cs.begin())>));
}

TEST(constexpr_string, cbegin)
{
    constexpr_string s{"hello"};
    EXPECT_EQ(*s.begin(), 'h');
    EXPECT_TRUE((std::is_same_v<decltype(s)::const_iterator, decltype(s.cbegin())>));
}

TEST(constexpr_string, end)
{
    constexpr_string s{"hello"};
    EXPECT_EQ(*(s.end() - 1), 'o');
    EXPECT_TRUE((std::is_same_v<decltype(s)::iterator, decltype(s.end())>));

    constexpr_string<5> const cs{s};
    static_assert(std::is_same_v<constexpr_string<5> const, decltype(cs)>);
    EXPECT_EQ(*(cs.end() - 1), 'o');
    EXPECT_TRUE((std::is_same_v<decltype(cs)::const_iterator, decltype(cs.end())>));
}

TEST(constexpr_string, cend)
{
    constexpr_string s{"hello"};
    EXPECT_EQ(*(s.cend() - 1), 'o');
    EXPECT_TRUE((std::is_same_v<decltype(s)::const_iterator, decltype(s.cend())>));
}

TEST(constexpr_string, swap)
{
    constexpr_string s1{"hello"};
    constexpr_string s2{"olleh"};
    {  // global function.
        swap(s1, s2);
        EXPECT_EQ(s1, constexpr_string{"olleh"});
        EXPECT_EQ(s2, constexpr_string{"hello"});
    }

    { // member function
        s1.swap(s2);
        EXPECT_EQ(s1, constexpr_string{"hello"});
        EXPECT_EQ(s2, constexpr_string{"olleh"});
    }
}

template <bool B>
struct bar
{
    bool get(){
        return B;
    }
};

TEST(constexpr_string, equality)
{
    EXPECT_TRUE(bar<(constexpr_string{"hello"}  == constexpr_string{"hello"})>{}.get());
    EXPECT_FALSE(bar<(constexpr_string{"hello"} == constexpr_string{"hell"})>{}.get());
    EXPECT_FALSE(bar<(constexpr_string{"hell"}  == constexpr_string{"hello"})>{}.get());
    EXPECT_FALSE(bar<(constexpr_string{"hella"} == constexpr_string{"hello"})>{}.get());
}

TEST(constexpr_string, inequality)
{
    EXPECT_FALSE(bar<(constexpr_string{"hello"} != constexpr_string{"hello"})>{}.get());
    EXPECT_TRUE(bar<(constexpr_string{"hello"}  != constexpr_string{"hell"})>{}.get());
    EXPECT_TRUE(bar<(constexpr_string{"hell"}   != constexpr_string{"hello"})>{}.get());
    EXPECT_TRUE(bar<(constexpr_string{"hella"}  != constexpr_string{"hello"})>{}.get());
}

template <std::size_t N>
constexpr constexpr_string<N> fill_constexpr_string(constexpr_string<N> s, char const val)
{
    for (auto && c : s)
    {
        c = val;
    }
    return s;
}

TEST(constexpr_string, compile_time_fill)
{
    EXPECT_TRUE(bar<(fill_constexpr_string(constexpr_string<4>{}, 'x') == constexpr_string{"xxxx"})>{}.get());
}
