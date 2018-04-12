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

#include <gtest/gtest.h>

#include <functional>
#include <string>

#include <seqan3/core/detail/strong_type.hpp>

using namespace std::literals;
using namespace seqan3;

struct pure_type : detail::strong_type<int, pure_type>
{
    using detail::strong_type<int, pure_type>::strong_type;
};
struct additive_type : detail::strong_type<int, additive_type, detail::strong_type_skill::additive>
{
    using detail::strong_type<int, additive_type, detail::strong_type_skill::additive>::strong_type;
};

struct multiplicative_type : detail::strong_type<int, multiplicative_type, detail::strong_type_skill::multiplicative>
{
    using detail::strong_type<int, multiplicative_type, detail::strong_type_skill::multiplicative>::strong_type;
};

struct bitwise_type : detail::strong_type<unsigned, bitwise_type, detail::strong_type_skill::bitwise_logic>
{
    using detail::strong_type<unsigned, bitwise_type, detail::strong_type_skill::bitwise_logic>::strong_type;
};

struct bitwise_shift_type : detail::strong_type<unsigned, bitwise_shift_type, detail::strong_type_skill::bitwise_shift>
{
    using detail::strong_type<unsigned, bitwise_shift_type, detail::strong_type_skill::bitwise_shift>::strong_type;
};

struct logic_type : detail::strong_type<bool, logic_type, detail::strong_type_skill::logic>
{
    using detail::strong_type<bool, logic_type, detail::strong_type_skill::logic>::strong_type;
};

struct inc_type : detail::strong_type<int, inc_type, detail::strong_type_skill::increment>
{
    using detail::strong_type<int, inc_type, detail::strong_type_skill::increment>::strong_type;
};

struct dec_type : detail::strong_type<int, dec_type, detail::strong_type_skill::decrement>
{
    using detail::strong_type<int, dec_type, detail::strong_type_skill::decrement>::strong_type;
};

struct lval_type : detail::strong_type<std::reference_wrapper<std::string>, lval_type>
{
    using detail::strong_type<std::reference_wrapper<std::string>, lval_type>::strong_type;
};

struct convertible_type : detail::strong_type<int, convertible_type, detail::strong_type_skill::convert>
{
    using detail::strong_type<int, convertible_type, detail::strong_type_skill::convert>::strong_type;
};

struct multi_skill_type : detail::strong_type<int,
                                              multi_skill_type,
                                              detail::strong_type_skill::additive  |
                                              detail::strong_type_skill::increment |
                                              detail::strong_type_skill::decrement |
                                              detail::strong_type_skill::convert>
{
    using detail::strong_type<int, multi_skill_type, detail::strong_type_skill::additive  |
                                                     detail::strong_type_skill::increment |
                                                     detail::strong_type_skill::decrement |
                                                     detail::strong_type_skill::convert>::strong_type;
};

TEST(strong_type, pure_type)
{
    {
        EXPECT_TRUE(std::is_default_constructible_v<pure_type>);
        EXPECT_TRUE(std::is_copy_constructible_v<pure_type>);
        EXPECT_TRUE(std::is_move_constructible_v<pure_type>);
        EXPECT_TRUE(std::is_copy_assignable_v<pure_type>);
        EXPECT_TRUE(std::is_move_assignable_v<pure_type>);
        EXPECT_TRUE(std::is_destructible_v<pure_type>);
    }

    { // From l-value
        int val = 1;
        pure_type p{val};
        EXPECT_EQ(p.get(), 1);
    }

    {  // From r-value
        pure_type p{10};
        EXPECT_EQ(p.get(), 10);
    }
}

TEST(strong_type, get)
{
    pure_type p1{1};
    pure_type const p2{p1};

    { // strong_type &
        EXPECT_EQ(p1.get(), 1);
    }

    { // const strong_type &
        EXPECT_EQ(p2.get(), 1);
    }

    { // strong_type &&
        auto && v = std::move(p1).get();
        EXPECT_TRUE((std::is_same_v<decltype(v), int &&>));
        EXPECT_EQ(v, 1);
    }

    { // const strong_type &&
        auto && v = std::move(p2).get();
        EXPECT_TRUE((std::is_same_v<decltype(v), int const &&>));
        EXPECT_EQ(v, 1);
    }
}

TEST(strong_type, additive_type)
{
    additive_type f1{10};
    additive_type f2{10};

    additive_type f3 = f1 + f2;
    EXPECT_EQ(f3.get(), 20);

    f3 = f2 - f3;
    EXPECT_EQ(f3.get(), -10);
}

TEST(strong_type, multiplicative_type)
{
    multiplicative_type f1{10};
    multiplicative_type f2{10};

    multiplicative_type f3 = f1 * f2;
    EXPECT_EQ(f3.get(), 100);

    f3 = f3 / f1;
    EXPECT_EQ(f3.get(), 10);

    f3 = f3 % f1;
    EXPECT_EQ(f3.get(), 0);
}

TEST(strong_type, bitwise_logic_type)
{
    bitwise_type f1{1};
    bitwise_type f2{2};

    bitwise_type f3 = f1 | f2;
    EXPECT_EQ(f3.get(), 3u);

    f3 = f3 & f1;
    EXPECT_EQ(f3.get(), 1u);

    f3 = ~f3;
    EXPECT_EQ(f3.get(), std::numeric_limits<unsigned>::max() - 1u);

    f3 = f3 ^ f2;
    EXPECT_EQ(f3.get(), std::numeric_limits<unsigned>::max() - 3u);
}

TEST(strong_type, bitwise_shift_type)
{
    bitwise_shift_type f1{1};
    bitwise_shift_type f2{2};

    bitwise_shift_type f3 = f2 << f1;
    EXPECT_EQ(f3.get(), 4u);

    f3 = f3 << 1;
    EXPECT_EQ(f3.get(), 8u);

    f3 = f3 >> f1;
    EXPECT_EQ(f3.get(), 4u);

    f3 = f3 >> 1;
    EXPECT_EQ(f3.get(), 2u);
}

TEST(strong_type, logic_type)
{
    logic_type f1{true};
    logic_type f2{false};

    bool f3 = f1 || f2;
    EXPECT_EQ(f3, true);

    f3 = f1 && f2;
    EXPECT_EQ(f3, false);

    f3 = !f2;
    EXPECT_EQ(f3, true);
}

TEST(strong_type, increment_type)
{
    inc_type f1{10};
    inc_type & ref = ++f1;
    EXPECT_EQ(f1.get(), 11);
    EXPECT_EQ((f1++).get(), 11);
    EXPECT_EQ(f1.get(), 12);
    ++ref;
    EXPECT_EQ(f1.get(), 13);
}

TEST(strong_type, decrement_type)
{
    dec_type f1{10};
    dec_type & ref = --f1;
    EXPECT_EQ(f1.get(), 9);
    EXPECT_EQ((f1--).get(), 9);
    EXPECT_EQ(f1.get(), 8);
    --ref;
    EXPECT_EQ(f1.get(), 7);
}

TEST(strong_type, lval_type)
{
    std::string s{"test"};
    lval_type f1{std::ref(s)};

    EXPECT_EQ(f1.get().get(), "test"s);
}

TEST(strong_type, convertible_type)
{
    convertible_type f1{1};
    int v{f1};
    EXPECT_EQ(v, f1.get());
}

TEST(strong_type, multi_skill_type)
{
    multi_skill_type f1{1};
    multi_skill_type f2{1};
    int v(++f1 - f2--);
    EXPECT_EQ(v, 1);
}
