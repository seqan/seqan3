// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <functional>
#include <string>

#include <seqan3/core/detail/strong_type.hpp>

using namespace std::literals;

namespace seqan3::detail
{
// Make operator| accessible by ADL for seqan3 namespace
using seqan3::operator|;
} // namespace seqan3::detail

struct pure_type : seqan3::detail::strong_type<int, pure_type>
{
    using seqan3::detail::strong_type<int, pure_type>::strong_type;
};
struct additive_type : seqan3::detail::strong_type<int, additive_type, seqan3::detail::strong_type_skill::additive>
{
    using seqan3::detail::strong_type<int, additive_type, seqan3::detail::strong_type_skill::additive>::strong_type;
};

struct multiplicative_type :
    seqan3::detail::strong_type<int, multiplicative_type, seqan3::detail::strong_type_skill::multiplicative>
{
    using seqan3::detail::strong_type<int, multiplicative_type, seqan3::detail::strong_type_skill::multiplicative>::
        strong_type;
};

struct bitwise_type :
    seqan3::detail::strong_type<unsigned, bitwise_type, seqan3::detail::strong_type_skill::bitwise_logic>
{
    using seqan3::detail::strong_type<unsigned, bitwise_type, seqan3::detail::strong_type_skill::bitwise_logic>::
        strong_type;
};

struct bitwise_shift_type :
    seqan3::detail::strong_type<unsigned, bitwise_shift_type, seqan3::detail::strong_type_skill::bitwise_shift>
{
    using seqan3::detail::strong_type<unsigned, bitwise_shift_type, seqan3::detail::strong_type_skill::bitwise_shift>::
        strong_type;
};

struct logic_type : seqan3::detail::strong_type<bool, logic_type, seqan3::detail::strong_type_skill::logic>
{
    using seqan3::detail::strong_type<bool, logic_type, seqan3::detail::strong_type_skill::logic>::strong_type;
};

struct inc_type : seqan3::detail::strong_type<int, inc_type, seqan3::detail::strong_type_skill::increment>
{
    using seqan3::detail::strong_type<int, inc_type, seqan3::detail::strong_type_skill::increment>::strong_type;
};

struct dec_type : seqan3::detail::strong_type<int, dec_type, seqan3::detail::strong_type_skill::decrement>
{
    using seqan3::detail::strong_type<int, dec_type, seqan3::detail::strong_type_skill::decrement>::strong_type;
};

struct lval_type : seqan3::detail::strong_type<std::reference_wrapper<std::string>, lval_type>
{
    using seqan3::detail::strong_type<std::reference_wrapper<std::string>, lval_type>::strong_type;
};

struct convertible_type : seqan3::detail::strong_type<int, convertible_type, seqan3::detail::strong_type_skill::convert>
{
    using seqan3::detail::strong_type<int, convertible_type, seqan3::detail::strong_type_skill::convert>::strong_type;
};

struct comp_type : seqan3::detail::strong_type<int, comp_type, seqan3::detail::strong_type_skill::comparable>
{
    using seqan3::detail::strong_type<int, comp_type, seqan3::detail::strong_type_skill::comparable>::strong_type;
};

struct multi_skill_type :
    seqan3::detail::strong_type<
        int,
        multi_skill_type,
        seqan3::detail::strong_type_skill::additive | seqan3::detail::strong_type_skill::increment
            | seqan3::detail::strong_type_skill::decrement | seqan3::detail::strong_type_skill::convert>
{
    using seqan3::detail::strong_type<
        int,
        multi_skill_type,
        seqan3::detail::strong_type_skill::additive | seqan3::detail::strong_type_skill::increment
            | seqan3::detail::strong_type_skill::decrement | seqan3::detail::strong_type_skill::convert>::strong_type;
};

TEST(strong_type, concept)
{
    EXPECT_TRUE(seqan3::detail::derived_from_strong_type<pure_type &>);
    EXPECT_TRUE(seqan3::detail::derived_from_strong_type<pure_type const &>);
    EXPECT_TRUE(seqan3::detail::derived_from_strong_type<pure_type &&>);
    EXPECT_TRUE(seqan3::detail::derived_from_strong_type<pure_type const &&>);
}

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

    { // From r-value
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

TEST(strong_type, comparable_type)
{
    comp_type f1{1};
    comp_type f2{1};
    comp_type f3{42};

    EXPECT_EQ(f1, f2);
    EXPECT_NE(f1, f3);
    EXPECT_NE(f2, f3);
}

TEST(strong_type, multi_skill_type)
{
    multi_skill_type f1{1};
    multi_skill_type f2{1};
    int v(++f1 - f2--);
    EXPECT_EQ(v, 1);
}
