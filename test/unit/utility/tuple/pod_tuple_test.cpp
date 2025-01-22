// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/utility/tuple/pod_tuple.hpp>

using tuple_t = seqan3::pod_tuple<int, long, float>;

TEST(pod_tuple, concepts)
{
    EXPECT_TRUE(seqan3::trivial<tuple_t>);
    EXPECT_TRUE(seqan3::standard_layout<tuple_t>);
}

// default/zero construction
TEST(pod_tuple, ctr)
{
    [[maybe_unused]] tuple_t t1;
}

// aggregate initialization
TEST(pod_tuple, aggr)
{
    [[maybe_unused]] tuple_t t1{4, 7l, 3.0f};
    [[maybe_unused]] tuple_t t2(4, 7l, 3.0f);
}

// zero initialization
TEST(pod_tuple, zro)
{
    tuple_t t1{0, 0, 0};
    tuple_t t2{};

    EXPECT_EQ(t1, t2);
}

// copy construction
TEST(pod_tuple, cp_ctr)
{
    tuple_t t1{4, 7l, 3.0f};
    tuple_t t2{t1};
    tuple_t t3(t1);
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

// move construction
TEST(pod_tuple, mv_ctr)
{
    tuple_t t0{4, 7l, 3.0f};
    tuple_t t1{4, 7l, 3.0f};
    tuple_t t2{std::move(t1)};
    EXPECT_EQ(t2, t0);
    tuple_t t3(std::move(t2));
    EXPECT_EQ(t3, t0);
}

// copy assignment
TEST(pod_tuple, cp_assgn)
{
    tuple_t t1{4, 7l, 3.0f};
    tuple_t t2;
    tuple_t t3;

    t2 = t1;
    t3 = t1;
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

// move assignment
TEST(pod_tuple, mv_assgn)
{
    tuple_t t0{4, 7l, 3.0f};
    tuple_t t1{4, 7l, 3.0f};
    tuple_t t2;
    tuple_t t3;
    t2 = std::move(t1);
    EXPECT_EQ(t2, t0);
    t3 = std::move(t2);
    EXPECT_EQ(t3, t0);
}

// swap
TEST(pod_tuple, swap)
{
    tuple_t t0{4, 7l, 3.0f};
    tuple_t t1{4, 7l, 3.0f};
    tuple_t t2{};
    tuple_t t3{};

    std::swap(t1, t2);
    EXPECT_EQ(t2, t0);
    EXPECT_EQ(t1, t3);
}

// get<1>
TEST(pod_tuple, get_i)
{
    tuple_t t0{4, 7l, 3.0f};

    static_assert(std::is_same_v<decltype(seqan3::get<0>(t0)), int &>);
    static_assert(std::is_same_v<decltype(seqan3::get<1>(t0)), long &>);
    static_assert(std::is_same_v<decltype(seqan3::get<2>(t0)), float &>);
    EXPECT_EQ(seqan3::get<0>(t0), 4);
    EXPECT_EQ(seqan3::get<1>(t0), 7l);
    EXPECT_EQ(seqan3::get<2>(t0), 3.0f);
}

// std::get<1>
TEST(pod_tuple, stdget_i)
{
    tuple_t t0{4, 7l, 3.0f};

    static_assert(std::is_same_v<decltype(std::get<0>(t0)), int &>);
    static_assert(std::is_same_v<decltype(std::get<1>(t0)), long &>);
    static_assert(std::is_same_v<decltype(std::get<2>(t0)), float &>);
    EXPECT_EQ(std::get<0>(t0), 4);
    EXPECT_EQ(std::get<1>(t0), 7l);
    EXPECT_EQ(std::get<2>(t0), 3.0f);
}

// structured bindings
TEST(pod_tuple, struct_binding)
{
    tuple_t t0{4, 7l, 3.0f};
    auto [i, l, f] = t0;

    EXPECT_EQ(i, 4);
    EXPECT_EQ(l, 7l);
    EXPECT_EQ(f, 3.0f);
}

// get<type>
TEST(pod_tuple, get_type)
{
    using pt = tuple_t;
    using ptc = pt const;
    pt t0{4, 7l, 3.0f};
    ptc t1{4, 7l, 3.0f};

    static_assert(std::is_same_v<decltype(seqan3::get<int>(t0)), int &>);
    static_assert(std::is_same_v<decltype(seqan3::get<long>(t0)), long &>);
    static_assert(std::is_same_v<decltype(seqan3::get<float>(t0)), float &>);

    static_assert(std::is_same_v<decltype(seqan3::get<int>(t1)), int const &>);
    static_assert(std::is_same_v<decltype(seqan3::get<long>(t1)), long const &>);
    static_assert(std::is_same_v<decltype(seqan3::get<float>(t1)), float const &>);

    static_assert(std::is_same_v<decltype(seqan3::get<int>(pt{4, 7l, 3.0f})), int &&>);
    static_assert(std::is_same_v<decltype(seqan3::get<long>(pt{4, 7l, 3.0f})), long &&>);
    static_assert(std::is_same_v<decltype(seqan3::get<float>(pt{4, 7l, 3.0f})), float &&>);

    static_assert(std::is_same_v<decltype(seqan3::get<int>(ptc{4, 7l, 3.0f})), int const &&>);
    static_assert(std::is_same_v<decltype(seqan3::get<long>(ptc{4, 7l, 3.0f})), long const &&>);
    static_assert(std::is_same_v<decltype(seqan3::get<float>(ptc{4, 7l, 3.0f})), float const &&>);

    EXPECT_EQ(seqan3::get<int>(t0), 4);
    EXPECT_EQ(seqan3::get<long>(t0), 7l);
    EXPECT_EQ(seqan3::get<float>(t0), 3.0f);

    EXPECT_EQ(seqan3::get<int>(t1), 4);
    EXPECT_EQ(seqan3::get<long>(t1), 7l);
    EXPECT_EQ(seqan3::get<float>(t1), 3.0f);

    EXPECT_EQ(seqan3::get<int>(pt{4, 7l, 3.0f}), 4);
    EXPECT_EQ(seqan3::get<long>(pt{4, 7l, 3.0f}), 7l);
    EXPECT_EQ(seqan3::get<float>(pt{4, 7l, 3.0f}), 3.0f);

    EXPECT_EQ(seqan3::get<int>(ptc{4, 7l, 3.0f}), 4);
    EXPECT_EQ(seqan3::get<long>(ptc{4, 7l, 3.0f}), 7l);
    EXPECT_EQ(seqan3::get<float>(ptc{4, 7l, 3.0f}), 3.0f);
}

// std::get<type>
TEST(pod_tuple, stdget_type)
{
    using pt = tuple_t;
    using ptc = pt const;
    pt t0{4, 7l, 3.0f};
    ptc t1{4, 7l, 3.0f};

    static_assert(std::is_same_v<decltype(std::get<int>(t0)), int &>);
    static_assert(std::is_same_v<decltype(std::get<long>(t0)), long &>);
    static_assert(std::is_same_v<decltype(std::get<float>(t0)), float &>);

    static_assert(std::is_same_v<decltype(std::get<int>(t1)), int const &>);
    static_assert(std::is_same_v<decltype(std::get<long>(t1)), long const &>);
    static_assert(std::is_same_v<decltype(std::get<float>(t1)), float const &>);

    static_assert(std::is_same_v<decltype(std::get<int>(pt{4, 7l, 3.0f})), int &&>);
    static_assert(std::is_same_v<decltype(std::get<long>(pt{4, 7l, 3.0f})), long &&>);
    static_assert(std::is_same_v<decltype(std::get<float>(pt{4, 7l, 3.0f})), float &&>);

    static_assert(std::is_same_v<decltype(std::get<int>(ptc{4, 7l, 3.0f})), int const &&>);
    static_assert(std::is_same_v<decltype(std::get<long>(ptc{4, 7l, 3.0f})), long const &&>);
    static_assert(std::is_same_v<decltype(std::get<float>(ptc{4, 7l, 3.0f})), float const &&>);

    EXPECT_EQ(std::get<int>(t0), 4);
    EXPECT_EQ(std::get<long>(t0), 7l);
    EXPECT_EQ(std::get<float>(t0), 3.0f);

    EXPECT_EQ(std::get<int>(t1), 4);
    EXPECT_EQ(std::get<long>(t1), 7l);
    EXPECT_EQ(std::get<float>(t1), 3.0f);

    EXPECT_EQ(std::get<int>(pt{4, 7l, 3.0f}), 4);
    EXPECT_EQ(std::get<long>(pt{4, 7l, 3.0f}), 7l);
    EXPECT_EQ(std::get<float>(pt{4, 7l, 3.0f}), 3.0f);

    EXPECT_EQ(std::get<int>(ptc{4, 7l, 3.0f}), 4);
    EXPECT_EQ(std::get<long>(ptc{4, 7l, 3.0f}), 7l);
    EXPECT_EQ(std::get<float>(ptc{4, 7l, 3.0f}), 3.0f);
}

// std::tuple_element
TEST(pod_tuple, tuple_element)
{
    using pt = tuple_t;

    static_assert(std::is_same_v<std::tuple_element_t<0, pt>, int>);
    static_assert(std::is_same_v<std::tuple_element_t<1, pt>, long>);
    static_assert(std::is_same_v<std::tuple_element_t<2, pt>, float>);
    static_assert(std::tuple_size_v<pt> == 3);
}

// type deduction
TEST(pod_tuple, type_deduce)
{
    seqan3::pod_tuple t0{4, 7l, 3.0f};
    using pt = decltype(t0);

    static_assert(std::is_same_v<std::tuple_element_t<0, pt>, int>);
    static_assert(std::is_same_v<std::tuple_element_t<1, pt>, long>);
    static_assert(std::is_same_v<std::tuple_element_t<2, pt>, float>);
}

// comparison operators
TEST(pod_tuple, cmp)
{
    tuple_t t0{4, 6l, 4.0f};
    tuple_t t1{4, 7l, 3.0f};
    tuple_t t2{4, 7l, 4.0f};

    EXPECT_LT(t0, t1);
    EXPECT_LE(t0, t1);
    EXPECT_LE(t1, t1);
    EXPECT_EQ(t1, t1);
    EXPECT_GE(t1, t1);
    EXPECT_GE(t2, t1);
    EXPECT_GT(t2, t1);
}
