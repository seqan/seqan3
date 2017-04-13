// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
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
// ============================================================================

#include <type_traits>

#include <gtest/gtest.h>

#include <seqan3/core/pod_tuple.hpp>

using namespace seqan3;

// default/zero construction
TEST(pod_tuple_ctr, ctr)
{
    pod_tuple<int, long, float> t1;
}

// aggregate initialization
TEST(pod_tuple_aggr, aggr)
{
    pod_tuple<int, long, float> t1;
    pod_tuple<int, long, float> t2{4, 7l, 3.0f};
    EXPECT_NE(t1, t2);
}

// zero initialization
TEST(pod_tuple_zro, zro)
{
    pod_tuple<int, long, float> t1{0, 0, 0};
    pod_tuple<int, long, float> t2{};

    EXPECT_EQ(t1, t2);
}

// copy construction
TEST(pod_tuple_cp_ctr, cp_ctr)
{
    pod_tuple<int, long, float> t1{4, 7l, 3.0f};
    pod_tuple<int, long, float> t2{t1};
    pod_tuple<int, long, float> t3(t1);
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

// move construction
TEST(pod_tuple_mv_ctr, mv_ctr)
{
    pod_tuple<int, long, float> t0{4, 7l, 3.0f};
    pod_tuple<int, long, float> t1{4, 7l, 3.0f};
    pod_tuple<int, long, float> t2{std::move(t1)};
    EXPECT_EQ(t2, t0);
    pod_tuple<int, long, float> t3(std::move(t2));
    EXPECT_EQ(t3, t0);
}

// copy assignment
TEST(pod_tuple_cp_assgn, cp_assgn)
{
    pod_tuple<int, long, float> t1{4, 7l, 3.0f};
    pod_tuple<int, long, float> t2;
    pod_tuple<int, long, float> t3;

    t2 = t1;
    t3 = t1;
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

// move assignment
TEST(pod_tuple_mv_assgn, mv_assgn)
{
    pod_tuple<int, long, float> t0{4, 7l, 3.0f};
    pod_tuple<int, long, float> t1{4, 7l, 3.0f};
    pod_tuple<int, long, float> t2;
    pod_tuple<int, long, float> t3;
    t2 = std::move(t1);
    EXPECT_EQ(t2, t0);
    t3 = std::move(t2);
    EXPECT_EQ(t3, t0);
}

// swap
TEST(pod_tuple_swap, swap)
{
    pod_tuple<int, long, float> t0{4, 7l, 3.0f};
    pod_tuple<int, long, float> t1{4, 7l, 3.0f};
    pod_tuple<int, long, float> t2{};
    pod_tuple<int, long, float> t3{};

    std::swap(t1, t2);
    EXPECT_EQ(t2, t0);
    EXPECT_EQ(t1, t3);
}

// get<1>
TEST(pod_tuple_get_i, get_i)
{
    pod_tuple<int, long, float> t0{4, 7l, 3.0f};

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
    pod_tuple<int, long, float> t0{4, 7l, 3.0f};

    static_assert(std::is_same_v<decltype(std::get<0>(t0)), int &>);
    static_assert(std::is_same_v<decltype(std::get<1>(t0)), long &>);
    static_assert(std::is_same_v<decltype(std::get<2>(t0)), float &>);
    EXPECT_EQ(std::get<0>(t0), 4);
    EXPECT_EQ(std::get<1>(t0), 7l);
    EXPECT_EQ(std::get<2>(t0), 3.0f);
}

// structured bindings
TEST(pod_tuple_struct_binding, struct_binding)
{
    pod_tuple<int, long, float> t0{4, 7l, 3.0f};
    auto [ i, l, f ] = t0;

    EXPECT_EQ(i, 4);
    EXPECT_EQ(l, 7l);
    EXPECT_EQ(f, 3.0f);
}

// get<type>
TEST(pod_tuple_get_type, get_type)
{
    using pt = pod_tuple<int, long, float>;
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
TEST(pod_tuple_get_type, stdget_type)
{
    using pt = pod_tuple<int, long, float>;
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
TEST(pod_tuple_tuple_element, tuple_element)
{
    using pt = pod_tuple<int, long, float>;

    static_assert(std::is_same_v<std::tuple_element_t<0, pt>, int>);
    static_assert(std::is_same_v<std::tuple_element_t<1, pt>, long>);
    static_assert(std::is_same_v<std::tuple_element_t<2, pt>, float>);
    static_assert(std::tuple_size_v<pt> == 3);
}

// type deduction
TEST(pod_tuple_type_deduce, type_deduce)
{
    pod_tuple t0{4, 7l, 3.0f};
    using pt = decltype(t0);

    static_assert(std::is_same_v<std::tuple_element_t<0, pt>, int>);
    static_assert(std::is_same_v<std::tuple_element_t<1, pt>, long>);
    static_assert(std::is_same_v<std::tuple_element_t<2, pt>, float>);
}

// comparison operators
TEST(pod_tuple_cmp, cmp)
{
    pod_tuple<int, long, float> t0{4, 6l, 4.0f};
    pod_tuple<int, long, float> t1{4, 7l, 3.0f};
    pod_tuple<int, long, float> t2{4, 7l, 4.0f};

    EXPECT_LT(t0, t1);
    EXPECT_LE(t0, t1);
    EXPECT_LE(t1, t1);
    EXPECT_EQ(t1, t1);
    EXPECT_GE(t1, t1);
    EXPECT_GE(t2, t1);
    EXPECT_GT(t2, t1);
}
