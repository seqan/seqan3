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

#include <gtest/gtest.h>

#include <seqan3/alphabet/alphabet.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/illumina18.hpp>
#include <seqan3/alphabet/quality/composition.hpp>

using namespace seqan3;

// default/zero construction
TEST(quality_composition_ctr, ctr)
{
    quality_composition<dna4, illumina18> t1;
}

// aggregate initialization
TEST(quality_composition_aggr, aggr)
{
    quality_composition<dna4, illumina18> t1;
    quality_composition<dna4, illumina18> t2{dna4::C, 7};
    EXPECT_NE(t1, t2);
}

// zero initialization
TEST(quality_composition_zro, zro)
{
    quality_composition<dna4, illumina18> t1{dna4::A, 0};
    quality_composition<dna4, illumina18> t2{};

    EXPECT_EQ(t1, t2);
}

// copy construction
TEST(quality_composition_cp_ctr, cp_ctr)
{
    quality_composition<dna4, illumina18> t1{dna4::C, 7};
    quality_composition<dna4, illumina18> t2{t1};
    quality_composition<dna4, illumina18> t3(t1);
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

// move construction
TEST(quality_composition_mv_ctr, mv_ctr)
{
    quality_composition<dna4, illumina18> t0{dna4::C, 7};
    quality_composition<dna4, illumina18> t1{dna4::C, 7};
    quality_composition<dna4, illumina18> t2{std::move(t1)};
    EXPECT_EQ(t2, t0);
    quality_composition<dna4, illumina18> t3(std::move(t2));
    EXPECT_EQ(t3, t0);
}

// copy assignment
TEST(quality_composition_cp_assgn, cp_assgn)
{
    quality_composition<dna4, illumina18> t1{dna4::C, 7};
    quality_composition<dna4, illumina18> t2;
    quality_composition<dna4, illumina18> t3;

    t2 = t1;
    t3 = t1;
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

// move assignment
TEST(quality_composition_mv_assgn, mv_assgn)
{
    quality_composition<dna4, illumina18> t0{dna4::C, 7};
    quality_composition<dna4, illumina18> t1{dna4::C, 7};
    quality_composition<dna4, illumina18> t2;
    quality_composition<dna4, illumina18> t3;
    t2 = std::move(t1);
    EXPECT_EQ(t2, t0);
    t3 = std::move(t2);
    EXPECT_EQ(t3, t0);
}

// swap
TEST(quality_composition_swap, swap)
{
    quality_composition<dna4, illumina18> t0{dna4::C, 7};
    quality_composition<dna4, illumina18> t1{dna4::C, 7};
    quality_composition<dna4, illumina18> t2{};
    quality_composition<dna4, illumina18> t3{};

    std::swap(t1, t2);
    EXPECT_EQ(t2, t0);
    EXPECT_EQ(t1, t3);
}

// std::get<1>
TEST(quality_composition_get_i, get_i)
{
    quality_composition<dna4, illumina18> t0{dna4::C, 7};

    static_assert(std::is_same_v<decltype(std::get<0>(t0)), dna4 &>);
    static_assert(std::is_same_v<decltype(std::get<1>(t0)), illumina18 &>);

    EXPECT_EQ(std::get<0>(t0), dna4::C);
    EXPECT_EQ(std::get<1>(t0), illumina18{7});
}

// structured bindings
TEST(quality_composition_struct_binding, struct_binding)
{
    quality_composition<dna4, illumina18> t0{dna4::C, 7};
    auto [ i, l ] = t0;

    static_assert(std::is_same_v<decltype(i), dna4>);
    //WARNING HERE BE DRAGONS:
//     static_assert(std::is_same_v<decltype(l), illumina18>);

    EXPECT_EQ(i, dna4::C);
//     EXPECT_EQ(l, illumina18{7});
}

// std::get<type> NOT IMPLEMENTED YET
// TEST(quality_composition_get_type, get_type)
// {
//     quality_composition<dna4, illumina18> t0{dna4::C, 7};
//
//     EXPECT_EQ(std::get<int>(t0), 4);
//     EXPECT_EQ(std::get<long>(t0), 7l);
//     EXPECT_EQ(std::get<float>(t0), 3.0f);
// }

// std::tuple_element
TEST(quality_composition_tuple_element, tuple_element)
{
    using pt = quality_composition<dna4, illumina18>;

    static_assert(std::is_same_v<std::tuple_element_t<0, pt>, dna4>);
    static_assert(std::is_same_v<std::tuple_element_t<1, pt>, illumina18>);
    //TODO generalize this, too:
//     static_assert(std::tuple_size_v<pt> == 2);
}

// type deduction
// TEST(quality_composition_type_deduce, type_deduce)
// {
//     quality_composition t0{dna4::C, illumina18{7}};
//     using pt = decltype(t0);
//
//     static_assert(std::is_same_v<std::tuple_element_t<0, pt>, dna4>);
//     static_assert(std::is_same_v<std::tuple_element_t<1, pt>, illumina18>);
//     static_assert(std::tuple_size_v<pt> == 2);
// }

// comparison operators illumina18 comparison broken?
TEST(quality_composition_cmp, cmp)
{
    quality_composition<dna4, illumina18> t0{dna4::C, 6};
    quality_composition<dna4, illumina18> t1{dna4::C, 7};
    quality_composition<dna4, illumina18> t2{dna4::G, 7};

    EXPECT_LT(t0, t1);
    EXPECT_LE(t0, t1);
    EXPECT_LE(t1, t1);
    EXPECT_EQ(t1, t1);
    EXPECT_GE(t1, t1);
    EXPECT_GE(t2, t1);
    EXPECT_GT(t2, t1);
}
