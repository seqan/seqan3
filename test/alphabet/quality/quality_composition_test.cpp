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

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/illumina18.hpp>
#include <seqan3/alphabet/quality/composition.hpp>

using namespace seqan3;

/************** TUPLE INHERITANCE **********************/

// default/zero construction
TEST(quality_composition, ctr)
{
    quality_composition<dna4, illumina18> t1;
}

// aggregate initialization
TEST(quality_composition, aggr)
{
    quality_composition<dna4, illumina18> t1;
    quality_composition<dna4, illumina18> t2{dna4::C, 7};
    EXPECT_NE(t1, t2);
}

// zero initialization
TEST(quality_composition, zro)
{
    quality_composition<dna4, illumina18> t1{dna4::A, 0};
    quality_composition<dna4, illumina18> t2{};

    EXPECT_EQ(t1, t2);
}

// copy construction
TEST(quality_composition, cp_ctr)
{
    quality_composition<dna4, illumina18> t1{dna4::C, 7};
    quality_composition<dna4, illumina18> t2{t1};
    quality_composition<dna4, illumina18> t3(t1);
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

// move construction
TEST(quality_composition, mv_ctr)
{
    quality_composition<dna4, illumina18> t0{dna4::C, 7};
    quality_composition<dna4, illumina18> t1{dna4::C, 7};
    quality_composition<dna4, illumina18> t2{std::move(t1)};
    EXPECT_EQ(t2, t0);
    quality_composition<dna4, illumina18> t3(std::move(t2));
    EXPECT_EQ(t3, t0);
}

// copy assignment
TEST(quality_composition, cp_assgn)
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
TEST(quality_composition, mv_assgn)
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
TEST(quality_composition, swap)
{
    quality_composition<dna4, illumina18> t0{dna4::C, 7};
    quality_composition<dna4, illumina18> t1{dna4::C, 7};
    quality_composition<dna4, illumina18> t2{};
    quality_composition<dna4, illumina18> t3{};

    std::swap(t1, t2);
    EXPECT_EQ(t2, t0);
    EXPECT_EQ(t1, t3);
}

// get<1>
TEST(quality_composition, get_i)
{
    quality_composition<dna4, illumina18> t0{dna4::C, 7};

    static_assert(std::is_same_v<decltype(seqan3::get<0>(t0)), dna4 &>);
    static_assert(std::is_same_v<decltype(seqan3::get<1>(t0)), illumina18 &>);

    EXPECT_EQ(seqan3::get<0>(t0), dna4::C);
    EXPECT_EQ(seqan3::get<1>(t0), illumina18{7});
}

// std::get<1>
TEST(quality_composition, stdget_i)
{
    quality_composition<dna4, illumina18> t0{dna4::C, 7};

    static_assert(std::is_same_v<decltype(std::get<0>(t0)), dna4 &>);
    static_assert(std::is_same_v<decltype(std::get<1>(t0)), illumina18 &>);

    EXPECT_EQ(std::get<0>(t0), dna4::C);
    EXPECT_EQ(std::get<1>(t0), illumina18{7});
}

// structured bindings
TEST(quality_composition, struct_binding)
{
    quality_composition<dna4, illumina18> t0{dna4::C, 7};
    auto [ i, l ] = t0;

    static_assert(std::is_same_v<decltype(i), dna4>);
    static_assert(std::is_same_v<decltype(l), illumina18>);

    EXPECT_EQ(i, dna4::C);
    EXPECT_EQ(l, illumina18{7});
}

// get<type>
TEST(quality_composition, get_type)
{
    quality_composition<dna4, illumina18> t0{dna4::C, 7};

    EXPECT_EQ(seqan3::get<dna4>(t0), dna4::C);
    EXPECT_EQ(seqan3::get<illumina18>(t0), illumina18{7});
}

// std::get<type>
TEST(quality_composition, stdget_type)
{
    quality_composition<dna4, illumina18> t0{dna4::C, 7};

    EXPECT_EQ(std::get<dna4>(t0), dna4::C);
    EXPECT_EQ(std::get<illumina18>(t0), illumina18{7});;
}

// std::tuple_element
TEST(quality_composition, tuple_element)
{
    using pt = quality_composition<dna4, illumina18>;

    static_assert(std::is_same_v<std::tuple_element_t<0, pt>, dna4>);
    static_assert(std::is_same_v<std::tuple_element_t<1, pt>, illumina18>);
    static_assert(std::tuple_size_v<pt> == 2);
}

// type deduction
TEST(quality_composition, type_deduce)
{
    quality_composition t0{dna4::C, illumina18{7}};
    using pt = decltype(t0);

    static_assert(std::is_same_v<std::tuple_element_t<0, pt>, dna4>);
    static_assert(std::is_same_v<std::tuple_element_t<1, pt>, illumina18>);
    static_assert(std::tuple_size_v<pt> == 2);
}

// explicit cast to element
TEST(quality_composition, cast_to_element)
{
    quality_composition<dna4, illumina18> t0{dna4::C, 7};

    auto d = static_cast<dna4>(t0);
    auto q = static_cast<illumina18>(t0);
    static_assert(std::is_same_v<decltype(d), dna4>);
    static_assert(std::is_same_v<decltype(q), illumina18>);

    EXPECT_EQ(d, dna4::C);
    EXPECT_EQ(q, illumina18{7});
}

// comparison operators illumina18 comparison broken?
TEST(quality_composition, cmp)
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

/************** ALPHABET and QUALITY concept **********************/

TEST(quality_composition, rank_type)
{
    EXPECT_TRUE((std::is_same_v<underlying_rank_t<quality_composition<dna4, illumina18>>,
                               uint8_t>));
}

TEST(quality_composition, char_type)
{
    EXPECT_TRUE((std::is_same_v<underlying_char_t<quality_composition<dna4, illumina18>>,
                               underlying_char_t<dna4>>));
}

TEST(quality_composition, phred_type)
{
    EXPECT_TRUE((std::is_same_v<underlying_phred_t<quality_composition<dna4, illumina18>>,
                               underlying_phred_t<illumina18>>));
}

TEST(quality_composition, alphabet_size_v)
{
    EXPECT_EQ((alphabet_size_v<quality_composition<dna4, illumina18>>),
              (alphabet_size_v<dna4> * alphabet_size_v<illumina18>));
}

TEST(quality_composition, to_rank)
{
    quality_composition<dna4, illumina18> t0{dna4::C, 6};
    EXPECT_EQ(to_rank(std::get<0>(t0)), 1);
    EXPECT_EQ(to_rank(std::get<1>(t0)), 6);
    EXPECT_EQ(to_rank(t0),
              to_rank(std::get<0>(t0)) +
              alphabet_size_v<dna4> * to_rank(std::get<1>(t0)));
}

TEST(quality_composition, assign_rank)
{
    using type = quality_composition<dna4, illumina18>;

    type t0{};

    for (underlying_rank_t<type> i = 0; i < alphabet_size_v<type>; ++i)
    {
        assign_rank(t0, i);
        EXPECT_EQ(to_rank(t0), i);
    }
}

TEST(quality_composition, to_char)
{
    quality_composition<dna4, illumina18> t0{dna4::C, 6};
    EXPECT_EQ(to_char(std::get<0>(t0)), 'C');
    EXPECT_EQ(to_char(std::get<1>(t0)), '!' + 6);
    EXPECT_EQ(to_char(t0), 'C');
}

TEST(quality_composition, assign_char)
{
    using type = quality_composition<dna4, illumina18>;

    type t0{dna4::C, 17};
    char qchar = to_char(std::get<1>(t0));

    assign_char(t0, 'A');
    EXPECT_EQ(to_char(t0), 'A');
    EXPECT_EQ(to_char(std::get<1>(t0)), qchar);
    assign_char(t0, 'C');
    EXPECT_EQ(to_char(t0), 'C');
    EXPECT_EQ(to_char(std::get<1>(t0)), qchar);
    assign_char(t0, 'G');
    EXPECT_EQ(to_char(t0), 'G');
    EXPECT_EQ(to_char(std::get<1>(t0)), qchar);
    assign_char(t0, 'T');
    EXPECT_EQ(to_char(t0), 'T');
    EXPECT_EQ(to_char(std::get<1>(t0)), qchar);
    assign_char(t0, 'N');
    EXPECT_EQ(to_char(t0), 'A');
    EXPECT_EQ(to_char(std::get<1>(t0)), qchar);
}

TEST(quality_composition, to_phred)
{
    quality_composition<dna4, illumina18> t0{dna4::C, 6};
    EXPECT_EQ(to_phred(std::get<1>(t0)), 6);
    EXPECT_EQ(to_phred(t0), 6);
}

TEST(quality_composition, assign_phred)
{
    using type = quality_composition<dna4, illumina18>;

    type t0{dna4::C, 17};
    char schar = to_char(t0);

    assign_phred(t0, 12);
    EXPECT_EQ(to_phred(t0), 12);
    EXPECT_EQ(to_char(t0), schar);
    assign_phred(t0, 37);
    EXPECT_EQ(to_phred(t0), 37);
    EXPECT_EQ(to_char(t0), schar);
}

TEST(quality_composition, outstream)
{
    quality_composition<dna4, illumina18> t0{dna4::C, 6};
    std::stringstream s;
    s << t0;
    t0 = dna4::A;
    s << t0;

    EXPECT_EQ(s.str(), "CA");
}
