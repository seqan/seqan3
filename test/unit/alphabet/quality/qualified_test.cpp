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
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>

using namespace seqan3;

/************** TUPLE INHERITANCE **********************/

// default/zero construction
TEST(qualified, ctr)
{
    [[maybe_unused]] qualified<dna4, phred42> t1{};
}

// aggregate initialization
TEST(qualified, aggr)
{
    qualified<dna4, phred42> t1{};
    qualified<dna4, phred42> t2{dna4::C, 7};
    EXPECT_NE(t1, t2);
}

// zero initialization
TEST(qualified, zro)
{
    qualified<dna4, phred42> t1{dna4::A, 0};
    qualified<dna4, phred42> t2{};

    EXPECT_EQ(t1, t2);
}

// copy construction
TEST(qualified, cp_ctr)
{
    qualified<dna4, phred42> t1{dna4::C, 7};
    qualified<dna4, phred42> t2{t1};
    qualified<dna4, phred42> t3(t1);
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

// move construction
TEST(qualified, mv_ctr)
{
    qualified<dna4, phred42> t0{dna4::C, 7};
    qualified<dna4, phred42> t1{dna4::C, 7};
    qualified<dna4, phred42> t2{std::move(t1)};
    EXPECT_EQ(t2, t0);
    qualified<dna4, phred42> t3(std::move(t2));
    EXPECT_EQ(t3, t0);
}

// copy assignment
TEST(qualified, cp_assgn)
{
    qualified<dna4, phred42> t1{dna4::C, 7};
    qualified<dna4, phred42> t2{};
    qualified<dna4, phred42> t3{};

    t2 = t1;
    t3 = t1;
    EXPECT_EQ(t1, t2);
    EXPECT_EQ(t2, t3);
}

// move assignment
TEST(qualified, mv_assgn)
{
    qualified<dna4, phred42> t0{dna4::C, 7};
    qualified<dna4, phred42> t1{dna4::C, 7};
    qualified<dna4, phred42> t2{};
    qualified<dna4, phred42> t3{};
    t2 = std::move(t1);
    EXPECT_EQ(t2, t0);
    t3 = std::move(t2);
    EXPECT_EQ(t3, t0);
}

// swap
TEST(qualified, swap)
{
    qualified<dna4, phred42> t0{dna4::C, 7};
    qualified<dna4, phred42> t1{dna4::C, 7};
    qualified<dna4, phred42> t2{};
    qualified<dna4, phred42> t3{};

    std::swap(t1, t2);
    EXPECT_EQ(t2, t0);
    EXPECT_EQ(t1, t3);
}

// get<1>
TEST(qualified, get_i)
{
    qualified<dna4, phred42> t0{dna4::C, 7};

    static_assert(std::is_same_v<decltype(seqan3::get<0>(t0)), dna4 &>);
    static_assert(std::is_same_v<decltype(seqan3::get<1>(t0)), phred42 &>);

    EXPECT_EQ(seqan3::get<0>(t0), dna4::C);
    EXPECT_EQ(seqan3::get<1>(t0), phred42{7});
}

// std::get<1>
TEST(qualified, stdget_i)
{
    qualified<dna4, phred42> t0{dna4::C, 7};

    static_assert(std::is_same_v<decltype(std::get<0>(t0)), dna4 &>);
    static_assert(std::is_same_v<decltype(std::get<1>(t0)), phred42 &>);

    EXPECT_EQ(std::get<0>(t0), dna4::C);
    EXPECT_EQ(std::get<1>(t0), phred42{7});
}

// structured bindings
TEST(qualified, struct_binding)
{
    qualified<dna4, phred42> t0{dna4::C, 7};
    auto [ i, l ] = t0;

    static_assert(std::is_same_v<decltype(i), dna4>);
    static_assert(std::is_same_v<decltype(l), phred42>);

    EXPECT_EQ(i, dna4::C);
    EXPECT_EQ(l, phred42{7});
}

// get<type>
TEST(qualified, get_type)
{
    qualified<dna4, phred42> t0{dna4::C, 7};

    EXPECT_EQ(seqan3::get<dna4>(t0), dna4::C);
    EXPECT_EQ(seqan3::get<phred42>(t0), phred42{7});
}

// std::get<type>
TEST(qualified, stdget_type)
{
    qualified<dna4, phred42> t0{dna4::C, 7};

    EXPECT_EQ(std::get<dna4>(t0), dna4::C);
    EXPECT_EQ(std::get<phred42>(t0), phred42{7});;
}

// std::tuple_element
TEST(qualified, tuple_element)
{
    using pt = qualified<dna4, phred42>;

    static_assert(std::is_same_v<std::tuple_element_t<0, pt>, dna4>);
    static_assert(std::is_same_v<std::tuple_element_t<1, pt>, phred42>);
    static_assert(std::tuple_size_v<pt> == 2);
}

// type deduction
TEST(qualified, type_deduce)
{
    qualified t0{dna4::C, phred42{7}};
    using pt = decltype(t0);

    static_assert(std::is_same_v<std::tuple_element_t<0, pt>, dna4>);
    static_assert(std::is_same_v<std::tuple_element_t<1, pt>, phred42>);
    static_assert(std::tuple_size_v<pt> == 2);
}

// explicit cast to element
TEST(qualified, cast_to_element)
{
    qualified<dna4, phred42> t0{dna4::C, 7};

    auto d = static_cast<dna4>(t0);
    auto q = static_cast<phred42>(t0);
    static_assert(std::is_same_v<decltype(d), dna4>);
    static_assert(std::is_same_v<decltype(q), phred42>);

    EXPECT_EQ(d, dna4::C);
    EXPECT_EQ(q, phred42{7});
}

// comparison operators phred42 comparison broken?
TEST(qualified, cmp)
{
    qualified<dna4, phred42> t0{dna4::C, 6};
    qualified<dna4, phred42> t1{dna4::C, 7};
    qualified<dna4, phred42> t2{dna4::G, 7};

    EXPECT_LT(t0, t1);
    EXPECT_LE(t0, t1);
    EXPECT_LE(t1, t1);
    EXPECT_EQ(t1, t1);
    EXPECT_GE(t1, t1);
    EXPECT_GE(t2, t1);
    EXPECT_GT(t2, t1);
}

/************** ALPHABET and QUALITY concept **********************/

TEST(qualified, rank_type)
{
    EXPECT_TRUE((std::is_same_v<underlying_rank_t<qualified<dna4, phred42>>,
                               uint8_t>));
}

TEST(qualified, char_type)
{
    EXPECT_TRUE((std::is_same_v<underlying_char_t<qualified<dna4, phred42>>,
                               underlying_char_t<dna4>>));
}

TEST(qualified, phred_type)
{
    EXPECT_TRUE((std::is_same_v<underlying_phred_t<qualified<dna4, phred42>>,
                               underlying_phred_t<phred42>>));
}

TEST(qualified, alphabet_size_v)
{
    EXPECT_EQ((alphabet_size_v<qualified<dna4, phred42>>),
              (alphabet_size_v<dna4> * alphabet_size_v<phred42>));
}

TEST(qualified, to_rank)
{
    qualified<dna4, phred42> t0{dna4::C, 6};
    EXPECT_EQ(to_rank(std::get<0>(t0)), 1);
    EXPECT_EQ(to_rank(std::get<1>(t0)), 6);
    EXPECT_EQ(to_rank(t0),
              to_rank(std::get<0>(t0)) +
              alphabet_size_v<dna4> * to_rank(std::get<1>(t0)));
}

TEST(qualified, assign_rank)
{
    using type = qualified<dna4, phred42>;

    type t0{};

    for (underlying_rank_t<type> i = 0; i < alphabet_size_v<type>; ++i)
    {
        assign_rank(t0, i);
        EXPECT_EQ(to_rank(t0), i);
    }
}

TEST(qualified, to_char)
{
    qualified<dna4, phred42> t0{dna4::C, 6};
    EXPECT_EQ(to_char(std::get<0>(t0)), 'C');
    EXPECT_EQ(to_char(std::get<1>(t0)), '!' + 6);
    EXPECT_EQ(to_char(t0), 'C');
}

TEST(qualified, assign_char)
{
    using type = qualified<dna4, phred42>;

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

TEST(qualified, to_phred)
{
    qualified<dna4, phred42> t0{dna4::C, 6};
    EXPECT_EQ(to_phred(std::get<1>(t0)), 6);
    EXPECT_EQ(to_phred(t0), 6);
}

TEST(qualified, assign_phred)
{
    using type = qualified<dna4, phred42>;

    type t0{dna4::C, 17};
    char schar = to_char(t0);

    assign_phred(t0, 12);
    EXPECT_EQ(to_phred(t0), 12);
    EXPECT_EQ(to_char(t0), schar);
    assign_phred(t0, 37);
    EXPECT_EQ(to_phred(t0), 37);
    EXPECT_EQ(to_char(t0), schar);
}

TEST(qualified, outstream)
{
    qualified<dna4, phred42> t0{dna4::C, 6};
    std::stringstream s;
    s << t0;
    t0 = dna4::A;
    s << t0;

    EXPECT_EQ(s.str(), "CA");
}
