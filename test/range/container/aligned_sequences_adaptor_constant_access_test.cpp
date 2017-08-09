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

#include <sstream>
#include <vector>

#include <gtest/gtest.h>


#include <sdsl/bit_vectors.hpp>
#include "sdsl/rank_support_v5.hpp"
#include "sdsl/select_support_mcl.hpp"
#include "sdsl/util.hpp"

#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/container/aligned_sequence_adaptor_constant_access.hpp>
#include <seqan3/range/detail/random_access_iterator.hpp>

using namespace seqan3;
using sequence_type = aligned_sequence_adaptor_constant_access<std::vector<gapped<dna4>>>;

TEST(aligned_sequences_test, bit_vector_mv)
{
    // bit vector test
    sdsl::bit_vector b1(2, 0), b2(1, 1);
    b1 = std::move(b2); // mv assign
    EXPECT_EQ(1u, b1[0]);
    EXPECT_EQ(1u, b1.size());

    sdsl::bit_vector b3(std::move(b1)); // mv construct
    EXPECT_EQ(1u, b3[0]);
    EXPECT_EQ(1u, b3.size());

    // rank support test
    b1 = sdsl::bit_vector(2, 0);
    b2 = sdsl::bit_vector(1, 1);
    sdsl::rank_support_v5<> rs1(&b1);
    sdsl::rank_support_v5<> rs2(&b2);
    // mv assign
    sdsl::rank_support_v5<> rs3 = std::move(rs1);
    EXPECT_EQ(2u, rs3.size());
    EXPECT_EQ(0u, rs3.rank(0));
    EXPECT_EQ(0u, rs3.rank(1));
    // mv construct
    sdsl::rank_support_v5<> rs4(std::move(rs2));
    EXPECT_EQ(1u, rs4.size());
    EXPECT_EQ(1u, rs4.rank(1));
}

// default (cp, mv) (de)constructors
TEST(aligned_sequences_test, constructors)
{
    sequence_type o, p, q;
    sequence_type r; // but not r()
    sequence_type s(r);
    sequence_type t = q;
    sequence_type u(std::move(p));
    sequence_type v = std::move(o);
    sequence_type * seq_ptr = &p;
    seq_ptr->sequence_type::~sequence_type();
}

// constructors and assignment
TEST(aligned_sequences_test, constructor_by_symbol)
{
    aligned_sequence_adaptor_constant_access<std::vector<gapped<dna4>>> as(1024, gap::GAP);
    auto it = as.begin();
    EXPECT_EQ(gapped<dna4>(gap::GAP), it[0]);
    EXPECT_EQ(gapped<dna4>(gap::GAP), it[1023]);

    aligned_sequence_adaptor_constant_access<std::vector<gapped<dna4>>> as2(1024, dna4::T);
    it = as2.begin();
    EXPECT_EQ(gapped<dna4>(dna4::T), it[0]);
    EXPECT_EQ(gapped<dna4>(dna4::T), it[1023]);
}

// constructor by sequence and random access
TEST(aligned_sequences_test, constructor_by_sequence)
{
    std::vector<gapped<dna4>> sequence = {dna4::A, dna4::C, gap::GAP, dna4::G, dna4::T};
    aligned_sequence_adaptor_constant_access<std::vector<gapped<dna4>>> as(sequence);
    aligned_sequence_adaptor_constant_access<std::vector<gapped<dna4>>> as2(as.begin(), as.end());
    // test as
    detail::random_access_iterator<sequence_type> it = as.begin();
    EXPECT_EQ(gapped<dna4>{dna4::A}, *it);
    it += 2;
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *it);
    // test as2
    it = as2.begin();
    EXPECT_EQ(gapped<dna4>{dna4::A}, *it);
    it += 2;
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *it);
}

// constructor by initializer list
TEST(aligned_sequences_test, constructor_by_init_List)
{
    aligned_sequence_adaptor_constant_access<std::vector<gapped<dna4>>> as({dna4::A, gap::GAP});
    detail::random_access_iterator<sequence_type> it = as.begin();
    EXPECT_EQ(gapped<dna4>{dna4::A}, *it++);
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *it);
}

// move constructors
TEST(aligned_sequences_test, constructor_move)
{
    // move construction with empty sequence
    aligned_sequence_adaptor_constant_access<std::vector<gapped<dna4>>> as_base;
    aligned_sequence_adaptor_constant_access<std::vector<gapped<dna4>>> as_derived(std::move(as_base));
    // move construction with non-empty sequence
    std::vector<gapped<dna4>> seq = {dna4::A, dna4::C, dna4::G, gap::GAP, dna4::T};
    aligned_sequence_adaptor_constant_access<std::vector<gapped<dna4>>> as_base2(seq);
    aligned_sequence_adaptor_constant_access<std::vector<gapped<dna4>>> as_derived2(std::move(as_base2));
    auto it2 = as_derived2.begin();
    EXPECT_EQ(gapped<dna4>{dna4::A}, *it2++);
    EXPECT_EQ(gapped<dna4>{dna4::C}, *it2++);
    EXPECT_EQ(gapped<dna4>{dna4::G}, *it2++);
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *it2++);
    EXPECT_EQ(gapped<dna4>{dna4::T}, *it2);

    // move construction via assignment
    std::vector<gapped<dna4>> seq2 = {dna4::T, dna4::A, gap::GAP};
    aligned_sequence_adaptor_constant_access<std::vector<gapped<dna4>>> as_base3(seq2), as_derived3(seq);
    as_derived3 = std::move(as_base3);
    auto it3 = as_derived3.begin();
    auto it3_end = as_derived3.end();
    EXPECT_EQ(gapped<dna4>{dna4::T}, *it3++);
    EXPECT_EQ(gapped<dna4>{dna4::A}, *it3++);
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *it3++);
    EXPECT_EQ(it3, it3_end);
}

// Destructor
TEST(aligned_sequences_test, destructor)
{
    std::vector<gapped<dna4>> seq = {dna4::T, dna4::A};
    using aligned_sequence_t = aligned_sequence_adaptor_constant_access<std::vector<gapped<dna4>>>;
    aligned_sequence_t s(seq);
    aligned_sequence_t * s_ptr = new aligned_sequence_t(seq);//&s;
    // aligned_sequences_adaptor_constant_access_test.cpp:172:12: error: 'using aligned_sequence_t = struct seqan3::aligned_sequence_adaptor_constant_access<std::vector<seqan3::gapped<seqan3::dna4> > >' is private within this context
    // s_ptr->aligned_sequence_t::~aligned_sequence_t();
    s_ptr->aligned_sequence_adaptor_constant_access<std::vector<gapped<dna4>>>::~aligned_sequence_adaptor_constant_access<std::vector<gapped<dna4>>>();
}

/*
// Container concept functions
TEST(aligned_sequences_test, container_concepts_iterators)
{
    std::vector<gapped<dna4>> seq = {gap::GAP, dna4::T, dna4::A}, seq2 = {dna4::C};
    aligned_sequence_adaptor_constant_access<std::vector<gapped<dna4>>> s(seq);

    // begin and end iterators
    auto it = s.begin();
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *it);
    it = s.end();
    EXPECT_EQ(gapped<dna4>{dna4::A}, *(--it));

    // const begin and const end iterators
    aligned_sequence_adaptor_constant_access<const std::vector<gapped<dna4>>> s_const(seq);
    auto it_const = s_const.cbegin();
    EXPECT_EQ(gapped<dna4>{dna4::T}, *it_const);
    it_const = s_const.cend();
    EXPECT_EQ(gapped<dna4>{dna4::A}, *(--it_const));
}
*/

/*
TEST(aligned_sequences_test, container_concepts_boolean)
{
    std::vector<gapped<dna4>> seq = {gap::GAP, dna4::T, dna4::A}, seq2 = {dna4::C};
    aligned_sequence_adaptor_constant_access<std::vector<dna4>> s(seq), t(seq);
    EXPECT_TRUE(s == t);
    aligned_sequence_adaptor_constant_access<std::vector<dna4>> u(seq2);
    EXPECT_TRUE(t != u);
}

TEST(aligned_sequences_test, container_concepts_swap)
{
    std::vector<gapped<dna4>> seq = {gap::GAP, dna4::T, dna4::A}, seq2 = {dna4::C};
    aligned_sequence_adaptor_constant_access<std::vector<dna4>> t(seq), u(seq2);
    // swap, size, max_size, empty
    t.swap(u);
    auto it = t.begin();
    EXPECT_EQ(gapped<dna4>>{dna4::C}, *it);
    it = u.begin();
    EXPECT_EQ(gapped<dna4>>{gap::GAP}, *it++);
    EXPECT_EQ(gapped<dna4>>{dna4::T}, *it++);
    EXPECT_EQ(gapped<dna4>>{dna4::A}, *it);
    swap(t, u); // friend
    it = t.begin();
    EXPECT_EQ(gapped<dna4>>{dna4::T}, *iterator++);
    EXPECT_EQ(gapped<dna4>>{dna4::A}, *iterator);
    iterator = u.begin();
    EXPECT_EQ(gapped<dna4>>{dna4::C}, *iterator);

    aligned_sequence_adaptor_constant_access<std::vector<dna4>>::size_type max_size = s.max_size();
    EXPECT_FALSE(s.empty());
    aligned_sequence_adaptor_constant_access<std::vector<gapped<dna4>>>> s_empty;
    EXPECT_TRUE(s_empty.empty());
}
*/

/*
// Sequence concept functions
TEST(aligned_sequences_test, sequence_concepts)
{
    using aligned_sequence = aligned_sequence_adaptor_constant_access<std::vector<dna4>>;
    using size_type = aligned_sequence::size_type;
    // construction by repeated value
    size_type size = 1024;
    aligned_sequence(size, dna4::A);
    // construction by sequence given by two iterators
    std::vector<dna4> source = {dna4::A, dna4::C, dna4::G, dna4::T, dna4::G};
    using iterator = detail::random_access_iterator<std::vector<dna4>>;
    iterator begin(source, 0), end(source, source.size());
    aligned_sequence(begin, end);
    // construction by initializer_list
    aligned_sequence s{dna4::A, dna4::A, dna4::C, dna4::T};
    // TODO: query inner sequence
    // assignment via iterators
    s.assign(iterator(source, 0), iterator(source, source.size()));
    // assignment with initializer_list
    s.assign({dna4::T, dna4::C, dna4::G, dna4::A});
    // assignment with value inserted x times
    s.assign(128, dna4::G);
    // insert by value
    aligned_sequence t{dna4::A};
    t.insert(dna4::T);
}
*/
