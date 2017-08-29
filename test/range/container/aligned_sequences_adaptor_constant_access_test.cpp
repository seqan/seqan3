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

#include <gtest/gtest-spi.h>
#include "gtest/gtest-death-test.h"
#include "gtest/gtest.h"
#include "gtest/internal/gtest-filepath.h"


#include <sdsl/bit_vectors.hpp>
#include "sdsl/rank_support_v5.hpp"
#include "sdsl/select_support_mcl.hpp"
#include "sdsl/util.hpp"

#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/container/aligned_sequence_adaptor_constant_access.hpp>
#include <seqan3/range/detail/random_access_iterator.hpp>

using namespace seqan3;
using sequence_t = aligned_sequence_adaptor_constant_access<std::vector<gapped<dna4>>>;

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
    sequence_t o, p, q;
    sequence_t r; // but not r()
    sequence_t s(r);
    sequence_t t = q;
    sequence_t u(std::move(p));
    sequence_t v = std::move(o);
    sequence_t * seq_ptr = &p;
    seq_ptr->sequence_t::~sequence_t();
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

// constructor by sequence, iterators, and initializer_list
TEST(aligned_sequences_test, constructor_by_sequence)
{
    std::vector<gapped<dna4>> sequence = {dna4::A, dna4::C, gap::GAP, dna4::G, dna4::T};
    aligned_sequence_adaptor_constant_access<std::vector<gapped<dna4>>> as(sequence);
    aligned_sequence_adaptor_constant_access<std::vector<gapped<dna4>>> as2(as.begin(), as.end());
    aligned_sequence_adaptor_constant_access<std::vector<gapped<dna4>>> as3{dna4::A, dna4::C, gap::GAP, dna4::G, dna4::T};
    // test as
    detail::random_access_iterator<sequence_t> it = as.begin();
    EXPECT_EQ(gapped<dna4>{dna4::A}, *it);
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *(it+2));
    // test as2
    it = as2.begin();
    EXPECT_EQ(gapped<dna4>{dna4::A}, *it);
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *(it+2));
    // test as3
    it = as3.begin();
    EXPECT_EQ(gapped<dna4>{dna4::A}, *it);
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *(it+2));
}

// constructor by initializer list
TEST(aligned_sequences_test, constructor_by_init_List)
{
    aligned_sequence_adaptor_constant_access<std::vector<gapped<dna4>>> as({dna4::A, gap::GAP});
    detail::random_access_iterator<sequence_t> it = as.begin();
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


// Container concept functions
TEST(aligned_sequences_test, container_concepts_iterators)
{
    std::vector<gapped<dna4>> seq = {gap::GAP, dna4::T, dna4::A};
    aligned_sequence_adaptor_constant_access<std::vector<gapped<dna4>>> s(seq);

    // (const) begin and end iterators of non-const sequence
    auto it = s.begin();
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *it);
    it = s.end();
    EXPECT_EQ(gapped<dna4>{dna4::A}, *(--it));

    // (const) begin and end iterators of non-const sequence
    const aligned_sequence_adaptor_constant_access<std::vector<gapped<dna4>>> s_const(seq);
    auto it_const = s_const.cbegin();
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *it_const);
    it_const = s_const.cend();
    EXPECT_EQ(gapped<dna4>{dna4::A}, *(--it_const));
}

TEST(aligned_sequences_test, container_concepts_boolean)
{
    std::vector<gapped<dna4>> seq = {gap::GAP, dna4::T, dna4::A}, seq2 = {dna4::C};
    aligned_sequence_adaptor_constant_access<std::vector<gapped<dna4>>> s(seq), t(seq);
    EXPECT_TRUE(s == t);
    aligned_sequence_adaptor_constant_access<std::vector<gapped<dna4>>> u(seq2);
    EXPECT_TRUE(t != u);
}

// swap, size, max_size, empty
TEST(aligned_sequences_test, container_concepts_swap)
{
    std::vector<gapped<dna4>> seq = {gap::GAP, dna4::T, dna4::A}, seq2 = {dna4::C};
    aligned_sequence_adaptor_constant_access<std::vector<gapped<dna4>>> t(seq), u(seq2);
    t.swap(u);
    auto it = t.begin();
    EXPECT_EQ(gapped<dna4>{dna4::C}, *it);
    it = u.begin();
    std::cout << it[0] << ", " << it[1] << ", " << it[2] << ", " << std::endl;
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *it++);
    EXPECT_EQ(gapped<dna4>{dna4::T}, *it++);
    EXPECT_EQ(gapped<dna4>{dna4::A}, *it);
    swap(t, u); // friend
    it = t.begin();
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *it++);
    EXPECT_EQ(gapped<dna4>{dna4::T}, *it++);
    EXPECT_EQ(gapped<dna4>{dna4::A}, *it);
    it = u.begin();
    EXPECT_EQ(gapped<dna4>{dna4::C}, *it);

    aligned_sequence_adaptor_constant_access<std::vector<gapped<dna4>>>::size_type max_size = t.max_size();
    EXPECT_TRUE(max_size > 0);
    EXPECT_FALSE(t.empty());
    aligned_sequence_adaptor_constant_access<std::vector<gapped<dna4>>> s_empty;
    EXPECT_TRUE(s_empty.empty());
}

// Sequence concept functions
TEST(aligned_sequences_test, sequence_concepts_assign)
{
    using aligned_sequence = aligned_sequence_adaptor_constant_access<std::vector<gapped<dna4>>>;
    std::vector<gapped<dna4>> seq1 = {dna4::A, gap::GAP, dna4::A, dna4::C, dna4::T};
    std::vector<gapped<dna4>> seq2 = {gap::GAP, dna4::G, gap::GAP};

    aligned_sequence s{seq1}, s2{seq1};
    aligned_sequence t{}, t2{};
    aligned_sequence u{seq2};

    // TODO: query inner sequence
    // assignment via iterators of non-empty sequence to empty sequence
    t.assign(s.begin(), s.end());
    EXPECT_EQ(5u, t.size());
    auto it = t.begin();
    for (unsigned int i = 0; i < 5u; ++i)
        EXPECT_EQ(seq1[i], it[i]);

    // assignment via iterators of non-empty to non-empty sequence
    s.assign(u.begin(), u.begin()+2);
    EXPECT_EQ(2u, s.size());
    it = s.begin();
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *it++);
    EXPECT_EQ(gapped<dna4>{dna4::G}, *it);

    // assignment with initializer_list of non-empty to empty sequence
    std::initializer_list<gapped<dna4>> l{dna4::T, dna4::C, dna4::G, dna4::A, gap::GAP}, l2{gap::GAP, gap::GAP};
    t2.assign(l);
    EXPECT_EQ(5u, t2.size());
    it = t2.begin();
    auto it_l = l.begin();
    for (unsigned int i = 0; i < 5u; ++i)
        EXPECT_EQ(*it_l++, *it++);

    // assignment with initializer_list of non-empty to empty sequence
    t2.assign(l2);
    EXPECT_EQ(2u, t2.size());
    it = t2.begin();
    it_l = l2.begin();
    EXPECT_EQ(*it_l++, *it++);
    EXPECT_EQ(*it_l, *it);

    // assignment with value inserted x times
    s.assign(128, dna4::G);
    EXPECT_EQ(128u, s.size());
    it = s.begin();
    EXPECT_EQ(gapped<dna4>{dna4::G}, *it);
    EXPECT_EQ(gapped<dna4>{dna4::G}, *(it+127));

    s.assign(64, gap::GAP);
    EXPECT_EQ(64u, s.size());
    it = s.begin();
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *it);
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *(it+63));
}

// TODO: test return values of insert ops
TEST(aligned_sequences_test, sequence_concepts_insert)
{
    // insert by position iterator into empty sequence
    sequence_t s{};
    s.insert(s.begin(), dna4::C);
    EXPECT_EQ(1u, s.size());
    EXPECT_EQ(gapped<dna4>{dna4::C}, *s.begin());

    // insert by position iterator into non-empty sequence
    std::initializer_list<gapped<dna4>> l{gap::GAP, gap::GAP, dna4::T, dna4::A};
    std::initializer_list<gapped<dna4>> l_post{dna4::C, gap::GAP, gap::GAP, dna4::T, dna4::A};
    sequence_t t(l);
    t.insert(t.begin(), dna4::C);
    EXPECT_EQ(5u, t.size());
    auto l_it = l_post.begin();
    for (auto it = t.begin(); l_it != l.end(); ++l_it, ++ it)
        EXPECT_EQ(*l_it, *it);

    // insert after last element
    sequence_t w(l);
    w.insert(w.end(), dna4::G);
    EXPECT_EQ(5u, w.size());

    // insert by position iterator and moved element into empty sequence
    sequence_t u{};
    u.insert(u.begin(), std::move(dna4::C));
    EXPECT_EQ(1u, u.size());
    EXPECT_EQ(gapped<dna4>{dna4::C}, *s.begin());

    // insert by position iterator and moved element into non-empty sequence
    std::initializer_list<gapped<dna4>> l_post2{gap::GAP, dna4::C, gap::GAP, dna4::T, dna4::A};
    sequence_t v(l);
    v.insert(v.begin()+1, std::move(dna4::C));
    EXPECT_EQ(5u, v.size());
    auto l_it2 = l_post2.begin();
    for (auto it = v.begin(); l_it2 != l_post2.end(); l_it2++, it++)
        EXPECT_EQ(*l_it2, *it);

    // bulk insert into empty sequence
    sequence_t x{};
    x.insert(x.begin(), 16, gap::GAP);
    EXPECT_EQ(16u, x.size());
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *(x.begin()));
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *(x.end()-1));
    // bulk insert into non-empty sequence
    sequence_t y{l};
    y.insert(y.begin()+2, 16, dna4::T);
    EXPECT_EQ(20u, y.size());
    auto it = y.begin();
    EXPECT_EQ(gapped<dna4>{gap::GAP}, y[0]);
    EXPECT_EQ(gapped<dna4>{gap::GAP}, y[1]);
    EXPECT_EQ(gapped<dna4>{dna4::T}, y[2]);
    EXPECT_EQ(gapped<dna4>{dna4::T}, y[17]);
    EXPECT_EQ(gapped<dna4>{dna4::T}, y[18]);
    EXPECT_EQ(gapped<dna4>{dna4::A}, y[19]);

}

// TODO: test return values of erase ops
// TODO: test size_type erase
TEST(aligned_sequences_test, sequence_concepts_erase_one)
{
    std::initializer_list<gapped<dna4>> l{gap::GAP, gap::GAP, dna4::T, dna4::A};
    sequence_t s{l};
    // erase 1st element
    auto it = s.erase(s.begin());
    EXPECT_EQ(it, s.begin());
    EXPECT_EQ(3u, s.size());
    auto it_l = l.begin()+1;
    for (it = s.begin(); it_l != l.end(); ++it, ++it_l)
        EXPECT_EQ(*it_l, *it);
    // erase last element
    it = s.erase(s.end()-1);
    EXPECT_EQ(2u, s.size());
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *s.begin());
    EXPECT_EQ(gapped<dna4>{dna4::T}, *(s.begin()+1));
    // erase midth element
    std::initializer_list<gapped<dna4>> l2{gap::GAP, dna4::T, dna4::A};
    sequence_t t{l2};
    it = t.erase(t.begin()+1);
    EXPECT_EQ(2u, t.size());
    EXPECT_EQ(it, t.begin()+1);
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *t.begin());
    EXPECT_EQ(gapped<dna4>{dna4::A}, *(t.begin()+1));
}

// erase elements in range between two iterators
TEST(aligned_sequences_test, sequence_concepts_erase_range)
{
    // range length = 1
    std::initializer_list<gapped<dna4>> l{gap::GAP, gap::GAP, dna4::T, dna4::A};
    sequence_t s{l};

    // erase 1st element
    auto it = s.erase(s.begin(), s.begin() + 1);
    EXPECT_EQ(it, s.begin());
    EXPECT_EQ(3u, s.size());
    auto it_l = l.begin()+1;
    for (it = s.begin(); it_l != l.end(); ++it, ++it_l){
        std::cout << (*it) << std::endl;
        EXPECT_EQ(*it_l, *it);
    }

    // erase last element
    it = s.erase(s.end()-1, s.end());
    EXPECT_EQ(2u, s.size());
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *s.begin());
    EXPECT_EQ(gapped<dna4>{dna4::T}, *(s.begin()+1));

    // erase midth element
    std::initializer_list<gapped<dna4>> l2{gap::GAP, dna4::T, dna4::A};
    sequence_t t{l2};
    it = t.erase(t.begin()+1, t.begin()+2);
    EXPECT_EQ(2u, t.size());
    EXPECT_EQ(it, t.begin()+1);
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *t.begin());
    EXPECT_EQ(gapped<dna4>{dna4::A}, *(t.begin()+1));

    // range length > 1 in the middle
    sequence_t u{128, dna4::T};
    u.erase(u.begin()+16, u.begin()+32);
    EXPECT_EQ(112u, u.size());

    // erase range mid to end
    sequence_t v{dna4::A, dna4::C, dna4::G, dna4::T, dna4::A, gap::GAP, dna4::G};
    it = v.erase(v.begin()+2, v.end());
    EXPECT_EQ(it, v.end());
    EXPECT_EQ(2u, v.size());
    it = v.begin();
    EXPECT_EQ(gapped<dna4>{dna4::A}, *it++);
    EXPECT_EQ(gapped<dna4>{dna4::C}, *it);

    // erase range begining to mid
    sequence_t w{gap::GAP, dna4::A, dna4::C, gap::GAP, dna4::G, dna4::A};
    it = w.erase(w.begin(), w.begin()+3);
    EXPECT_EQ(it, w.begin());
    EXPECT_EQ(3u, w.size());
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *it++);
    EXPECT_EQ(gapped<dna4>{dna4::G}, *it++);
    EXPECT_EQ(gapped<dna4>{dna4::A}, *it);
}

TEST(aligned_sequences_test, sequence_concepts_push_back)
{
    // case 1: into empty sequence a gap or alphabet symbol
    sequence_t s{};
    s.push_back(gapped<dna4>{gap::GAP});
    EXPECT_EQ(1u, s.size());
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *s.begin());
    s.erase(s.begin(), s.end());
    EXPECT_EQ(0u, s.size());
    s.push_back(gapped<dna4>{dna4::A});
    EXPECT_EQ(1u, s.size());
    EXPECT_EQ(gapped<dna4>{dna4::A}, *s.begin());
    // case 2 i): into non-empty sequence, old state: only alphabet symbols
    std::initializer_list<gapped<dna4>> l1{dna4::T, dna4::A};
    sequence_t t{l1};
    t.push_back(dna4::C);
    EXPECT_EQ(3u, t.size());
    EXPECT_EQ(gapped<dna4>{dna4::C}, *(t.end()-1));
    t.erase(t.end()-1);
    EXPECT_EQ(2u, t.size());
    t.push_back(gap::GAP);
    EXPECT_EQ(3u, t.size());
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *(t.end()-1));
    // case 2 ii): into non-empty sequence, old state: only gap symbols
    std::initializer_list<gapped<dna4>> l2{gap::GAP, gap::GAP};
    sequence_t u{l2};
    u.push_back(dna4::A);
    EXPECT_EQ(3u, u.size());
    EXPECT_EQ(gapped<dna4>{dna4::A}, *(u.end()-1));
    u.erase(u.end()-1);
    EXPECT_EQ(2u, u.size());
    u.push_back(gap::GAP);
    EXPECT_EQ(3u, u.size());
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *(u.end()-1));

    // case 2 iii): into non-empty sequence, old state: alphabet and gap symbol
    std::initializer_list<gapped<dna4>> l3{gap::GAP, dna4::T, gap::GAP, dna4::A};
    sequence_t v{l3};
    v.push_back(dna4::C);
    EXPECT_EQ(5u, v.size());
    EXPECT_EQ(gapped<dna4>{dna4::C}, *(v.end()-1));
    v.erase(v.end()-1);
    EXPECT_EQ(4u, v.size());
    v.push_back(gap::GAP);
    EXPECT_EQ(5u, v.size());
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *(v.end()-1));
}

TEST(aligned_sequences_test, sequence_concepts_pop_back)
{
    // case 1: pop last element from sequence made out of gaps only
    std::initializer_list<gapped<dna4>> l1{gap::GAP, gap::GAP, gap::GAP};
    sequence_t s{l1};
    s.pop_back();
    EXPECT_EQ(2u, s.size());
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *(s.end()-1));

    // case 2: pop last element from sequence made out of alphabet symbols only
    std::initializer_list<gapped<dna4>> l2{dna4::A, dna4::C, dna4::G};
    sequence_t t{l2};
    t.pop_back();
    EXPECT_EQ(2u, t.size());
    EXPECT_EQ(gapped<dna4>{dna4::C}, *(t.end()-1));

    // case 3: pop last element from mixed sequence
    std::initializer_list<gapped<dna4>> l3{gap::GAP, dna4::T, gap::GAP, dna4::A};
    sequence_t u{l3};
    u.pop_back();
    EXPECT_EQ(3u, u.size());
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *(u.end()-1));

    // case 4: throw assertion when trying to pop_back on empty sequence
    sequence_t v{};
    // TODO: test assertion
    //ASSERT_EXIT(v.pop_back(), testing::ExitedWithCode(6), "");
}

TEST(aligned_sequences_test, sequence_concepts_clear)
{
    // clear empty sequence
    sequence_t s{};
    s.clear();
    EXPECT_EQ(0u, s.size());

    // clear non-empty sequence
    sequence_t t{gap::GAP, dna4::T, gap::GAP, gap::GAP};
    t.clear();
    EXPECT_EQ(0u, t.size());
}

TEST(aligned_sequences_test, sequence_concepts_front)
{
    // case 1: front on empty sequence
    sequence_t s{};

    // case 2: front element is gap
    sequence_t t{gap::GAP, dna4::T};
    EXPECT_EQ(gapped<dna4>{gap::GAP}, t.front());

    // case 3: front element is alphabet symbol
    sequence_t u{dna4::A, dna4::T};
    EXPECT_EQ(gapped<dna4>{dna4::A}, u.front());
}

TEST(aligned_sequences_test, sequence_concepts_back)
{
    sequence_t s{dna4::A};
    EXPECT_EQ(gapped<dna4>{dna4::A}, s.back());
    s.push_back(gap::GAP);
    EXPECT_EQ(gapped<dna4>{gap::GAP}, s.back());
}

TEST(aligned_sequences_test, get_underlying_sequence)
{
    sequence_t s{};
    EXPECT_EQ(0u, s.get_underlying_sequence().size());
    s.push_back(gap::GAP);
    EXPECT_EQ(0u, s.get_underlying_sequence().size());
    s.push_back(dna4::A);
    EXPECT_EQ(1u, s.get_underlying_sequence().size());
    EXPECT_EQ(gapped<dna4>{dna4::A}, s.get_underlying_sequence().at(0));
}

TEST(aligned_sequences_test, insert_gap)
{
    // case 1: into empty sequence
    sequence_t s{};
    s.insert_gap(0);
    EXPECT_EQ(1u, s.size());
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *(s.begin()));

    // case 2: insert into non-empty sequence
    s.push_back(dna4::A);
    EXPECT_EQ(2u, s.size());
    s.insert_gap(2, 2);
    EXPECT_EQ(4u, s.size());
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *(s.end()-1));

}

TEST(aligned_sequences_test, map_to_aligned_position)
{
    sequence_t s{gap::GAP, dna4::A, gap::GAP, dna4::T};
    EXPECT_EQ(1u, s.map_to_aligned_position(0));
    EXPECT_EQ(3u, s.map_to_aligned_position(1));
}

TEST(aligned_sequences_test, map_to_underlying_position)
{
    sequence_t s{gap::GAP, dna4::A, gap::GAP, gap::GAP, dna4::T};
    std::array<sequence_t::difference_type, 5> upos{ {-1, 0, 0, 0, 1} };
    auto it = upos.begin();
    for (sequence_t::size_type   i = 0; i < s.size(); ++i)
        EXPECT_EQ(*(it++), s.map_to_underlying_position(i));
}

TEST(aligned_sequences_test, random_access_operators)
{
    // case 1: []-operator on gap postion
    sequence_t s{gap::GAP, dna4::A, gap::GAP, gap::GAP, dna4::T};
    EXPECT_EQ(gapped<dna4>{gap::GAP}, s[0]);
    EXPECT_EQ(gapped<dna4>{dna4::T}, s[4]);
    // case 2: at()-operator
    EXPECT_EQ(gapped<dna4>{gap::GAP}, s.at(0));
    EXPECT_EQ(gapped<dna4>{dna4::T}, s.at(4));
}
