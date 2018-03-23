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
#include <seqan3/range/container/gap_decorator_anchor_list.hpp>
#include <seqan3/range/detail/random_access_iterator.hpp>

using namespace seqan3;
using container_t = std::vector<dna4>;
using sequence_t = gap_decorator_anchor_list<container_t>;

/*
// default (cp, mv) (de)constructors
TEST(gap_anchor_list_test, constructors)
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
TEST(gap_anchor_list_test, constructor_by_sequence)
{
    // case 1: empty sequence
    container_t sequence_empty{};
    sequence_t as_empty(&sequence_empty);
    EXPECT_EQ(0u, as_empty.size());
    // case 2: non-empty sequence
    container_t sequence{dna4::A, dna4::C, dna4::G, dna4::T};
    sequence_t as(&sequence);
    auto it = as.begin();
    EXPECT_EQ(sequence.size(), as.size());
    EXPECT_EQ(gapped<dna4>{dna4::A}, it[0]);
    EXPECT_EQ(gapped<dna4>{dna4::T}, it[3]);
}

// move constructors
TEST(gap_anchor_list_test, constructor_move)
{
    // case 1.1: move construction with empty sequence
    sequence_t as_base;
    sequence_t as_derived(std::move(as_base));
    // case 1.2: move construction with non-empty sequence
    container_t seq = {dna4::A, dna4::C, dna4::G, dna4::T};
    sequence_t as_base2(&seq);
    sequence_t as_derived2(std::move(as_base2));
    auto it2 = as_derived2.begin();
    EXPECT_EQ(gapped<dna4>{dna4::A}, *it2++);
    EXPECT_EQ(gapped<dna4>{dna4::C}, *it2++);
    EXPECT_EQ(gapped<dna4>{dna4::G}, *it2++);
    EXPECT_EQ(gapped<dna4>{dna4::T}, *it2);

    // case 2:  move construction via assignment
    container_t seq2 = {dna4::T, dna4::A};
    sequence_t as_base3(&seq2), as_derived3(&seq);
    as_derived3 = std::move(as_base3);
    auto it3 = as_derived3.begin();
    auto it3_end = as_derived3.end();
    EXPECT_EQ(gapped<dna4>{dna4::T}, *it3++);
    EXPECT_EQ(gapped<dna4>{dna4::A}, *it3++);
    EXPECT_EQ(it3, it3_end);
}

// Destructor
TEST(gap_anchor_list_test, destructor)
{
    container_t seq = {dna4::T, dna4::A};
    using aligned_sequence_t = sequence_t;
    aligned_sequence_t s(&seq);
    aligned_sequence_t * s_ptr = new aligned_sequence_t(&seq);//&s;
    s_ptr->gap_decorator_anchor_list<container_t>::~gap_decorator_anchor_list<container_t>();
}

// Container concept functions begin(), end()
TEST(gap_anchor_list_test, container_concepts_iterators)
{
    container_t seq = {dna4::T, dna4::A};
    sequence_t s(&seq);
    auto it = s.begin();
    EXPECT_EQ(gapped<dna4>{dna4::T}, *it);
    it = s.end();
    EXPECT_EQ(gapped<dna4>{dna4::A}, *(--it));
}

TEST(gap_anchor_list_test, container_concepts_boolean)
{
    container_t seq = {dna4::T, dna4::A}, seq2 = {dna4::C};
    sequence_t s(&seq), t(&seq);
    EXPECT_TRUE(s == t);
    sequence_t u(&seq2);
    EXPECT_TRUE(t != u);
}

// swap, size, max_size, empty
TEST(gap_anchor_list_test, container_concepts_swap)
{
    container_t seq = {dna4::T, dna4::A}, seq2 = {dna4::C};
    sequence_t t(&seq), u(&seq2);
    t.swap(u);
    auto it = t.begin();
    EXPECT_EQ(gapped<dna4>{dna4::C}, *it);
    it = u.begin();
    EXPECT_EQ(gapped<dna4>{dna4::T}, *it++);
    EXPECT_EQ(gapped<dna4>{dna4::A}, *it);
    swap(t, u); // friend
    it = t.begin();
    EXPECT_EQ(gapped<dna4>{dna4::T}, *it++);
    EXPECT_EQ(gapped<dna4>{dna4::A}, *it);
    it = u.begin();
    EXPECT_EQ(gapped<dna4>{dna4::C}, *it);

    sequence_t::size_type max_size = t.max_size();
    EXPECT_TRUE(max_size > 0);
    EXPECT_FALSE(t.empty());
    container_t s_empty{};
    sequence_t seq_empty(&s_empty);
    EXPECT_TRUE(seq_empty.empty());
}
*/
// gap insertion
TEST(gap_anchor_list_test, insert_gap)
{
    // case 1.1: insert 1 gap by position iterator into empty sequence
    sequence_t s{};
    auto it = s.insert_gap(s.begin());  //auto it = s.insert_gap(s.begin());
    EXPECT_EQ(1u, s.size());
    EXPECT_EQ(s.begin(), it);
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *s.begin());

    // case 1.2: insert by position iterator into non-empty sequence
    container_t seq{dna4::C, dna4::T, dna4::A};
    container_t l{seq};
    std::vector<gapped<dna4>> l_post{gap::GAP, dna4::C, gap::GAP, gap::GAP, dna4::T, dna4::A};
    sequence_t t(&seq);
    it = t.insert_gap(t.begin());
    for (unsigned int i = 0; i < t.size(); ++i)
        std::cout << "t[" << i << "] = " << t[i] << std::endl;
    EXPECT_EQ(t.begin(), it);
    EXPECT_EQ(4u, t.size());
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *it++);
    EXPECT_EQ(gapped<dna4>{dna4::C}, *it++);
    EXPECT_EQ(gapped<dna4>{dna4::T}, *it++);
    EXPECT_EQ(gapped<dna4>{dna4::A}, *it);

/*    t.insert_gap(t.begin()+2);
    it = t.insert_gap(t.begin()+2);
    EXPECT_EQ(t.begin() + 2, it);
    EXPECT_EQ(6u, t.size());
    [[maybe_unused]] auto l_it = l_post.begin();
    for (auto it = t.begin(); l_it != l_post.end(); ++l_it, ++it)
        EXPECT_EQ(*l_it, *it);

    // case 2: insert after last element
    sequence_t w(seq);
    it = w.insert_gap(w.end());
    EXPECT_EQ(4u, w.size());
    EXPECT_EQ(w.end()-1, it);

    // case 3: insert by position index
    // case 3.1: into empty sequence
    sequence_t z{};
    bool b = z.insert_gap(1);
    EXPECT_EQ(false, b);
    EXPECT_EQ(0u, z.size());
    b = z.insert_gap(0);
    EXPECT_EQ(true, b);
    EXPECT_EQ(1u, z.size());
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *z.begin());
    // case 3.2: into non-empty sequence
    container_t seq_u{dna4::A, dna4::C, dna4::G};
    sequence_t u{&seq_u};
    b = u.insert_gap(1);
    EXPECT_EQ(true, b);
    EXPECT_EQ(4u, u.size());
    it = u.begin();
    EXPECT_EQ(gapped<dna4>{dna4::A}, *it++);
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *it++);
    EXPECT_EQ(gapped<dna4>{dna4::C}, *it++);
    EXPECT_EQ(gapped<dna4>{dna4::G}, *it++);
    // case 3.3: into non-empty sequence past-the-end
    b = u.insert_gap(4u);
    EXPECT_EQ(true, b);
    EXPECT_EQ(5u, u.size());
    EXPECT_EQ(gapped<dna4>{dna4::G}, *(u.begin()+3));
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *(u.end()-1));
    */
}

/*
// insertion of multiple gaps
TEST(gap_anchor_list_test, insert_gaps)
{
    // case 1.1: insert 2 gaps by position iterator into empty sequence
    sequence_t s{};
    auto it = s.insert_gap(s.begin(), 2);
    EXPECT_EQ(2u, s.size());
    EXPECT_EQ(s.begin(), it);
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *s.begin());
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *(s.end()-1));

    // case 1.2: insert by position iterator into non-empty sequence
    container_t seq{dna4::C, dna4::T, dna4::A};
    container_t l{seq};
    std::vector<gapped<dna4>> l_post{gap::GAP, gap::GAP, dna4::C, gap::GAP, gap::GAP, dna4::T, dna4::A};
    sequence_t t(&seq);
    it = t.insert_gap(t.begin(), 2);
    EXPECT_EQ(t.begin(), it);
    EXPECT_EQ(5u, t.size());
    bool b = t.insert_gap(3, 2);
    EXPECT_EQ(true, b);
    EXPECT_EQ(7u, t.size());
    auto l_it = l_post.begin();
    for (auto it = t.begin(); l_it != l_post.end(); ++l_it, ++it)
        EXPECT_EQ(*l_it, *it);

    // case 2: insert after last element
    sequence_t w(seq);
    it = w.insert_gap(w.end(), 3);
    EXPECT_EQ(6u, w.size());
    EXPECT_EQ(w.begin()+3, it);
    EXPECT_EQ(gapped<dna4>{dna4::A}, *(w.begin()+2));
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *(w.end()-1));

    // case 3: insert by position index
    // case 3.1: into empty sequence
    sequence_t z{};
    b = z.insert_gap(1, 2);
    EXPECT_EQ(false, b);
    EXPECT_EQ(0u, z.size());
    b = z.insert_gap(0, 128u);
    EXPECT_EQ(true, b);
    EXPECT_EQ(128u, z.size());
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *z.begin());
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *(z.end()-1));

    // case 3.2: into non-empty sequence
    container_t u_vec{dna4::A, dna4::C, dna4::G};
    sequence_t u{&u_vec};
    b = u.insert_gap(1, 16u);
    EXPECT_EQ(true, b);
    EXPECT_EQ(19u, u.size());
    it = u.begin();
    EXPECT_EQ(gapped<dna4>{dna4::A}, *it);
    EXPECT_EQ(gapped<dna4>{gap::GAP}, it[1]);
    EXPECT_EQ(gapped<dna4>{gap::GAP}, it[16]);
    EXPECT_EQ(gapped<dna4>{dna4::C}, it[17]);
    // case 3.3: into non-empty sequence past-the-end
    b = u.insert_gap(19u, 3);
    EXPECT_EQ(true, b);
    EXPECT_EQ(22u, u.size());
    EXPECT_EQ(gapped<dna4>{dna4::G}, it[18]);
    EXPECT_EQ(gapped<dna4>{gap::GAP}, it[19]);
}

// erase single elements given by single iterator or position index
TEST(gap_anchor_list_test, erase_gap)
{
    container_t seq{dna4::T, dna4::A};
    std::vector<gapped<dna4>> aseq{gap::GAP, dna4::T, dna4::A};
    sequence_t s{&seq};
    s.insert_gap(0, 2);
    auto it = s.begin();
    EXPECT_EQ(4u, s.size());
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *it++);
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *it++);
    EXPECT_EQ(gapped<dna4>{dna4::T}, *it++);
    EXPECT_EQ(gapped<dna4>{dna4::A}, *it);

    // case 1: erase 1st element by iterator
    // case 1.1: erase 1st gap
    it = s.erase_gap(s.begin());

    EXPECT_EQ(it, s.begin());
    EXPECT_EQ(3u, s.size());
    auto it_aseq = aseq.begin();
    for (; it != s.end(); ++it, ++it_aseq)
        EXPECT_EQ(*it_aseq, *it);
    // case 1.2: erase midth gap
    s.insert_gap(s.begin()+2, 1);   // -T-A
    EXPECT_EQ(4u, s.size());
    it = s.erase_gap(s.begin()+2);        // -TA
    EXPECT_EQ(3u, s.size());
    EXPECT_EQ(s.begin()+2, it);
    it = s.begin();
    EXPECT_EQ(gapped<dna4>{gap::GAP}, it[0]);
    EXPECT_EQ(gapped<dna4>{dna4::T}, it[1]);
    EXPECT_EQ(gapped<dna4>{dna4::A}, it[2]);
    // case 1.3: erase last element
    s.insert_gap(3u, 1);   // -TA-
    EXPECT_EQ(4u, s.size());
    it = s.erase_gap(s.end()-1); // -TA
    EXPECT_EQ(3u, s.size());
    EXPECT_EQ(it, s.end());
    EXPECT_EQ(gapped<dna4>{dna4::A}, *(it-1));
    // case 2: erase single gaps by position index
    // case 2.1: erase 1st gap
    container_t vec_t{dna4::A, dna4::C};
    sequence_t t{&vec_t};
    t.insert_gap(0);            // -AC
    bool b = t.erase_gap(0u);    // AC
    EXPECT_EQ(2u, t.size());
    EXPECT_EQ(true, b);
    it = t.begin();
    EXPECT_EQ(gapped<dna4>{dna4::A}, it[0]);
    EXPECT_EQ(gapped<dna4>{dna4::C}, it[1]);
    // case 2.2: erase midth gap or try beyond current size
    t.insert_gap(1);            // A-C
    EXPECT_EQ(3u, t.size());
    b = t.erase_gap(3u);
    EXPECT_FALSE(b);
    b = t.erase_gap(1);    // AC
    EXPECT_EQ(2u, t.size());
    EXPECT_EQ(gapped<dna4>{dna4::A}, it[0]);
    EXPECT_EQ(gapped<dna4>{dna4::C}, it[1]);
    // case 2.3: erase last gap
    t.insert_gap(2u);       // AC-
    EXPECT_EQ(3u, t.size());
    b = t.erase_gap(2u);    // AC
    EXPECT_EQ(gapped<dna4>{dna4::A}, it[0]);
    EXPECT_EQ(gapped<dna4>{dna4::C}, it[1]);
}

// erase elements in range between two iterators or position indices
TEST(gap_anchor_list_test, erase_gaps_iterators)
{
    std::vector<gapped<dna4>> gseq{gap::GAP, gap::GAP, dna4::T, gap::GAP, gap::GAP, dna4::A, gap::GAP, gap::GAP};
    container_t vec_s{dna4::T, dna4::A};
    sequence_t s{&vec_s};    // TA
    s.insert_gap(0u, 2);    // --TA
    s.insert_gap(3u, 2);    // --T--A
    s.insert_gap(6u, 2);    // --T--A--
    EXPECT_EQ(8u, s.size());

    // case 1: erase range given by iterators
    // case 1.1: erase front range gaps
    auto it = s.erase_gap(s.begin(), s.begin() + 3);    // T--A--
    EXPECT_EQ(it, s.begin());
    EXPECT_EQ(6u, s.size());
    auto it_l = gseq.begin()+2;
    for (it = s.begin(); it_l != gseq.end(); ++it, ++it_l)
        EXPECT_EQ(*it_l, *it);

    // case 1.2: erase mid range gaps
    it = s.erase_gap(s.begin()+1, s.end()-1);   // TA-
    EXPECT_EQ(3u, s.size());
    EXPECT_EQ(it, s.begin()+1);
    EXPECT_EQ(gapped<dna4>{dna4::A}, it[0]);
    EXPECT_EQ(gapped<dna4>{gap::GAP}, it[1]);

    // case 1.3: erase last gaps
    it = s.erase_gap(s.begin()+2, s.end());     // TA
    EXPECT_EQ(2u, s.size());
    EXPECT_EQ(it, s.end());
    EXPECT_EQ(gapped<dna4>{dna4::T}, *s.begin());
    EXPECT_EQ(gapped<dna4>{dna4::A}, *(s.begin()+1));
}

// case 2: erase range by position indices
TEST(gap_anchor_list_test, erase_gaps_indices)
{
    std::vector<gapped<dna4>> gseq{gap::GAP, gap::GAP, dna4::T, gap::GAP, gap::GAP, dna4::A, gap::GAP, gap::GAP};
    container_t vec_s{dna4::T, dna4::A};
    sequence_t s{&vec_s};
    s.insert_gap(0u, 2);
    s.insert_gap(3u, 2);
    s.insert_gap(6u, 2);
    EXPECT_EQ(8u, s.size());

    // case 1: erase range given by position indices
    // case 1.1: erase front range gaps
    bool b = s.erase_gap(0u, 3u);    // T--A--
    EXPECT_EQ(true, b);
    EXPECT_EQ(6u, s.size());
    auto it_l = gseq.begin()+2;
    for (auto it = s.begin(); it_l != gseq.end(); ++it, ++it_l)
        EXPECT_EQ(*it_l, *it);

    // case 1.2: erase mid range gaps
    b = s.erase_gap(1u, 4u);   // TA--
    EXPECT_EQ(4u, s.size());
    EXPECT_EQ(true, b);
    auto it = s.begin();
    EXPECT_EQ(gapped<dna4>{dna4::A}, it[1]);
    EXPECT_EQ(gapped<dna4>{gap::GAP}, it[2]);

    // case 1.3: erase last gaps
    b = s.erase_gap(2u, 4u);     // TA
    EXPECT_EQ(2u, s.size());
    EXPECT_EQ(true, b);
    EXPECT_EQ(gapped<dna4>{dna4::T}, it[0]);
    EXPECT_EQ(gapped<dna4>{dna4::A}, it[1]);
}

TEST(gap_anchor_list_test, sequence_concepts_push_back)
{
    // case 1: into empty sequence a gap or alphabet symbol
    container_t vec{};
    sequence_t s{&vec};
    s.push_back();
    EXPECT_EQ(1u, s.size());
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *s.begin());

    // case 2: into non-empty sequence
    container_t t_vec{dna4::A};
    sequence_t t{&t_vec};
    t.push_back();
    EXPECT_EQ(2u, t.size());
    EXPECT_EQ(gapped<dna4>{dna4::A}, *(t.begin()));
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *(t.begin()+1));
}

TEST(gap_anchor_list_test, sequence_concepts_pop_back)
{
    // case 1: pop last element from sequence made out of gaps only
    container_t vec{};
    sequence_t s{&vec};
    s.insert_gap(0, 4);
    EXPECT_EQ(4u, s.size());
    bool b = s.pop_back();
    EXPECT_EQ(3u, s.size());
    EXPECT_EQ(gapped<dna4>{gap::GAP}, *(s.end()-1));

    // case 2: pop last element from sequence with gaps
    container_t vec_t{dna4::A, dna4::A};
    sequence_t t{&vec_t};
    t.insert_gap(0u);   // -AA
    b = t.pop_back();   // -AA
    EXPECT_EQ(3u, t.size());
    EXPECT_FALSE(b);
    t.push_back();      // -AA-
    EXPECT_EQ(4u, t.size());
    b = t.pop_back();
    EXPECT_TRUE(b);
    EXPECT_EQ(3u, t.size());
}

TEST(gap_anchor_list_test, sequence_concepts_clear)
{
    // case 1.1: clear empty sequence
    container_t vec{};
    sequence_t s{&vec};
    s.clear();
    EXPECT_EQ(0u, s.size());

    // case 1.2: clear non-empty sequence with no gaps
    container_t vec_t{dna4::T};
    sequence_t t{&vec_t};
    t.clear();
    EXPECT_EQ(1u, t.size());

    // case 1.3: clear non-empty sequence with gaps
    t.insert_gap(0, 2); // --T
    t.insert_gap(3, 2); // --T--
    t.clear();
    EXPECT_EQ(1u, t.size());
    EXPECT_EQ(gapped<dna4>{dna4::T}, *t.begin());
}


TEST(gap_anchor_list_test, sequence_concepts_front)
{
    // case 1: front on empty sequence
    container_t v{};
    sequence_t s{&v};
    //should throw assertion s.front();

    // case 2: front element is alphabet symbol
    container_t seq_u{dna4::A, dna4::T};
    sequence_t u{&seq_u};
    EXPECT_EQ(gapped<dna4>{dna4::A}, u.front());

    // case 3: front element is gap
    u.insert_gap(0);
    EXPECT_EQ(gapped<dna4>{gap::GAP}, u.front());
}

TEST(gap_anchor_list_test, sequence_concepts_back)
{
    container_t v{dna4::A};
    sequence_t s{&v};
    EXPECT_EQ(gapped<dna4>{dna4::A}, s.back());
    s.push_back();
    EXPECT_EQ(gapped<dna4>{gap::GAP}, s.back());
}

TEST(gap_anchor_list_test, get_underlying_sequence)
{
    container_t v;
    sequence_t s{&v};
    EXPECT_EQ(0u, s.get_underlying_sequence().size());
    s.push_back();
    EXPECT_EQ(0u, s.get_underlying_sequence().size());
    container_t seq{dna4::A};
    s.set_underlying_sequence(seq);
    EXPECT_EQ(1u, s.get_underlying_sequence().size());
    EXPECT_EQ(gapped<dna4>{dna4::A}, s.get_underlying_sequence().at(0));
}

TEST(gap_anchor_list_test, set_underlying_sequence)
{
    // case 1: set empty sequence
    container_t v{};
    sequence_t s{v};
    EXPECT_EQ(0u, s.get_underlying_sequence().size());
    // case 2: set non-empty sequence
    container_t seq{dna4::A, dna4::A, dna4::A, dna4::A};
    s.set_underlying_sequence(seq);
    EXPECT_EQ(seq.size(), s.get_underlying_sequence().size());
    EXPECT_EQ(gapped<dna4>{dna4::A}, s.get_underlying_sequence().at(0));
}

TEST(gap_anchor_list_test, map_to_aligned_position)
{
    container_t s_vec{dna4::A, dna4::T};
    sequence_t s{s_vec};
    s.insert_gap(0); s.insert_gap(2); // -A-T
    EXPECT_EQ(1u, s.map_to_aligned_position(0));
    EXPECT_EQ(3u, s.map_to_aligned_position(1));
}

TEST(gap_anchor_list_test, map_to_underlying_position)
{
    container_t s_vec{dna4::A, dna4::T};
    sequence_t s{s_vec};
    s.insert_gap(0); s.insert_gap(2, 2); // -A--T
    std::array<sequence_t::difference_type, 5> upos{ {-1, 0, 0, 0, 1} };
    auto it = upos.begin();
    for (sequence_t::size_type   i = 0; i < s.size(); ++i)
        EXPECT_EQ(*(it++), s.map_to_underlying_position(i));
}

TEST(gap_anchor_list_test, random_access_operators)
{
    // []-operator on gap postion
    container_t s_vec{dna4::A, dna4::T};
    sequence_t s{s_vec};
    s.insert_gap(0); s.insert_gap(2, 2);
    EXPECT_EQ(gapped<dna4>{gap::GAP}, s[0]);
    EXPECT_EQ(gapped<dna4>{dna4::T}, s[4]);
    // at()-operator
    EXPECT_EQ(gapped<dna4>{gap::GAP}, s.at(0));
    EXPECT_EQ(gapped<dna4>{dna4::T}, s.at(4));
}
*/
