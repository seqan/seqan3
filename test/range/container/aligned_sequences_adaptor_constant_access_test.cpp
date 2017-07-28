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

#include <seqan3/alphabet/gap/gapped_alphabet.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/container/aligned_sequence_adaptor_constant_access.hpp>

using namespace seqan3;

//using gapped_alphabet = gapped_alphabet<dna4>;

TEST(aligned_sequences_test, constructor_empty)
{
    aligned_sequence_adaptor_constant_access<std::vector<gapped_alphabet<dna4>>> s;
}

/*
TEST(aligned_sequences_test, constructor_sequence)
{
    std::vector<dna4> seq = {dna4::A};
    aligned_sequence_adaptor_constant_access<std::vector<dna4>> s(seq);
}


// copy constructors
TEST(aligned_sequences_test, constructor_cpy)
{
    // copy construction with empty sequence
    aligned_sequence_adaptor_constant_access<std::vector<dna4>> as_base;
    aligned_sequence_adaptor_constant_access<std::vector<dna4>> as_derived(as_base);
    // copy construction with non-empty sequence
    std::vector<dna4> seq = {dna4::A, dna4::C, dna4::G, dna4::T};
    aligned_sequence_adaptor_constant_access<std::vector<dna4>> as_base2(seq);
    aligned_sequence_adaptor_constant_access<std::vector<dna4>> as_derived2(as_base2);
    // copy construction via assignment
    aligned_sequence_adaptor_constant_access<std::vector<dna4>> as_derived3 = as_base2;
}

// move constructors
TEST(aligned_sequences_test, constructor_move)
{
    // move construction with empty sequence
    aligned_sequence_adaptor_constant_access<std::vector<dna4>> as_base;
    aligned_sequence_adaptor_constant_access<std::vector<dna4>> as_derived(std::move(as_base));
    // move construction with non-empty sequence
    std::vector<dna4> seq = {dna4::A, dna4::C, dna4::G, dna4::T};
    std::vector<dna4> seq2 = {dna4::T, dna4::A};
    aligned_sequence_adaptor_constant_access<std::vector<dna4>> as_base2(seq);
    aligned_sequence_adaptor_constant_access<std::vector<dna4>> as_derived2(std::move(as_base2));
    auto it2 = as_derived2.begin();
    EXPECT_EQ(dna4::A, *it2++);
    EXPECT_EQ(dna4::C, *it2++);
    EXPECT_EQ(dna4::G, *it2++);
    EXPECT_EQ(dna4::T, *it2);
    // move construction via assignment
    aligned_sequence_adaptor_constant_access<std::vector<dna4>> as_base3(seq), as_derived3(seq2);
    as_derived3 = std::move(as_base3);
    auto it3 = as_derived3.begin();
    EXPECT_EQ(dna4::A, *it3++);
    EXPECT_EQ(dna4::C, *it3++);
    EXPECT_EQ(dna4::G, *it3++);
    EXPECT_EQ(dna4::T, *it3);
}

// Destructor
TEST(aligned_sequences_test, destructor)
{
    std::vector<dna4> seq = {dna4::T, dna4::A};
    using aligned_sequence_t = aligned_sequence_adaptor_constant_access<std::vector<dna4>>;
    aligned_sequence_t s(seq);
    aligned_sequence_t * s_ptr = &s;
    s_ptr->aligned_sequence_t::~aligned_sequence_t();
}


// Container concept functions
TEST(aligned_sequences_test, container_concepts)
{
    std::vector<dna4> seq = {dna4::T, dna4::A}, seq2 = {dna4::C};
    aligned_sequence_adaptor_constant_access<std::vector<dna4>> s(seq);

    // begin and end iterators
    auto iterator = s.begin();
    EXPECT_EQ(dna4::T, *iterator);
    iterator = s.end();
    EXPECT_EQ(dna4::A, *(--iterator));

    // const begin and const end iterators
    auto iterator_const = s.cbegin();
    EXPECT_EQ(dna4::T, *iterator_const);
    iterator_const = s.cend();
    EXPECT_EQ(dna4::A, *(--iterator_const));

    // boolean operators == and !=
    aligned_sequence_adaptor_constant_access<std::vector<dna4>> t(seq);
    EXPECT_TRUE(s == t);
    aligned_sequence_adaptor_constant_access<std::vector<dna4>> u(seq2);
    EXPECT_TRUE(t != u);

    // swap, size, max_size, empty
    t.swap(u);
    iterator = t.begin();
    EXPECT_EQ(dna4::C, *iterator);
    iterator = u.begin();
    EXPECT_EQ(dna4::T, *iterator++);
    EXPECT_EQ(dna4::A, *iterator);
    swap(t, u); // friend
    iterator = t.begin();
    EXPECT_EQ(dna4::T, *iterator++);
    EXPECT_EQ(dna4::A, *iterator);
    iterator = u.begin();
    EXPECT_EQ(dna4::C, *iterator);

    aligned_sequence_adaptor_constant_access<std::vector<dna4>>::size_type max_size = s.max_size();

    EXPECT_FALSE(s.empty());
    aligned_sequence_adaptor_constant_access<std::vector<dna4>> s_empty;
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
