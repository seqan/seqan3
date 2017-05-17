// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2017, Knut Reinert, FU Berlin
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
// ==========================================================================

#include <gtest/gtest.h>

#include <variant>

#include <seqan3/alphabet/gap/gap.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/union_alphabet.hpp>

using namespace seqan3;

TEST(union_alphabet_test, default_constructor)
{
    using alphabet_t = union_alphabet<dna4, gap>;
    constexpr alphabet_t letter1{};

    EXPECT_EQ(letter1._value, 0);
}

TEST(union_alphabet_test, initialize_from_component_alphabet)
{
    using alphabet_t = union_alphabet<dna4, dna5, gap>;
    using variant_t = std::variant<dna4, dna5, gap>;

    constexpr variant_t variant0{dna4::A};
    constexpr alphabet_t letter0{dna4::A};

    constexpr variant_t variant1 = dna4::C;
    constexpr alphabet_t letter1 = dna4::C;

    constexpr variant_t variant2 = {dna4::G};
    constexpr alphabet_t letter2 = {dna4::G};

    constexpr variant_t variant3 = static_cast<variant_t>(dna4::T);
    constexpr alphabet_t letter3 = static_cast<alphabet_t>(dna4::T);

    constexpr variant_t variant4 = {static_cast<variant_t>(dna5::A)};
    constexpr alphabet_t letter4 = {static_cast<alphabet_t>(dna5::A)};

    variant_t variant5{dna5::C};
    alphabet_t letter5{dna5::C};

    variant_t variant6 = dna5::G;
    alphabet_t letter6 = dna5::G;

    variant_t variant7 = {dna5::T};
    alphabet_t letter7 = {dna5::T};

    variant_t variant8 = static_cast<variant_t>(dna5::N);
    alphabet_t letter8 = static_cast<alphabet_t>(dna5::N);

    variant_t variant9 = {static_cast<variant_t>(gap::GAP)};
    alphabet_t letter9 = {static_cast<alphabet_t>(gap::GAP)};

    EXPECT_EQ(letter0.to_rank(), 0);
    EXPECT_EQ(letter1.to_rank(), 1);
    EXPECT_EQ(letter2.to_rank(), 2);
    EXPECT_EQ(letter3.to_rank(), 3);
    EXPECT_EQ(letter4.to_rank(), 4);
    EXPECT_EQ(letter5.to_rank(), 5);
    EXPECT_EQ(letter6.to_rank(), 6);
    EXPECT_EQ(letter7.to_rank(), 7);
    EXPECT_EQ(letter8.to_rank(), 8);
    EXPECT_EQ(letter9.to_rank(), 9);
}

TEST(union_alphabet_test, initialise_from_same_component_alphabet)
{
    using alphabet_t = union_alphabet<dna4, dna4>;

    constexpr alphabet_t letter0{std::in_place_index_t<0>{}, dna4::A};
    constexpr alphabet_t letter1{std::in_place_index_t<0>{}, dna4::C};
    constexpr alphabet_t letter2{std::in_place_index_t<0>{}, dna4::G};
    constexpr alphabet_t letter3{std::in_place_index_t<0>{}, dna4::T};
    constexpr alphabet_t letter4{std::in_place_index_t<1>{}, dna4::A};
    constexpr alphabet_t letter5{std::in_place_index_t<1>{}, dna4::C};
    constexpr alphabet_t letter6{std::in_place_index_t<1>{}, dna4::G};
    constexpr alphabet_t letter7{std::in_place_index_t<1>{}, dna4::T};

    EXPECT_EQ(letter0.to_rank(), 0);
    EXPECT_EQ(letter1.to_rank(), 1);
    EXPECT_EQ(letter2.to_rank(), 2);
    EXPECT_EQ(letter3.to_rank(), 3);
    EXPECT_EQ(letter4.to_rank(), 4);
    EXPECT_EQ(letter5.to_rank(), 5);
    EXPECT_EQ(letter6.to_rank(), 6);
    EXPECT_EQ(letter7.to_rank(), 7);
}

TEST(union_alphabet_test, copy_constructor)
{
    using alphabet_t = union_alphabet<dna4, gap>;
    constexpr alphabet_t letter1{dna4::T};
    alphabet_t letter2{letter1};

    EXPECT_EQ(letter1._value, 3);
    EXPECT_EQ(letter2._value, 3);
}

TEST(union_alphabet_test, move_constructor)
{
    using alphabet_t = union_alphabet<dna4, gap>;
    alphabet_t letter1{dna4::G};
    alphabet_t letter2{std::move(letter1)};

    EXPECT_EQ(letter2._value, 2);
}

TEST(union_alphabet_test, assign_from_component_alphabet)
{
    using alphabet_t = union_alphabet<dna4, dna5, gap>;
    using variant_t = std::variant<dna4, dna5, gap>;
    alphabet_t letter{};
    variant_t variant{};

    variant = dna4::A;
    letter = dna4::A;
    EXPECT_EQ(variant.index(), 0);
    EXPECT_EQ(letter.to_rank(), 0);

    variant = {dna4::C};
    letter = {dna4::C};
    EXPECT_EQ(variant.index(), 0);
    EXPECT_EQ(letter.to_rank(), 1);

    variant = static_cast<variant_t>(dna4::G);
    letter = static_cast<alphabet_t>(dna4::G);
    EXPECT_EQ(variant.index(), 0);
    EXPECT_EQ(letter.to_rank(), 2);

    variant = {static_cast<variant_t>(dna4::T)};
    letter = {static_cast<alphabet_t>(dna4::T)};
    EXPECT_EQ(letter.to_rank(), 3);

    letter = dna5::A;
    EXPECT_EQ(letter.to_rank(), 4);

    letter = dna5::C;
    EXPECT_EQ(letter.to_rank(), 5);

    letter = dna5::G;
    EXPECT_EQ(letter.to_rank(), 6);

    letter = dna5::T;
    EXPECT_EQ(letter.to_rank(), 7);

    letter = dna5::N;
    EXPECT_EQ(letter.to_rank(), 8);

    letter = gap::GAP;
    EXPECT_EQ(letter.to_rank(), 9);
}

TEST(union_alphabet_test, copy_assignment)
{
    using alphabet_t = union_alphabet<dna4, gap>;
    constexpr alphabet_t letter1{dna4::T};
    alphabet_t letter2, letter3;
    letter2 = letter1;
    letter3 = {letter1};

    EXPECT_EQ(letter1.to_rank(), 3);
    EXPECT_EQ(letter2.to_rank(), 3);
    EXPECT_EQ(letter3.to_rank(), 3);
}

TEST(union_alphabet_test, move_assignment)
{
    using alphabet_t = union_alphabet<dna4, gap>;
    alphabet_t letter1 = dna4::G;
    alphabet_t letter2{std::move(letter1)};

    EXPECT_EQ(letter2.to_rank(), 2);
}

TEST(union_alphabet_test, single_union)
{
    using alphabet_t = union_alphabet<dna4>;
    constexpr alphabet_t letter1{};

    EXPECT_EQ(letter1.to_rank(), 0);
}

TEST(union_alphabet_test, fulfills_concepts)
{
    static_assert(std::is_pod_v<union_alphabet<dna5, dna5>>);
    static_assert(std::is_trivial_v<union_alphabet<dna5, dna5>>);
    static_assert(std::is_trivially_copyable_v<union_alphabet<dna5, dna5>>);
    static_assert(std::is_standard_layout_v<union_alphabet<dna5, dna5>>);
    static_assert(alphabet_concept<union_alphabet<dna5, dna5>>);
}

TEST(union_alphabet_test, rank_type)
{
    using alphabet1_t = union_alphabet<dna4, dna5, gap>;
    using alphabet2_t = union_alphabet<gap, dna5, dna4>;
    using alphabet3_t = union_alphabet<gap>;

    constexpr auto expect1 = std::is_same_v<alphabet1_t::rank_type, uint8_t>;
    constexpr auto expect2 = std::is_same_v<alphabet2_t::rank_type, uint8_t>;
    constexpr auto expect3 = std::is_same_v<alphabet3_t::rank_type, bool>;

    EXPECT_TRUE(expect1);
    EXPECT_TRUE(expect2);
    EXPECT_TRUE(expect3);
}

TEST(union_alphabet_test, from_and_to_rank)
{
    using alphabet_t = union_alphabet<dna4, dna5, gap>;
    alphabet_t letter{};

    EXPECT_EQ(letter.assign_rank(0).to_rank(), 0);
    EXPECT_EQ(letter.assign_rank(1).to_rank(), 1);
    EXPECT_EQ(letter.assign_rank(2).to_rank(), 2);
    EXPECT_EQ(letter.assign_rank(3).to_rank(), 3);
    EXPECT_EQ(letter.assign_rank(4).to_rank(), 4);
    EXPECT_EQ(letter.assign_rank(5).to_rank(), 5);
    EXPECT_EQ(letter.assign_rank(6).to_rank(), 6);
    EXPECT_EQ(letter.assign_rank(7).to_rank(), 7);
    EXPECT_EQ(letter.assign_rank(8).to_rank(), 8);
    EXPECT_EQ(letter.assign_rank(9).to_rank(), 9);
}

TEST(union_alphabet_test, to_char)
{
    using alphabet_t = union_alphabet<dna4, dna5, gap>;
    alphabet_t letter{};

    EXPECT_EQ(letter.assign_rank(0).to_char(), 'A');
    EXPECT_EQ(letter.assign_rank(1).to_char(), 'C');
    EXPECT_EQ(letter.assign_rank(2).to_char(), 'G');
    EXPECT_EQ(letter.assign_rank(3).to_char(), 'T');
    EXPECT_EQ(letter.assign_rank(4).to_char(), 'A');
    EXPECT_EQ(letter.assign_rank(5).to_char(), 'C');
    EXPECT_EQ(letter.assign_rank(6).to_char(), 'G');
    EXPECT_EQ(letter.assign_rank(7).to_char(), 'T');
    EXPECT_EQ(letter.assign_rank(8).to_char(), 'N');
    EXPECT_EQ(letter.assign_rank(9).to_char(), '-');

//     EXPECT_EQ(letter.assign_rank(10).to_char(), static_cast<char>(0));
}

TEST(union_alphabet_test, relations)
{
    using alphabet_t = union_alphabet<dna4, dna5, gap>;
    alphabet_t letter1{};
    alphabet_t letter2{};

    EXPECT_EQ(letter1.assign_rank(0), letter2.assign_rank(0));
    EXPECT_EQ(letter1.assign_rank(5), letter2.assign_rank(5));
    EXPECT_NE(letter1.assign_rank(1), letter2.assign_rank(5));
    EXPECT_GT(letter1.assign_rank(2), letter2.assign_rank(1));
    EXPECT_GE(letter1.assign_rank(2), letter2.assign_rank(1));
    EXPECT_GE(letter1.assign_rank(2), letter2.assign_rank(2));
    EXPECT_LT(letter1.assign_rank(1), letter2.assign_rank(2));
    EXPECT_LE(letter1.assign_rank(1), letter2.assign_rank(2));
    EXPECT_LE(letter1.assign_rank(2), letter2.assign_rank(2));
}
