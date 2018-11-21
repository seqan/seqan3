// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
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

#include <seqan3/alphabet/composition/union_composition.hpp>
#include <seqan3/alphabet/gap/gap.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>

using namespace seqan3;

// These test case only test seqan3::union_composition specific functions/properties that are not offered by the general
// `seqan3::Alphabet` interface. Those common interface function of `seqan3::union_composition` will be tested
// in `alphabet/alphabet_test.hpp`.

TEST(union_composition_test, initialise_from_component_alphabet)
{
    dna5 l(rna5::A);

    using alphabet_t = union_composition<dna4, dna5, gap>;
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

TEST(union_composition_test, initialise_from_component_alphabet_subtype)
{
    using alphabet_t = union_composition<dna4, dna5, gap>;
    using variant_t = std::variant<dna4, dna5, gap>;

    variant_t variant0{rna4::A};
    alphabet_t letter0{rna4::A};

    variant_t variant1 = rna4::C;
    alphabet_t letter1 = rna4::C;

    variant_t variant2 = {rna4::G};
    alphabet_t letter2 = {rna4::G};

    variant_t variant3 = static_cast<variant_t>(rna4::T);
    alphabet_t letter3 = static_cast<alphabet_t>(rna4::T);

    variant_t variant4 = {static_cast<variant_t>(rna5::A)};
    alphabet_t letter4 = {static_cast<alphabet_t>(rna5::A)};

    variant_t variant5{rna5::C};
    alphabet_t letter5{rna5::C};

    variant_t variant6 = rna5::G;
    alphabet_t letter6 = rna5::G;

    variant_t variant7 = {rna5::T};
    alphabet_t letter7 = {rna5::T};

    variant_t variant8 = static_cast<variant_t>(rna5::N);
    alphabet_t letter8 = static_cast<alphabet_t>(rna5::N);

    EXPECT_EQ(letter0.to_rank(), 0);
    EXPECT_EQ(letter1.to_rank(), 1);
    EXPECT_EQ(letter2.to_rank(), 2);
    EXPECT_EQ(letter3.to_rank(), 3);
    EXPECT_EQ(letter4.to_rank(), 4);
    EXPECT_EQ(letter5.to_rank(), 5);
    EXPECT_EQ(letter6.to_rank(), 6);
    EXPECT_EQ(letter7.to_rank(), 7);
    EXPECT_EQ(letter8.to_rank(), 8);
}

TEST(union_composition_test, assign_from_component_alphabet)
{
    using alphabet_t = union_composition<dna4, dna5, gap>;
    using variant_t = std::variant<dna4, dna5, gap>;
    alphabet_t letter{};
    variant_t variant{};

    variant = dna4::A;
    letter = dna4::A;
    EXPECT_EQ(variant.index(), 0u);
    EXPECT_EQ(letter.to_rank(), 0);

    variant = {dna4::C};
    letter = {dna4::C};
    EXPECT_EQ(variant.index(), 0u);
    EXPECT_EQ(letter.to_rank(), 1);

    variant = static_cast<variant_t>(dna4::G);
    letter = static_cast<alphabet_t>(dna4::G);
    EXPECT_EQ(variant.index(), 0u);
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

TEST(union_composition_test, assign_from_component_alphabet_subtype)
{
    using alphabet_t = union_composition<dna4, dna5, gap>;
    using variant_t = std::variant<dna4, dna5, gap>;
    alphabet_t letter{};
    variant_t variant{};

    variant = rna4::A;
    letter = rna4::A;
    EXPECT_EQ(variant.index(), 0u);
    EXPECT_EQ(letter.to_rank(), 0);

    variant = {rna4::C};
    letter = {rna4::C};
    EXPECT_EQ(variant.index(), 0u);
    EXPECT_EQ(letter.to_rank(), 1);

    variant = static_cast<variant_t>(rna4::G);
    letter = static_cast<alphabet_t>(rna4::G);
    EXPECT_EQ(variant.index(), 0u);
    EXPECT_EQ(letter.to_rank(), 2);

    variant = {static_cast<variant_t>(rna4::T)};
    letter = {static_cast<alphabet_t>(rna4::T)};
    EXPECT_EQ(letter.to_rank(), 3);

    letter = rna5::A;
    EXPECT_EQ(letter.to_rank(), 4);

    letter = rna5::C;
    EXPECT_EQ(letter.to_rank(), 5);

    letter = rna5::G;
    EXPECT_EQ(letter.to_rank(), 6);

    letter = rna5::T;
    EXPECT_EQ(letter.to_rank(), 7);

    letter = rna5::N;
    EXPECT_EQ(letter.to_rank(), 8);
}

TEST(union_composition_test, compare_to_component_alphabet)
{
    using alphabet_t = union_composition<dna4, dna5>;

    constexpr alphabet_t letter0{dna4::G};

    EXPECT_EQ(letter0, dna4::G);
    EXPECT_NE(letter0, dna4::A);
    EXPECT_NE(letter0, dna5::A);

    EXPECT_EQ(dna4::G, letter0);
    EXPECT_NE(dna4::A, letter0);
    EXPECT_NE(dna5::A, letter0);
}

TEST(union_composition_test, compare_to_component_alphabet_subtype)
{
    using alphabet_t = union_composition<dna4, dna5>;

    constexpr alphabet_t letter0{dna4::G};

    EXPECT_EQ(letter0, rna4::G);
    EXPECT_NE(letter0, rna4::A);
    EXPECT_NE(letter0, rna5::A);

    EXPECT_EQ(rna4::G, letter0);
    EXPECT_NE(rna4::A, letter0);
    EXPECT_NE(rna5::A, letter0);
}

TEST(union_composition_test, fulfills_concepts)
{
    EXPECT_TRUE((Alphabet<union_composition<dna5, gap>>));
}

TEST(union_composition_test, rank_type)
{
    using alphabet1_t = union_composition<dna4, dna5, gap>;
    using alphabet2_t = union_composition<gap, dna5, dna4>;
    using alphabet3_t = union_composition<char, gap>;

    EXPECT_TRUE((std::is_same_v<alphabet1_t::rank_type, uint8_t>));
    EXPECT_TRUE((std::is_same_v<alphabet2_t::rank_type, uint8_t>));
    EXPECT_TRUE((std::is_same_v<alphabet3_t::rank_type, uint16_t>));
}

TEST(union_composition_test, value_size)
{
    using alphabet1_t = union_composition<dna4, dna5, gap>;
    using alphabet2_t = union_composition<gap, dna5, dna4>;
    using alphabet3_t = union_composition<char, gap>;

    EXPECT_TRUE((std::is_same_v<decltype(alphabet1_t::value_size), const uint8_t>));
    EXPECT_TRUE((std::is_same_v<decltype(alphabet2_t::value_size), const uint8_t>));
    EXPECT_TRUE((std::is_same_v<decltype(alphabet3_t::value_size), const uint16_t>));

    EXPECT_EQ(alphabet1_t::value_size, 10);
    EXPECT_EQ(alphabet2_t::value_size, 10);
    EXPECT_EQ(alphabet3_t::value_size, 257);
}

TEST(union_composition_test, convert_by_index)
{
    union_composition<dna4, dna5, gap> u;
    u = dna5::C;

    EXPECT_FALSE(u.is_alternative<0>());
    EXPECT_TRUE(u.is_alternative<1>());
    EXPECT_FALSE(u.is_alternative<2>());

    EXPECT_THROW(u.convert_to<0>(), std::bad_variant_access);
    EXPECT_NO_THROW(u.convert_to<1>());
    EXPECT_THROW(u.convert_to<2>(), std::bad_variant_access);

    dna5 out = u.convert_to<1>();
    EXPECT_EQ(out, dna5::C);

    u = gap::GAP;

    gap g = u.convert_unsafely_to<2>();
    EXPECT_EQ(g, gap::GAP);
}

TEST(union_composition_test, convert_by_type)
{
    union_composition<dna4, dna5, gap> u;
    u = dna5::C;

    EXPECT_FALSE(u.is_alternative<dna4>());
    EXPECT_TRUE(u.is_alternative<dna5>());
    EXPECT_FALSE(u.is_alternative<gap>());

    EXPECT_THROW(u.convert_to<dna4>(), std::bad_variant_access);
    EXPECT_NO_THROW(u.convert_to<dna5>());
    EXPECT_THROW(u.convert_to<gap>(), std::bad_variant_access);

    dna5 out = u.convert_to<dna5>();
    EXPECT_EQ(out, dna5::C);

    u = gap::GAP;
    gap g = u.convert_unsafely_to<gap>();
    EXPECT_EQ(g, gap::GAP);
}
