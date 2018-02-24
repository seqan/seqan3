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

#include <sstream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

using namespace seqan3;

// These test case only test seqan3::gapped specific functions/properties that are not offered by the general
// `seqan3::alphabet_concept` interface. Those common interface function of `seqan3::gapped` will be tested
// in `alphabet/alphabet_test.hpp`.

TEST(gapped_test, initialise_from_component_alphabet)
{
    using alphabet_t = gapped<dna4>;

    constexpr alphabet_t letter0{dna4::A};
    constexpr alphabet_t letter1 = dna4::C;
    constexpr alphabet_t letter2 = {dna4::G};
    constexpr alphabet_t letter3 = static_cast<alphabet_t>(dna4::T);

    alphabet_t letter4{dna4::A};
    alphabet_t letter5 = dna4::C;
    alphabet_t letter6 = {dna4::G};
    alphabet_t letter7 = static_cast<alphabet_t>(dna4::T);

    constexpr alphabet_t letter8{gap::GAP}; // letter3 = dna4::T; does not work
    alphabet_t letter9{gap::GAP};

    EXPECT_EQ(letter0.to_rank(), 0);
    EXPECT_EQ(letter1.to_rank(), 1);
    EXPECT_EQ(letter2.to_rank(), 2);
    EXPECT_EQ(letter3.to_rank(), 3);
    EXPECT_EQ(letter4.to_rank(), 0);
    EXPECT_EQ(letter5.to_rank(), 1);
    EXPECT_EQ(letter6.to_rank(), 2);
    EXPECT_EQ(letter7.to_rank(), 3);
    EXPECT_EQ(letter8.to_rank(), 4);
    EXPECT_EQ(letter9.to_rank(), 4);
}

TEST(gapped_test, assign_from_component_alphabet)
{
    using alphabet_t = gapped<dna4>;
    alphabet_t letter{};

    letter = dna4::A;
    EXPECT_EQ(letter.to_rank(), 0);

    letter = {dna4::C}; // letter = {dna4::C}; does not work
    EXPECT_EQ(letter.to_rank(), 1);

    letter = static_cast<alphabet_t>(dna4::G);
    EXPECT_EQ(letter.to_rank(), 2);

    letter = {static_cast<alphabet_t>(dna4::T)};
    EXPECT_EQ(letter.to_rank(), 3);

    letter = gap::GAP;
    EXPECT_EQ(letter.to_rank(), 4);
}

TEST(gapped_test, fulfills_concepts)
{
    using alphabet_t = gapped<dna4>;
    EXPECT_TRUE((std::is_pod_v<alphabet_t>));
    EXPECT_TRUE((std::is_trivial_v<alphabet_t>));
    EXPECT_TRUE((std::is_trivially_copyable_v<alphabet_t>));
    EXPECT_TRUE((std::is_standard_layout_v<alphabet_t>));
}

TEST(gapped_test, stream_operator)
{
    using alphabet_t = gapped<dna4>;

    auto letterA = alphabet_t{dna4::A};
    auto letterC = alphabet_t{dna4::C};
    auto letterG = alphabet_t{dna4::G};
    auto letterT = alphabet_t{dna4::T};
    auto letterGAP = alphabet_t{gap::GAP};

    std::stringstream ss;
    ss << letterA << letterT << letterG << letterGAP << letterC;
    EXPECT_EQ(ss.str(), "ATG-C");
}

TEST(gapped_test, alphabet_name)
{
    using alphabet_t = gapped<dna4>;

    EXPECT_EQ(alphabet_name<alphabet_t>::value.string(), std::string{"gapped_dna4"});
    EXPECT_EQ(alphabet_name_v<alphabet_t>.string(), std::string{"gapped_dna4"});
}
