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

#include <type_traits>
#include <utility>

#include <meta/meta.hpp>

#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/range/view/to_char.hpp>

#include <range/v3/view/generate_n.hpp>

using namespace seqan3;

TEST(align_pairwise, single_rng_lvalue)
{

    auto seq1 = "ACGTGATG"_dna4;
    auto seq2 = "AGTGATACT"_dna4;

    auto p = std::make_pair(seq1, seq2);

    {  // the score
        configuration cfg = align_cfg::edit | align_cfg::result{align_cfg::with_score};
        for (auto && res : align_pairwise(p, cfg))
        {
            EXPECT_EQ(res.get_score(), -4);
        }
    }

    {  // the trace
        configuration cfg = align_cfg::edit | align_cfg::result{align_cfg::with_trace};
        for (auto && res : align_pairwise(p, cfg))
        {
            EXPECT_EQ(res.get_score(), -4);
            auto [cmp1, cmp2] = res.get_end_coordinate();
            EXPECT_EQ((std::tie(cmp1, cmp2)), (std::tuple{7, 8}));
            auto && [gap1, gap2] = res.get_trace();
            EXPECT_EQ(std::string{gap1 | view::to_char}, "ACGTGATG--");
            EXPECT_EQ(std::string{gap2 | view::to_char}, "A-GTGATACT");
        }
    }
}

TEST(align_pairwise, single_view_lvalue)
{

    auto seq1 = "ACGTGATG"_dna4;
    auto seq2 = "AGTGATACT"_dna4;

    // auto p = std::make_pair(seq1, seq2);
    auto v = ranges::view::single(std::tie(seq1, seq2)) | ranges::view::bounded;

    {  // the score
        configuration cfg = align_cfg::edit | align_cfg::result{align_cfg::with_score};
        for (auto && res : align_pairwise(v, cfg))
        {
            EXPECT_EQ(res.get_score(), -4);
        }
    }
    {  // the trace
        configuration cfg = align_cfg::edit | align_cfg::result{align_cfg::with_trace};
        for (auto && res : align_pairwise(v, cfg))
        {
            EXPECT_EQ(res.get_score(), -4);
            auto && [gap1, gap2] = res.get_trace();
            EXPECT_EQ(std::string{gap1 | view::to_char}, "ACGTGATG--");
            EXPECT_EQ(std::string{gap2 | view::to_char}, "A-GTGATACT");
        }
    }
}

TEST(align_pairwise, multiple_rng_lvalue)
{

    auto seq1 = "ACGTGATG"_dna4;
    auto seq2 = "AGTGATACT"_dna4;

    auto p = std::make_pair(seq1, seq2);
    std::vector<decltype(p)> vec{10, p};

    configuration cfg = align_cfg::edit | align_cfg::result{align_cfg::with_trace};
    for (auto && res : align_pairwise(vec, cfg))
    {
        EXPECT_EQ(res.get_score(), -4);
        auto && [gap1, gap2] = res.get_trace();
        EXPECT_EQ(std::string{gap1 | view::to_char}, "ACGTGATG--");
        EXPECT_EQ(std::string{gap2 | view::to_char}, "A-GTGATACT");
    }
}
