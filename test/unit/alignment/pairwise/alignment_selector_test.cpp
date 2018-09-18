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

#include <list>
#include <type_traits>
#include <utility>
#include <vector>

#include <meta/meta.hpp>

#include <seqan3/alignment/configuration/align_config_edit.hpp>
#include <seqan3/alignment/configuration/align_config_output.hpp>
#include <seqan3/alignment/pairwise/alignment_selector.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/range/view/persist.hpp>
#include <seqan3/range/view/to_char.hpp>

using namespace seqan3;

TEST(alignment_selector, determine_result_type)
{
    using seq1_t = std::vector<dna4>;
    using seq2_t = std::list<dna4>;

    { // test case I
        auto cfg = align_cfg::edit;
        using _t = typename detail::determine_result_type<seq1_t, seq2_t, decltype(cfg)>::type;

        EXPECT_EQ(std::tuple_size_v<_t>, 2);
        EXPECT_TRUE((std::is_same_v<std::tuple_element_t<0, _t>, uint32_t>));
        EXPECT_TRUE((std::is_same_v<std::tuple_element_t<1, _t>, int32_t>));
    }

    { // test case II
        auto cfg = align_cfg::edit | align_cfg::output<align_result_key::score>;
        using _t = typename detail::determine_result_type<seq1_t, seq2_t, decltype(cfg)>::type;

        EXPECT_EQ(std::tuple_size_v<_t>, 2);
        EXPECT_TRUE((std::is_same_v<std::tuple_element_t<0, _t>, uint32_t>));
        EXPECT_TRUE((std::is_same_v<std::tuple_element_t<1, _t>, int32_t>));
    }

    { // test case III
        auto cfg = align_cfg::edit | align_cfg::output<align_result_key::trace>;
        using _t = typename detail::determine_result_type<seq1_t, seq2_t, decltype(cfg)>::type;

        using gapped_seq1_t = std::vector<gapped<dna4>>;
        using gapped_seq2_t = std::vector<gapped<dna4>>;

        EXPECT_EQ(std::tuple_size_v<_t>, 5);
        EXPECT_TRUE((std::is_same_v<std::tuple_element_t<0, _t>, uint32_t>));
        EXPECT_TRUE((std::is_same_v<std::tuple_element_t<1, _t>, int32_t>));
        EXPECT_TRUE((std::is_same_v<std::tuple_element_t<2, _t>, std::pair<size_t, size_t>>));
        EXPECT_TRUE((std::is_same_v<std::tuple_element_t<3, _t>, std::pair<size_t, size_t>>));
        EXPECT_TRUE((std::is_same_v<std::tuple_element_t<4, _t>, std::tuple<gapped_seq1_t, gapped_seq2_t>>));
    }
}

TEST(alignment_selector, select)
{
    using seq1_t = std::vector<dna4>;
    using seq2_t = std::vector<dna4>;
    using seq_t = std::pair<seq1_t, seq2_t>;

    { // test case I
        auto cfg = align_cfg::edit;
        using _t = detail::alignment_selector<seq_t, decltype(cfg)>;

        using res_t = typename detail::determine_result_type<seq1_t, seq2_t, decltype(cfg)>::type;

        seq_t s{};
        EXPECT_TRUE((std::is_same_v<res_t, typename _t::result_type>));
        EXPECT_TRUE((std::is_same_v<decltype(_t{cfg}.select(s)),
                                    std::function<res_t(res_t &)>>));
    }
}
