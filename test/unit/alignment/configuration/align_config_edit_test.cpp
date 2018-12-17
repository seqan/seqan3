// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
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
// ============================================================================

#include <gtest/gtest.h>

#include <seqan3/alignment/configuration/align_config_edit.hpp>

using namespace seqan3;

TEST(align_cfg_edit, is_global)
{
    EXPECT_EQ((std::is_same_v<std::remove_reference_t<decltype(get<align_cfg::mode>(align_cfg::edit).value)>,
                              detail::global_alignment_type>),
               true);
}

TEST(align_cfg_edit, is_hamming)
{
    auto scheme = get<align_cfg::scoring>(align_cfg::edit).value;

    for (unsigned i = 0; i < decltype(scheme)::matrix_size; ++i)
    {
        for (unsigned j = 0; j < decltype(scheme)::matrix_size; ++j)
        {
            if (i == j)
                EXPECT_EQ((scheme.score(assign_rank(dna15{}, i), assign_rank(dna15{}, j))), 0);
            else
                EXPECT_EQ((scheme.score(assign_rank(dna15{}, i), assign_rank(dna15{}, j))), -1);
        }
    }
}

TEST(align_cfg_edit, is_simple_gap)
{
    auto scheme = get<align_cfg::gap>(align_cfg::edit).value;
    EXPECT_EQ(scheme.get_gap_score(), -1);
    EXPECT_EQ(scheme.get_gap_open_score(), 0);
}
