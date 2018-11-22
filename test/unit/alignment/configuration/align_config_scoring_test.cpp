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

#include <seqan3/alignment/configuration/align_config_scoring.hpp>
#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/std/concepts>

using namespace seqan3;

template <typename t>
class align_confg_scoring_test : public ::testing::Test
{
public:

    using alph_t = typename t::alphabet_type;
};

using test_types = ::testing::Types<aminoacid_scoring_scheme<int8_t>, nucleotide_scoring_scheme<int8_t>>;
TYPED_TEST_CASE(align_confg_scoring_test, test_types);

TYPED_TEST(align_confg_scoring_test, construction)
{
    EXPECT_TRUE((std::Constructible<align_cfg::scoring<TypeParam>, TypeParam>));
}

TYPED_TEST(align_confg_scoring_test, configuration)
{
    using alph_t = typename TestFixture::alph_t;
    {
        align_cfg::scoring elem{TypeParam{}};
        configuration cfg{elem};

        EXPECT_EQ((get<align_cfg::scoring>(cfg).value.score(assign_char(alph_t{}, 'a'),
                                                            assign_char(alph_t{}, 'a'))), 0);
    }

    {
        configuration cfg{align_cfg::scoring{TypeParam{}}};

        EXPECT_EQ((get<align_cfg::scoring>(cfg).value.score(assign_char(alph_t{}, 'a'),
                                                            assign_char(alph_t{}, 'c'))), -1);
    }
}
