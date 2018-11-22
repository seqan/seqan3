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

#include <sstream>

#include <gtest/gtest.h>

#include <range/v3/algorithm/equal.hpp>

#include <seqan3/alignment/scoring/gap_scheme.hpp>
#include <seqan3/alignment/scoring/gap_scheme_concept.hpp>

#if SEQAN3_WITH_CEREAL
#include <seqan3/test/tmp_filename.hpp>

#include <fstream>

#include <cereal/archives/binary.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/portable_binary.hpp>
#include <cereal/archives/xml.hpp>
#endif // SEQAN3_WITH_CEREAL

using namespace seqan3;

TEST(gap_scheme, concept_check)
{
    EXPECT_TRUE((gap_scheme_concept<gap_scheme<>>));
    EXPECT_TRUE((gap_scheme_concept<gap_scheme<int32_t> const>));
    EXPECT_TRUE((gap_scheme_concept<gap_scheme<float> const &>));
}

TEST(gap_scheme, constructors_and_type_deduction_guides)
{

    {
        gap_scheme scheme;
        EXPECT_TRUE((std::is_same_v<decltype(scheme), gap_scheme<int8_t>>));
    }

    {
        gap_scheme scheme{};
        EXPECT_TRUE((std::is_same_v<decltype(scheme), gap_scheme<int8_t>>));
    }

    {
        gap_scheme scheme{gap_score{-2}, gap_open_score{-4}};
        EXPECT_TRUE((std::is_same_v<decltype(scheme), gap_scheme<int8_t>>));
    }

    {
        gap_scheme scheme{gap_score{-2}};
        EXPECT_TRUE((std::is_same_v<decltype(scheme), gap_scheme<int8_t>>));
    }

    {
        gap_scheme scheme{gap_score{-2.}, gap_open_score{-4.}};
        EXPECT_TRUE((std::is_same_v<decltype(scheme), gap_scheme<float>>));
    }

    {
        gap_scheme scheme{gap_score{-2.}};
        EXPECT_TRUE((std::is_same_v<decltype(scheme), gap_scheme<float>>));
    }
}

TEST(gap_scheme, member_types)
{
    gap_scheme scheme{};

    using score_t = typename decltype(scheme)::score_type;
    EXPECT_TRUE((std::is_same_v<score_t, int8_t>));
}

TEST(gap_scheme, get_gap_score)
{
    gap_scheme scheme{};
    EXPECT_EQ(scheme.get_gap_score(), 0);
    EXPECT_TRUE((std::is_same_v<typename decltype(scheme)::score_type &, decltype(scheme.get_gap_score())>));
}

TEST(gap_scheme, set_score_gap)
{
    gap_scheme scheme{};
    EXPECT_EQ(scheme.get_gap_score(), 0);
    scheme.get_gap_score() = -2;
    EXPECT_EQ(scheme.get_gap_score(), -2);
}

TEST(gap_scheme, get_gap_open_score)
{
    gap_scheme scheme{};
    EXPECT_EQ(scheme.get_gap_open_score(), 0);
    EXPECT_TRUE((std::is_same_v<typename decltype(scheme)::score_type &, decltype(scheme.get_gap_open_score())>));
}

TEST(gap_scheme, set_score_gap_open)
{
    gap_scheme scheme{};
    EXPECT_EQ(scheme.get_gap_open_score(), 0);
    scheme.get_gap_open_score() = -2;
    EXPECT_EQ(scheme.get_gap_open_score(), -2);
}

TEST(gap_scheme, set_linear)
{
    gap_scheme scheme{gap_score{-2}};
    EXPECT_EQ(scheme.get_gap_score(), -2);
    EXPECT_EQ(scheme.get_gap_open_score(), 0);

    scheme.set_linear(gap_score{-3});
    EXPECT_EQ(scheme.get_gap_score(), -3);
    EXPECT_EQ(scheme.get_gap_open_score(), 0);
}

TEST(gap_scheme, set_affine)
{
    gap_scheme scheme{gap_score{-2}, gap_open_score{-4}};
    EXPECT_EQ(scheme.get_gap_score(), -2);
    EXPECT_EQ(scheme.get_gap_open_score(), -4);

    scheme.set_affine(gap_score{-3}, gap_open_score{-6});
    EXPECT_EQ(scheme.get_gap_score(), -3);
    EXPECT_EQ(scheme.get_gap_open_score(), -6);
}

TEST(gap_scheme, score)
{
    gap_scheme scheme{gap_score{-2}};
    EXPECT_EQ(scheme.score(0), 0);
    EXPECT_EQ(scheme.score(2), -4);
    EXPECT_EQ(scheme.score(5), -10);

    scheme.set_affine(gap_score{-3}, gap_open_score{-6});
    EXPECT_EQ(scheme.score(0), 0);
    EXPECT_EQ(scheme.score(2), -12);
    EXPECT_EQ(scheme.score(5), -21);
}

#if SEQAN3_WITH_CEREAL
template <typename in_archive_t, typename out_archive_t, typename TypeParam>
void do_serialisation(TypeParam const l)
{
    // This makes sure the file is also deleted if an exception is thrown in one of the tests below
    // Generate unique file name.
    test::tmp_filename filename{"scoring_scheme_cereal_test"};

    {
        std::ofstream os{filename.get_path(), std::ios::binary};
        out_archive_t oarchive{os};
        oarchive(l);
    }

    {
        TypeParam in_l{};
        std::ifstream is{filename.get_path(), std::ios::binary};
        in_archive_t iarchive{is};
        iarchive(in_l);
        EXPECT_EQ(l, in_l);
    }
}

TEST(gap_scheme, serialisation)
{
    gap_scheme scheme1;
    scheme1.set_linear(gap_score{-3});

    do_serialisation<cereal::BinaryInputArchive,         cereal::BinaryOutputArchive>        (scheme1);
    do_serialisation<cereal::PortableBinaryInputArchive, cereal::PortableBinaryOutputArchive>(scheme1);
    do_serialisation<cereal::JSONInputArchive,           cereal::JSONOutputArchive>          (scheme1);
    do_serialisation<cereal::XMLInputArchive,            cereal::XMLOutputArchive>           (scheme1);

    scheme1.set_affine(gap_score{-3}, gap_open_score{-6});

    do_serialisation<cereal::BinaryInputArchive,         cereal::BinaryOutputArchive>        (scheme1);
    do_serialisation<cereal::PortableBinaryInputArchive, cereal::PortableBinaryOutputArchive>(scheme1);
    do_serialisation<cereal::JSONInputArchive,           cereal::JSONOutputArchive>          (scheme1);
    do_serialisation<cereal::XMLInputArchive,            cereal::XMLOutputArchive>           (scheme1);
}
#endif // SEQAN3_WITH_CEREAL
