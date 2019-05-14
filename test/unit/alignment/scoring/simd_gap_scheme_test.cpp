// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <sstream>

#include <gtest/gtest.h>

#include <range/v3/algorithm/equal.hpp>

#include <seqan3/alignment/scoring/simd_gap_scheme.hpp>
#include <seqan3/alignment/scoring/gap_scheme_concept.hpp>
#include <seqan3/core/simd/simd.hpp>
#include <seqan3/test/simd_utility.hpp>

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
    EXPECT_TRUE((gap_scheme_concept<simd_gap_scheme<simd_type_t<int32_t>>>));
    // TODO enable when we allow floats for simd computation.
    // EXPECT_TRUE((gap_scheme_concept<simd_gap_scheme<typename simd::simd_type<float>::type> const &>));
}

TEST(gap_scheme, member_types)
{
    using simd_t = simd_type_t<int32_t>;
    simd_gap_scheme<simd_t> scheme{};

    using score_t = typename decltype(scheme)::score_type;
    EXPECT_TRUE((std::is_same_v<score_t, simd_t>));
}

TEST(gap_scheme, get_gap_score)
{
    using simd_t = simd_type_t<int32_t>;

    simd_gap_scheme<simd_t> scheme{};
    SIMD_EQ(scheme.get_gap_score(), simd::fill<simd_t>(-1));
    EXPECT_TRUE((std::is_same_v<typename decltype(scheme)::score_type &, decltype(scheme.get_gap_score())>));
}

TEST(gap_scheme, set_score_gap)
{
    using simd_t = simd_type_t<int32_t>;

    simd_gap_scheme<simd_t> scheme{};
    SIMD_EQ(scheme.get_gap_score(), simd::fill<simd_t>(-1));
    scheme.get_gap_score() = simd::fill<simd_t>(-2);
    SIMD_EQ(scheme.get_gap_score(), simd::fill<simd_t>(-2));
}

TEST(gap_scheme, get_gap_open_score)
{
    using simd_t = simd_type_t<int32_t>;

    simd_gap_scheme<simd_t> scheme{};
    SIMD_EQ(scheme.get_gap_open_score(), simd::fill<simd_t>(0));
    EXPECT_TRUE((std::is_same_v<typename decltype(scheme)::score_type &, decltype(scheme.get_gap_open_score())>));
}

TEST(gap_scheme, set_score_gap_open)
{
    using simd_t = simd_type_t<int32_t>;

    simd_gap_scheme<simd_t> scheme{};
    SIMD_EQ(scheme.get_gap_open_score(), simd::fill<simd_t>(0));
    scheme.get_gap_open_score() = simd::fill<simd_t>(-2);
    SIMD_EQ(scheme.get_gap_open_score(), simd::fill<simd_t>(-2));
}

TEST(gap_scheme, set_linear)
{
    using simd_t = simd_type_t<int32_t>;

    simd_gap_scheme<simd_t> scheme{gap_score{-2}};
    SIMD_EQ(scheme.get_gap_score(), simd::fill<simd_t>(-2));
    SIMD_EQ(scheme.get_gap_open_score(), simd::fill<simd_t>(0));

    scheme.set_linear(gap_score{-3});
    SIMD_EQ(scheme.get_gap_score(), simd::fill<simd_t>(-3));
    SIMD_EQ(scheme.get_gap_open_score(), simd::fill<simd_t>(0));
}

TEST(gap_scheme, set_affine)
{
    using simd_t = simd_type_t<int32_t>;

    simd_gap_scheme<simd_t> scheme{gap_score{-2}, gap_open_score{-4}};
    SIMD_EQ(scheme.get_gap_score(), simd::fill<simd_t>(-2));
    SIMD_EQ(scheme.get_gap_open_score(), simd::fill<simd_t>(-4));

    scheme.set_affine(gap_score{-3}, gap_open_score{-6});
    SIMD_EQ(scheme.get_gap_score(), simd::fill<simd_t>(-3));
    SIMD_EQ(scheme.get_gap_open_score(), simd::fill<simd_t>(-6));
}

TEST(gap_scheme, score)
{
    using simd_t = simd_type_t<int32_t>;

    simd_gap_scheme<simd_t> scheme{gap_score{-2}};
    EXPECT_EQ(scheme.score(0), 0);
    EXPECT_EQ(scheme.score(2), -4);
    EXPECT_EQ(scheme.score(5), -10);

    scheme.set_affine(gap_score{-3}, gap_open_score{-6});
    EXPECT_EQ(scheme.score(0), 0);
    EXPECT_EQ(scheme.score(2), -12);
    EXPECT_EQ(scheme.score(5), -21);
}

// TODO Currently not working.
// #if SEQAN3_WITH_CEREAL
// template <typename in_archive_t, typename out_archive_t, typename TypeParam>
// void do_serialisation(TypeParam const l)
// {
//     // This makes sure the file is also deleted if an exception is thrown in one of the tests below
//     // Generate unique file name.
//     test::tmp_filename filename{"simd_gap_scheme_cereal_test"};
//
//     {
//         std::ofstream os{filename.get_path(), std::ios::binary};
//         out_archive_t oarchive{os};
//         oarchive(l);
//     }
//
//     {
//         TypeParam in_l{};
//         std::ifstream is{filename.get_path(), std::ios::binary};
//         in_archive_t iarchive{is};
//         iarchive(in_l);
//         EXPECT_EQ(l, in_l);
//     }
// }
//
// TEST(gap_scheme, serialisation)
// {
//     using simd_t = simd_type_t<int32_t>;
//
//     simd_gap_scheme<simd_t> scheme1;
//     scheme1.set_linear(gap_score{-3});
//
//     do_serialisation<cereal::BinaryInputArchive,         cereal::BinaryOutputArchive>        (scheme1);
//     do_serialisation<cereal::PortableBinaryInputArchive, cereal::PortableBinaryOutputArchive>(scheme1);
//     do_serialisation<cereal::JSONInputArchive,           cereal::JSONOutputArchive>          (scheme1);
//     do_serialisation<cereal::XMLInputArchive,            cereal::XMLOutputArchive>           (scheme1);
//
//     scheme1.set_affine(gap_score{-3}, gap_open_score{-6});
//
//     do_serialisation<cereal::BinaryInputArchive,         cereal::BinaryOutputArchive>        (scheme1);
//     do_serialisation<cereal::PortableBinaryInputArchive, cereal::PortableBinaryOutputArchive>(scheme1);
//     do_serialisation<cereal::JSONInputArchive,           cereal::JSONOutputArchive>          (scheme1);
//     do_serialisation<cereal::XMLInputArchive,            cereal::XMLOutputArchive>           (scheme1);
// }
// #endif // SEQAN3_WITH_CEREAL
