// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <meta/meta.hpp>

#include <seqan3/alphabet/aminoacid/all.hpp>
#include <seqan3/core/type_list.hpp>

using namespace seqan3;

template <typename T>
class aminoacid_conversion : public ::testing::Test
{};

using aminoacid_types = type_list<aa20, aa27>; // needed for some tests
using aminoacid_gtest_types = detail::transfer_template_args_onto_t<aminoacid_types, ::testing::Types>;

TYPED_TEST_CASE(aminoacid_conversion, aminoacid_gtest_types);

// conversion to any other amino acid type
TYPED_TEST(aminoacid_conversion, explicit_conversion)
{
    meta::for_each(aminoacid_types{}, [&] (auto && aa) constexpr
    {
        using out_type = std::decay_t<decltype(aa)>;
        EXPECT_EQ(static_cast<out_type>(TypeParam::A), out_type::A);
        EXPECT_EQ(static_cast<out_type>(TypeParam::C), out_type::C);
        EXPECT_EQ(static_cast<out_type>(TypeParam::D), out_type::D);
        EXPECT_EQ(static_cast<out_type>(TypeParam::E), out_type::E);
        EXPECT_EQ(static_cast<out_type>(TypeParam::F), out_type::F);
        EXPECT_EQ(static_cast<out_type>(TypeParam::G), out_type::G);
        EXPECT_EQ(static_cast<out_type>(TypeParam::H), out_type::H);
        EXPECT_EQ(static_cast<out_type>(TypeParam::I), out_type::I);
        EXPECT_EQ(static_cast<out_type>(TypeParam::K), out_type::K);
        EXPECT_EQ(static_cast<out_type>(TypeParam::L), out_type::L);
        EXPECT_EQ(static_cast<out_type>(TypeParam::M), out_type::M);
        EXPECT_EQ(static_cast<out_type>(TypeParam::N), out_type::N);
        EXPECT_EQ(static_cast<out_type>(TypeParam::P), out_type::P);
        EXPECT_EQ(static_cast<out_type>(TypeParam::Q), out_type::Q);
        EXPECT_EQ(static_cast<out_type>(TypeParam::R), out_type::R);
        EXPECT_EQ(static_cast<out_type>(TypeParam::S), out_type::S);
        EXPECT_EQ(static_cast<out_type>(TypeParam::T), out_type::T);
        EXPECT_EQ(static_cast<out_type>(TypeParam::V), out_type::V);
        EXPECT_EQ(static_cast<out_type>(TypeParam::W), out_type::W);
        EXPECT_EQ(static_cast<out_type>(TypeParam::Y), out_type::Y);
        if (std::is_same_v<TypeParam, aa27>)
        {
            EXPECT_EQ(static_cast<out_type>(TypeParam::B), out_type::B);
            EXPECT_EQ(static_cast<out_type>(TypeParam::J), out_type::J);
            EXPECT_EQ(static_cast<out_type>(TypeParam::O), out_type::O);
            EXPECT_EQ(static_cast<out_type>(TypeParam::U), out_type::U);
            EXPECT_EQ(static_cast<out_type>(TypeParam::X), out_type::X);
            EXPECT_EQ(static_cast<out_type>(TypeParam::Z), out_type::Z);
            EXPECT_EQ(static_cast<out_type>(TypeParam::TERMINATOR), out_type::TERMINATOR);
            EXPECT_EQ(static_cast<out_type>(TypeParam::UNKNOWN), out_type::UNKNOWN);
        }
        else if (std::is_same_v<TypeParam, aa20>)
        {
            EXPECT_EQ(static_cast<out_type>(TypeParam::B), out_type::D);
            EXPECT_EQ(static_cast<out_type>(TypeParam::J), out_type::L);
            EXPECT_EQ(static_cast<out_type>(TypeParam::O), out_type::L);
            EXPECT_EQ(static_cast<out_type>(TypeParam::U), out_type::C);
            EXPECT_EQ(static_cast<out_type>(TypeParam::X), out_type::S);
            EXPECT_EQ(static_cast<out_type>(TypeParam::Z), out_type::E);
            EXPECT_EQ(static_cast<out_type>(TypeParam::TERMINATOR), out_type::W);
            EXPECT_EQ(static_cast<out_type>(TypeParam::UNKNOWN), out_type::S);
        }
    });
}
