// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <array>
#include <vector>

#include <seqan3/core/concept/cereal.hpp>

#if SEQAN3_WITH_CEREAL
#include <cereal/archives/xml.hpp>
#include <cereal/archives/json.hpp>
#include <cereal/archives/binary.hpp>
#include <cereal/archives/portable_binary.hpp>

#include <cereal/types/array.hpp>
#endif

using namespace seqan3;

#if SEQAN3_WITH_CEREAL

TEST(cereal, cereal_output_archive)
{
    EXPECT_TRUE((cereal_output_archive<cereal::XMLOutputArchive>));
    EXPECT_TRUE((cereal_output_archive<cereal::JSONOutputArchive>));
    EXPECT_TRUE((cereal_output_archive<cereal::BinaryOutputArchive>));
    EXPECT_TRUE((cereal_output_archive<cereal::PortableBinaryOutputArchive>));
    EXPECT_FALSE((cereal_output_archive<cereal::XMLInputArchive>));
    EXPECT_FALSE((cereal_output_archive<cereal::JSONInputArchive>));
    EXPECT_FALSE((cereal_output_archive<cereal::BinaryInputArchive>));
    EXPECT_FALSE((cereal_output_archive<cereal::PortableBinaryInputArchive>));
}

TEST(cereal, cereal_input_archive)
{
    EXPECT_FALSE((cereal_input_archive<cereal::XMLOutputArchive>));
    EXPECT_FALSE((cereal_input_archive<cereal::JSONOutputArchive>));
    EXPECT_FALSE((cereal_input_archive<cereal::BinaryOutputArchive>));
    EXPECT_FALSE((cereal_input_archive<cereal::PortableBinaryOutputArchive>));
    EXPECT_TRUE((cereal_input_archive<cereal::XMLInputArchive>));
    EXPECT_TRUE((cereal_input_archive<cereal::JSONInputArchive>));
    EXPECT_TRUE((cereal_input_archive<cereal::BinaryInputArchive>));
    EXPECT_TRUE((cereal_input_archive<cereal::PortableBinaryInputArchive>));
}

TEST(cereal, cereal_archive)
{
    EXPECT_TRUE((cereal_archive<cereal::XMLOutputArchive>));
    EXPECT_TRUE((cereal_archive<cereal::JSONOutputArchive>));
    EXPECT_TRUE((cereal_archive<cereal::BinaryOutputArchive>));
    EXPECT_TRUE((cereal_archive<cereal::PortableBinaryOutputArchive>));
    EXPECT_TRUE((cereal_archive<cereal::XMLInputArchive>));
    EXPECT_TRUE((cereal_archive<cereal::JSONInputArchive>));
    EXPECT_TRUE((cereal_archive<cereal::BinaryInputArchive>));
    EXPECT_TRUE((cereal_archive<cereal::PortableBinaryInputArchive>));
}

TEST(cereal, cereal_text_archive)
{
    EXPECT_TRUE((cereal_text_archive<cereal::XMLOutputArchive>));
    EXPECT_TRUE((cereal_text_archive<cereal::JSONOutputArchive>));
    EXPECT_FALSE((cereal_text_archive<cereal::BinaryOutputArchive>));
    EXPECT_FALSE((cereal_text_archive<cereal::PortableBinaryOutputArchive>));
    EXPECT_TRUE((cereal_text_archive<cereal::XMLInputArchive>));
    EXPECT_TRUE((cereal_text_archive<cereal::JSONInputArchive>));
    EXPECT_FALSE((cereal_text_archive<cereal::BinaryInputArchive>));
    EXPECT_FALSE((cereal_text_archive<cereal::PortableBinaryInputArchive>));
}

struct my_struct{};

TEST(cereal, cerealisable)
{
    EXPECT_TRUE((cerealisable<int>));
    EXPECT_TRUE((cerealisable<float>));

    // my_struct does not define any serialise functions
    EXPECT_FALSE((cerealisable<my_struct>));

    // will be true, since <cereal/types/array.hpp> is included
    EXPECT_TRUE((cerealisable<std::array<int, 10>>));
    // is false, because <cereal/types/vector.hpp> is not included
    EXPECT_FALSE((cerealisable<std::vector<int>>));

    // recursive containers of cerealisable value types work
    EXPECT_TRUE((cerealisable<std::array<std::array<int, 10>, 10>>));
}

#endif
