// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
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

TEST(cereal, CerealOutputArchive)
{
    EXPECT_TRUE((CerealOutputArchive<cereal::XMLOutputArchive>));
    EXPECT_TRUE((CerealOutputArchive<cereal::JSONOutputArchive>));
    EXPECT_TRUE((CerealOutputArchive<cereal::BinaryOutputArchive>));
    EXPECT_TRUE((CerealOutputArchive<cereal::PortableBinaryOutputArchive>));
    EXPECT_FALSE((CerealOutputArchive<cereal::XMLInputArchive>));
    EXPECT_FALSE((CerealOutputArchive<cereal::JSONInputArchive>));
    EXPECT_FALSE((CerealOutputArchive<cereal::BinaryInputArchive>));
    EXPECT_FALSE((CerealOutputArchive<cereal::PortableBinaryInputArchive>));
}

TEST(cereal, CerealInputArchive)
{
    EXPECT_FALSE((CerealInputArchive<cereal::XMLOutputArchive>));
    EXPECT_FALSE((CerealInputArchive<cereal::JSONOutputArchive>));
    EXPECT_FALSE((CerealInputArchive<cereal::BinaryOutputArchive>));
    EXPECT_FALSE((CerealInputArchive<cereal::PortableBinaryOutputArchive>));
    EXPECT_TRUE((CerealInputArchive<cereal::XMLInputArchive>));
    EXPECT_TRUE((CerealInputArchive<cereal::JSONInputArchive>));
    EXPECT_TRUE((CerealInputArchive<cereal::BinaryInputArchive>));
    EXPECT_TRUE((CerealInputArchive<cereal::PortableBinaryInputArchive>));
}

TEST(cereal, CerealArchive)
{
    EXPECT_TRUE((CerealArchive<cereal::XMLOutputArchive>));
    EXPECT_TRUE((CerealArchive<cereal::JSONOutputArchive>));
    EXPECT_TRUE((CerealArchive<cereal::BinaryOutputArchive>));
    EXPECT_TRUE((CerealArchive<cereal::PortableBinaryOutputArchive>));
    EXPECT_TRUE((CerealArchive<cereal::XMLInputArchive>));
    EXPECT_TRUE((CerealArchive<cereal::JSONInputArchive>));
    EXPECT_TRUE((CerealArchive<cereal::BinaryInputArchive>));
    EXPECT_TRUE((CerealArchive<cereal::PortableBinaryInputArchive>));
}

TEST(cereal, CerealTextArchive)
{
    EXPECT_TRUE((CerealTextArchive<cereal::XMLOutputArchive>));
    EXPECT_TRUE((CerealTextArchive<cereal::JSONOutputArchive>));
    EXPECT_FALSE((CerealTextArchive<cereal::BinaryOutputArchive>));
    EXPECT_FALSE((CerealTextArchive<cereal::PortableBinaryOutputArchive>));
    EXPECT_TRUE((CerealTextArchive<cereal::XMLInputArchive>));
    EXPECT_TRUE((CerealTextArchive<cereal::JSONInputArchive>));
    EXPECT_FALSE((CerealTextArchive<cereal::BinaryInputArchive>));
    EXPECT_FALSE((CerealTextArchive<cereal::PortableBinaryInputArchive>));
}

struct my_struct{};

TEST(cereal, Cerealisable)
{
    EXPECT_TRUE((Cerealisable<int>));
    EXPECT_TRUE((Cerealisable<float>));

    // my_struct does not define any serialise functions
    EXPECT_FALSE((Cerealisable<my_struct>));

    // will be true, since <cereal/types/array.hpp> is included
    EXPECT_TRUE((Cerealisable<std::array<int, 10>>));
    // is false, because <cereal/types/vector.hpp> is not included
    EXPECT_FALSE((Cerealisable<std::vector<int>>));

    // recursive containers of cerealisable value types work
    EXPECT_TRUE((Cerealisable<std::array<std::array<int, 10>, 10>>));
}

#endif
