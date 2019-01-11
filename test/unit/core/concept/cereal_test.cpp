
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

TEST(cereal, cereal_input_archive_concept)
{
    EXPECT_FALSE((cereal_input_archive_concept<cereal::XMLOutputArchive>));
    EXPECT_FALSE((cereal_input_archive_concept<cereal::JSONOutputArchive>));
    EXPECT_FALSE((cereal_input_archive_concept<cereal::BinaryOutputArchive>));
    EXPECT_FALSE((cereal_input_archive_concept<cereal::PortableBinaryOutputArchive>));
    EXPECT_TRUE((cereal_input_archive_concept<cereal::XMLInputArchive>));
    EXPECT_TRUE((cereal_input_archive_concept<cereal::JSONInputArchive>));
    EXPECT_TRUE((cereal_input_archive_concept<cereal::BinaryInputArchive>));
    EXPECT_TRUE((cereal_input_archive_concept<cereal::PortableBinaryInputArchive>));
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

TEST(cereal, cereal_text_archive_concept)
{
    EXPECT_TRUE((cereal_text_archive_concept<cereal::XMLOutputArchive>));
    EXPECT_TRUE((cereal_text_archive_concept<cereal::JSONOutputArchive>));
    EXPECT_FALSE((cereal_text_archive_concept<cereal::BinaryOutputArchive>));
    EXPECT_FALSE((cereal_text_archive_concept<cereal::PortableBinaryOutputArchive>));
    EXPECT_TRUE((cereal_text_archive_concept<cereal::XMLInputArchive>));
    EXPECT_TRUE((cereal_text_archive_concept<cereal::JSONInputArchive>));
    EXPECT_FALSE((cereal_text_archive_concept<cereal::BinaryInputArchive>));
    EXPECT_FALSE((cereal_text_archive_concept<cereal::PortableBinaryInputArchive>));
}

struct my_struct{};

TEST(cereal, cerealisable_concept)
{
    EXPECT_TRUE((cerealisable_concept<int>));
    EXPECT_TRUE((cerealisable_concept<float>));

    // my_struct does not define any serialise functions
    EXPECT_FALSE((cerealisable_concept<my_struct>));

    // will be true, since <cereal/types/array.hpp> is included
    EXPECT_TRUE((cerealisable_concept<std::array<int, 10>>));
    // is false, because <cereal/types/vector.hpp> is not included
    EXPECT_FALSE((cerealisable_concept<std::vector<int>>));

    // recursive containers of cerealisable value types work
    EXPECT_TRUE((cerealisable_concept<std::array<std::array<int, 10>, 10>>));
}

#endif
