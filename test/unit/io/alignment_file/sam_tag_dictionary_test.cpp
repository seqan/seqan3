// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <type_traits>

#include <seqan3/io/alignment_file/sam_tag_dictionary.hpp>
#include <seqan3/std/concepts>

using seqan3::operator""_tag;

TEST(sam_tag, name_to_uint_conversion)
{
    // some example calculations
    EXPECT_EQ("NM"_tag, (uint16_t)('N' * 256 + 'M'));
    EXPECT_EQ("nm"_tag, (uint16_t)('n' * 256 + 'm'));
    EXPECT_EQ("N0"_tag, (uint16_t)('N' * 256 + '0'));
    EXPECT_EQ("N9"_tag, (uint16_t)('N' * 256 + '9'));
    EXPECT_EQ("AZ"_tag, (uint16_t)('A' * 256 + 'Z'));
    EXPECT_EQ("az"_tag, (uint16_t)('a' * 256 + 'z'));
    EXPECT_NE("NM"_tag, "nm"_tag); // case sensitivity
}

TEST(sam_tag_type, type_member_of_known_tags)
{
    // types according to the SAM specifications
    // (see https://samtools.github.io/hts-specs/SAMtags.pdf)
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"AM"_tag>::type, int32_t>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"AS"_tag>::type, int32_t>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"BC"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"BQ"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"BZ"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"CB"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"CC"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"CG"_tag>::type, std::vector<int32_t>>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"CM"_tag>::type, int32_t>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"CO"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"CP"_tag>::type, int32_t>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"CQ"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"CR"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"CS"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"CT"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"CY"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"E2"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"FI"_tag>::type, int32_t>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"FS"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"FZ"_tag>::type, std::vector<uint16_t>>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"H0"_tag>::type, int32_t>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"H1"_tag>::type, int32_t>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"H2"_tag>::type, int32_t>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"HI"_tag>::type, int32_t>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"IH"_tag>::type, int32_t>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"LB"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"MC"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"MD"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"MI"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"MQ"_tag>::type, int32_t>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"NH"_tag>::type, int32_t>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"NM"_tag>::type, int32_t>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"OC"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"OP"_tag>::type, int32_t>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"OQ"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"OX"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"PG"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"PQ"_tag>::type, int32_t>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"PT"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"PU"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"Q2"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"QT"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"QX"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"R2"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"RG"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"RT"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"RX"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"SA"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"SM"_tag>::type, int32_t>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"TC"_tag>::type, int32_t>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"U2"_tag>::type, std::string>));
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type<"UQ"_tag>::type, int32_t>));

    // test the short cut helper
    EXPECT_TRUE((std::same_as<seqan3::sam_tag_type_t<"AM"_tag>, seqan3::sam_tag_type<"AM"_tag>::type>));
}

TEST(sam_tag_dictionary, get_function_known_tag)
{
    seqan3::sam_tag_dictionary dict{};

    dict.get<"NM"_tag>() = 3;
    dict.get<"NM"_tag>() = 5; // overwrites previous
    dict.get<"CO"_tag>() = "comment";
    dict.get<"CG"_tag>() = std::vector<int32_t>{3, 4, 5};

    EXPECT_EQ(dict.get<"NM"_tag>(), 5);
    EXPECT_EQ(dict.get<"CO"_tag>(), "comment");
    EXPECT_EQ(dict.get<"CG"_tag>(), (std::vector<int32_t>{3, 4, 5}));
}

TEST(sam_tag_dictionary, get_function_unknown_tag)
{
    using variant_type = seqan3::sam_tag_dictionary::variant_type;

    seqan3::sam_tag_dictionary dict{};

    dict["nm"_tag] = 'a'; // lower case nm tag type is NOT defined so it can be either type
    dict["nm"_tag] = std::vector<int32_t>{3, 4, 5}; // overwrites previous
    dict["co"_tag] = std::string("comment");
    dict["cg"_tag] = std::vector<int32_t>{3, 4, 5};

    EXPECT_EQ(dict["nm"_tag], variant_type{(std::vector<int32_t>{3, 4, 5})});
    EXPECT_EQ(dict["co"_tag], variant_type{"comment"});
    EXPECT_EQ(dict["cg"_tag], variant_type{(std::vector<int32_t>{3, 4, 5})});
}

TEST(sam_tag_dictionary, get_function_const)
{
    seqan3::sam_tag_dictionary dict{};

    dict.get<"NM"_tag>() = 5;
    dict.get<"CO"_tag>() = "comment";
    dict.get<"CG"_tag>() = std::vector<int32_t>{3, 4, 5};

    auto const & dict2 = dict;

    EXPECT_EQ(dict2.get<"NM"_tag>(), 5);
    EXPECT_EQ(dict2.get<"CO"_tag>(), "comment");
    EXPECT_EQ(dict2.get<"CG"_tag>(), (std::vector<int32_t>{3, 4, 5}));
}

TEST(sam_tag_dictionary, get_function_temporary)
{
    seqan3::sam_tag_dictionary dict{};

    dict.get<"NM"_tag>() = 5;
    dict.get<"CO"_tag>() = "comment";
    dict.get<"CG"_tag>() = std::vector<int32_t>{3, 4, 5};

    EXPECT_TRUE((std::is_rvalue_reference_v<decltype(std::move(dict).get<"NM"_tag>())>));
    EXPECT_TRUE((std::is_rvalue_reference_v<decltype(std::move(dict).get<"CO"_tag>())>));
    EXPECT_TRUE((std::is_rvalue_reference_v<decltype(std::move(dict).get<"CG"_tag>())>));
}

TEST(sam_tag_dictionary, get_function_temporary_const)
{
    seqan3::sam_tag_dictionary dict{};

    dict.get<"NM"_tag>() = 5;
    dict.get<"CO"_tag>() = "comment";
    dict.get<"CG"_tag>() = std::vector<int32_t>{3, 4, 5};

    auto const & dict2 = dict;

    EXPECT_TRUE((std::is_rvalue_reference_v<decltype(std::move(dict2).get<"NM"_tag>())>));
    EXPECT_TRUE((std::is_rvalue_reference_v<decltype(std::move(dict2).get<"CO"_tag>())>));
    EXPECT_TRUE((std::is_rvalue_reference_v<decltype(std::move(dict2).get<"CG"_tag>())>));
}
