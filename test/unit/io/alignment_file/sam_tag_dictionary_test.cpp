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

#include <type_traits>

#include <seqan3/io/alignment_file/sam_tag_dictionary.hpp>
#include <seqan3/std/concepts>

using namespace seqan3;

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
    EXPECT_TRUE((std::Same<sam_tag_type<"AM"_tag>::type, int32_t>));
    EXPECT_TRUE((std::Same<sam_tag_type<"AS"_tag>::type, int32_t>));
    EXPECT_TRUE((std::Same<sam_tag_type<"BC"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"BQ"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"BZ"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"CB"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"CC"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"CG"_tag>::type, std::vector<int32_t>>));
    EXPECT_TRUE((std::Same<sam_tag_type<"CM"_tag>::type, int32_t>));
    EXPECT_TRUE((std::Same<sam_tag_type<"CO"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"CP"_tag>::type, int32_t>));
    EXPECT_TRUE((std::Same<sam_tag_type<"CQ"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"CR"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"CS"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"CT"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"CY"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"E2"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"FI"_tag>::type, int32_t>));
    EXPECT_TRUE((std::Same<sam_tag_type<"FS"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"FZ"_tag>::type, std::vector<uint16_t>>));
    EXPECT_TRUE((std::Same<sam_tag_type<"H0"_tag>::type, int32_t>));
    EXPECT_TRUE((std::Same<sam_tag_type<"H1"_tag>::type, int32_t>));
    EXPECT_TRUE((std::Same<sam_tag_type<"H2"_tag>::type, int32_t>));
    EXPECT_TRUE((std::Same<sam_tag_type<"HI"_tag>::type, int32_t>));
    EXPECT_TRUE((std::Same<sam_tag_type<"IH"_tag>::type, int32_t>));
    EXPECT_TRUE((std::Same<sam_tag_type<"LB"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"MC"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"MD"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"MI"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"MQ"_tag>::type, int32_t>));
    EXPECT_TRUE((std::Same<sam_tag_type<"NH"_tag>::type, int32_t>));
    EXPECT_TRUE((std::Same<sam_tag_type<"NM"_tag>::type, int32_t>));
    EXPECT_TRUE((std::Same<sam_tag_type<"OC"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"OP"_tag>::type, int32_t>));
    EXPECT_TRUE((std::Same<sam_tag_type<"OQ"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"OX"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"PG"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"PQ"_tag>::type, int32_t>));
    EXPECT_TRUE((std::Same<sam_tag_type<"PT"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"PU"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"Q2"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"QT"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"QX"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"R2"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"RG"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"RT"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"RX"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"SA"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"SM"_tag>::type, int32_t>));
    EXPECT_TRUE((std::Same<sam_tag_type<"TC"_tag>::type, int32_t>));
    EXPECT_TRUE((std::Same<sam_tag_type<"U2"_tag>::type, std::string>));
    EXPECT_TRUE((std::Same<sam_tag_type<"UQ"_tag>::type, int32_t>));

    // test the short cut helper
    EXPECT_TRUE((std::Same<sam_tag_type_t<"AM"_tag>, sam_tag_type<"AM"_tag>::type>));
}

TEST(sam_tag_dictionary, get_function_known_tag)
{
    sam_tag_dictionary dict{};

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
    using variant_type = std::variant<char, int32_t, float, std::string,
                                      std::vector<uint8_t>, std::vector<int8_t>,
                                      std::vector<uint16_t>, std::vector<int16_t>,
                                      std::vector<uint32_t>, std::vector<int32_t>,
                                      std::vector<float>>;
    sam_tag_dictionary dict{};

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
    sam_tag_dictionary dict{};

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
    sam_tag_dictionary dict{};

    dict.get<"NM"_tag>() = 5;
    dict.get<"CO"_tag>() = "comment";
    dict.get<"CG"_tag>() = std::vector<int32_t>{3, 4, 5};

    EXPECT_TRUE((std::is_rvalue_reference_v<decltype(std::move(dict).get<"NM"_tag>())>));
    EXPECT_TRUE((std::is_rvalue_reference_v<decltype(std::move(dict).get<"CO"_tag>())>));
    EXPECT_TRUE((std::is_rvalue_reference_v<decltype(std::move(dict).get<"CG"_tag>())>));
}

TEST(sam_tag_dictionary, get_function_temporary_const)
{
    sam_tag_dictionary dict{};

    dict.get<"NM"_tag>() = 5;
    dict.get<"CO"_tag>() = "comment";
    dict.get<"CG"_tag>() = std::vector<int32_t>{3, 4, 5};

    auto const & dict2 = dict;

    EXPECT_TRUE((std::is_rvalue_reference_v<decltype(std::move(dict2).get<"NM"_tag>())>));
    EXPECT_TRUE((std::is_rvalue_reference_v<decltype(std::move(dict2).get<"CO"_tag>())>));
    EXPECT_TRUE((std::is_rvalue_reference_v<decltype(std::move(dict2).get<"CG"_tag>())>));
}
