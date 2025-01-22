// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <algorithm>
#include <ranges>
#include <string>
#include <vector>

#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/nucleotide/rna15.hpp>
#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/alphabet/nucleotide/rna5.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/alphabet/views/translate.hpp>
#include <seqan3/test/expect_range_eq.hpp>

#include "../../range/iterator_test_template.hpp"
#include "../../range/range_test_template.hpp"

using seqan3::operator""_aa27;

template <seqan3::nucleotide_alphabet nucleotide_alphabet_t, seqan3::translation_frames translation_frames>
struct translate_single_view_test_fixture : public range_test_fixture
{
    using range_value_t = seqan3::aa27;
    using range_reference_t = range_value_t;

    using range_const_value_t = range_value_t;
    using range_const_reference_t = range_value_t;

    static constexpr bool input_range = true;
    static constexpr bool forward_range = true;
    static constexpr bool bidirectional_range = true;
    static constexpr bool random_access_range = true;
    static constexpr bool contiguous_range = false;
    static constexpr bool output_range = false;

    static constexpr bool common_range = true;
    static constexpr bool viewable_range = true;
    static constexpr bool view = true;
    static constexpr bool sized_range = true;
    static constexpr bool const_iterable_range = true;
    static constexpr bool size_member = true;
    static constexpr bool const_size_member = true;
    static constexpr bool subscript_member = true;

    seqan3::aa27_vector expected_range()
    {
        switch (translation_frames)
        {
        case seqan3::translation_frames::forward_frame0:
            return "TYVR"_aa27;
        case seqan3::translation_frames::reverse_frame0:
            return "YVRT"_aa27;
        case seqan3::translation_frames::forward_frame1:
            return "RTYV"_aa27;
        case seqan3::translation_frames::reverse_frame1:
            return "TYVR"_aa27;
        case seqan3::translation_frames::forward_frame2:
            return "VRT"_aa27;
        case seqan3::translation_frames::reverse_frame2:
            return "RTY"_aa27;
        }
        return {};
    }

    auto _underlying_range()
    {
        // string_view is a borrowed_range:
        return std::string_view{"ACGTACGTACGTA"} | seqan3::views::char_to<nucleotide_alphabet_t>;
    }

    auto range()
    {
        return _underlying_range() | seqan3::views::translate_single(translation_frames);
    }
};

using translate_single_view_test_fixtures_t =
    ::testing::Types<translate_single_view_test_fixture<seqan3::dna4, seqan3::translation_frames::forward_frame0>,
                     translate_single_view_test_fixture<seqan3::dna4, seqan3::translation_frames::reverse_frame0>,
                     translate_single_view_test_fixture<seqan3::dna4, seqan3::translation_frames::forward_frame1>,
                     translate_single_view_test_fixture<seqan3::dna4, seqan3::translation_frames::reverse_frame1>,
                     translate_single_view_test_fixture<seqan3::dna4, seqan3::translation_frames::forward_frame2>,
                     translate_single_view_test_fixture<seqan3::dna4, seqan3::translation_frames::reverse_frame2>>;

INSTANTIATE_TYPED_TEST_SUITE_P(translate_single_view_test, range_test, translate_single_view_test_fixtures_t, );
INSTANTIATE_TYPED_TEST_SUITE_P(translate_single_view_test, iterator_fixture, translate_single_view_test_fixtures_t, );

template <seqan3::nucleotide_alphabet nucleotide_alphabet_t, seqan3::translation_frames translation_frames>
struct translate_view_test_fixture : public range_test_fixture
{
    using _underlying_range_t = decltype(seqan3::views::char_to<nucleotide_alphabet_t>(std::string_view{}));

    using range_value_t = decltype(seqan3::views::translate_single(_underlying_range_t{}, translation_frames));
    using range_reference_t = range_value_t;

    using range_const_value_t = range_value_t;
    using range_const_reference_t = range_value_t;

    static constexpr bool input_range = true;
    static constexpr bool forward_range = true;
    static constexpr bool bidirectional_range = true;
    static constexpr bool random_access_range = true;
    static constexpr bool contiguous_range = false;
    static constexpr bool output_range = false;

    static constexpr bool common_range = true;
    static constexpr bool viewable_range = true;
    static constexpr bool view = true;
    static constexpr bool sized_range = true;
    static constexpr bool const_iterable_range = true;
    static constexpr bool size_member = true;
    static constexpr bool const_size_member = true;
    static constexpr bool subscript_member = true;

    template <typename range_value_t, typename expected_range_value_t>
    static void expect_range_value_equal(range_value_t && range_value, expected_range_value_t && expected_range_value)
    {
        EXPECT_RANGE_EQ(range_value, expected_range_value);
    }

    std::vector<seqan3::aa27_vector> expected_range()
    {
        switch (translation_frames)
        {
        case seqan3::translation_frames::forward_frame0:
            return {"TYVR"_aa27};
        case seqan3::translation_frames::reverse_frame0:
            return {"YVRT"_aa27};
        case seqan3::translation_frames::forward_frame1:
            return {"RTYV"_aa27};
        case seqan3::translation_frames::reverse_frame1:
            return {"TYVR"_aa27};
        case seqan3::translation_frames::forward_frame2:
            return {"VRT"_aa27};
        case seqan3::translation_frames::reverse_frame2:
            return {"RTY"_aa27};
        case seqan3::translation_frames::forward_reverse0:
            return {"TYVR"_aa27, "YVRT"_aa27};
        case seqan3::translation_frames::forward_reverse1:
            return {"RTYV"_aa27, "TYVR"_aa27};
        case seqan3::translation_frames::forward_reverse2:
            return {"VRT"_aa27, "RTY"_aa27};
        case seqan3::translation_frames::forward_frames:
            return {"TYVR"_aa27, "RTYV"_aa27, "VRT"_aa27};
        case seqan3::translation_frames::reverse_frames:
            return {"YVRT"_aa27, "TYVR"_aa27, "RTY"_aa27};
        case seqan3::translation_frames::six_frames:
            return {"TYVR"_aa27, "RTYV"_aa27, "VRT"_aa27, "YVRT"_aa27, "TYVR"_aa27, "RTY"_aa27};
        }

        return {};
    }

    auto _underlying_range()
    {
        // string_view is a borrowed_range:
        return std::string_view{"ACGTACGTACGTA"} | seqan3::views::char_to<nucleotide_alphabet_t>;
    }

    auto range()
    {
        return seqan3::views::translate(_underlying_range(), translation_frames);
    }
};

using translate_view_test_fixtures_t =
    ::testing::Types<translate_view_test_fixture<seqan3::dna4, seqan3::translation_frames::forward_frame0>,
                     translate_view_test_fixture<seqan3::dna4, seqan3::translation_frames::reverse_frame0>,
                     translate_view_test_fixture<seqan3::dna4, seqan3::translation_frames::forward_frame1>,
                     translate_view_test_fixture<seqan3::dna4, seqan3::translation_frames::reverse_frame1>,
                     translate_view_test_fixture<seqan3::dna4, seqan3::translation_frames::forward_frame2>,
                     translate_view_test_fixture<seqan3::dna4, seqan3::translation_frames::reverse_frame2>,
                     translate_view_test_fixture<seqan3::dna4, seqan3::translation_frames::forward_reverse0>,
                     translate_view_test_fixture<seqan3::dna4, seqan3::translation_frames::forward_reverse1>,
                     translate_view_test_fixture<seqan3::dna4, seqan3::translation_frames::forward_reverse2>,
                     translate_view_test_fixture<seqan3::dna4, seqan3::translation_frames::forward_frames>,
                     translate_view_test_fixture<seqan3::dna4, seqan3::translation_frames::reverse_frames>,
                     translate_view_test_fixture<seqan3::dna4, seqan3::translation_frames::six_frames>>;

INSTANTIATE_TYPED_TEST_SUITE_P(translate_view_test, range_test, translate_view_test_fixtures_t, );
INSTANTIATE_TYPED_TEST_SUITE_P(translate_view_test, iterator_fixture, translate_view_test_fixtures_t, );

template <typename T>
class translate_adaptor_test : public ::testing::Test
{};

// add all alphabets here
using nucleotide_types =
    ::testing::Types<seqan3::dna4, seqan3::dna5, seqan3::dna15, seqan3::rna4, seqan3::rna5, seqan3::rna15>;

TYPED_TEST_SUITE(translate_adaptor_test, nucleotide_types, );

TYPED_TEST(translate_adaptor_test, view_translate_single_exceptions)
{
    auto vec = std::string_view{"ACGTACGTACGTA"} | seqan3::views::char_to<TypeParam>;

    // Construct with multiple frames
    EXPECT_THROW(seqan3::views::translate_single(vec, seqan3::translation_frames::forward_frames),
                 std::invalid_argument);

    // Construct with default (empty) frame
    auto view = seqan3::views::translate_single(vec, seqan3::translation_frames{});
    EXPECT_THROW(view.size(), std::invalid_argument);
    EXPECT_THROW(view[0], std::invalid_argument);

    // Construct with default (empty) frame
    auto const const_view = seqan3::views::translate_single(vec, seqan3::translation_frames{});
    EXPECT_THROW(const_view.size(), std::invalid_argument);
    EXPECT_THROW(const_view[0], std::invalid_argument);
}

TYPED_TEST(translate_adaptor_test, view_translate_single)
{
    auto vec = std::string_view{"ACGTACGTACGTA"} | seqan3::views::char_to<TypeParam>;

    // default parameter translation_frames
    EXPECT_RANGE_EQ("TYVR"_aa27, vec | seqan3::views::translate_single);

    // default parameter translation_frames
    EXPECT_RANGE_EQ("TYVR"_aa27, vec | seqan3::views::translate_single());

    // single frame translation
    EXPECT_RANGE_EQ("TYVR"_aa27, vec | seqan3::views::translate_single(seqan3::translation_frames::forward_frame0));

    // function syntax
    EXPECT_RANGE_EQ("TYVR"_aa27, seqan3::views::translate_single(vec, seqan3::translation_frames::forward_frame0));

    // combinability
    EXPECT_RANGE_EQ("CMHA"_aa27,
                    vec | seqan3::views::complement
                        | seqan3::views::translate_single(seqan3::translation_frames::forward_frame0));

    // combinability
    EXPECT_RANGE_EQ("AHMC"_aa27,
                    vec | seqan3::views::complement
                        | seqan3::views::translate_single(seqan3::translation_frames::forward_frame0)
                        | std::views::reverse);
}

TYPED_TEST(translate_adaptor_test, view_translate)
{
    auto vec = std::string_view{"ACGTACGTACGTA"} | seqan3::views::char_to<TypeParam>;

    // default parameter translation_frames
    auto v1 = vec | seqan3::views::translate;
    EXPECT_EQ(v1.size(), 6u);
    EXPECT_RANGE_EQ(v1[0], "TYVR"_aa27);
    EXPECT_RANGE_EQ(v1[1], "RTYV"_aa27);
    EXPECT_RANGE_EQ(v1[2], "VRT"_aa27);
    EXPECT_RANGE_EQ(v1[3], "YVRT"_aa27);
    EXPECT_RANGE_EQ(v1[4], "TYVR"_aa27);
    EXPECT_RANGE_EQ(v1[5], "RTY"_aa27);

    // default parameter translation_frames
    auto v2 = vec | seqan3::views::translate();
    EXPECT_EQ(v2.size(), 6u);
    EXPECT_RANGE_EQ(v2[0], "TYVR"_aa27);
    EXPECT_RANGE_EQ(v2[1], "RTYV"_aa27);
    EXPECT_RANGE_EQ(v2[2], "VRT"_aa27);
    EXPECT_RANGE_EQ(v2[3], "YVRT"_aa27);
    EXPECT_RANGE_EQ(v2[4], "TYVR"_aa27);
    EXPECT_RANGE_EQ(v2[5], "RTY"_aa27);

    // single frame translation
    auto v3 = vec | seqan3::views::translate(seqan3::translation_frames::forward_frame0);
    EXPECT_EQ(v3.size(), 1u);
    EXPECT_RANGE_EQ(v3[0], "TYVR"_aa27);

    // six frame translation
    auto v4 = vec | seqan3::views::translate(seqan3::translation_frames::six_frames);
    EXPECT_EQ(v4.size(), 6u);
    EXPECT_RANGE_EQ(v4[0], "TYVR"_aa27);
    EXPECT_RANGE_EQ(v4[1], "RTYV"_aa27);
    EXPECT_RANGE_EQ(v4[2], "VRT"_aa27);
    EXPECT_RANGE_EQ(v4[3], "YVRT"_aa27);
    EXPECT_RANGE_EQ(v4[4], "TYVR"_aa27);
    EXPECT_RANGE_EQ(v4[5], "RTY"_aa27);

    // user-defined frame combination
    seqan3::translation_frames frames =
        seqan3::translation_frames::forward_frame0 | seqan3::translation_frames::forward_frame2;
    auto v5 = vec | seqan3::views::translate(frames);
    EXPECT_EQ(v5.size(), 2u);
    EXPECT_RANGE_EQ(v5[0], "TYVR"_aa27);
    EXPECT_RANGE_EQ(v5[1], "VRT"_aa27);

    // function syntax
    auto v6 = seqan3::views::translate(vec, seqan3::translation_frames::forward_reverse0);
    EXPECT_EQ(v6.size(), 2u);
    EXPECT_RANGE_EQ(v6[0], "TYVR"_aa27);
    EXPECT_RANGE_EQ(v6[1], "YVRT"_aa27);

    // combinability
    auto v7 = vec | seqan3::views::complement | seqan3::views::translate(seqan3::translation_frames::forward_reverse0);
    EXPECT_EQ(v7.size(), 2u);
    EXPECT_RANGE_EQ(v7[0], "CMHA"_aa27);
    EXPECT_RANGE_EQ(v7[1], "MHAC"_aa27);

    // combinability
    auto v8 = vec | seqan3::views::complement | seqan3::views::translate(seqan3::translation_frames::forward_reverse0)
            | std::views::take(1);
    EXPECT_EQ(v8.size(), 1u);
    EXPECT_RANGE_EQ(v8[0], "CMHA"_aa27);

    // combinability and function syntax
    auto v9 =
        seqan3::detail::view_translate(seqan3::views::complement(vec), seqan3::translation_frames::forward_reverse0);
    EXPECT_EQ(v9.size(), 2u);
    EXPECT_RANGE_EQ(v9[0], "CMHA"_aa27);
    EXPECT_RANGE_EQ(v9[1], "MHAC"_aa27);
}

// https://github.com/seqan/seqan3/issues/1339
TYPED_TEST(translate_adaptor_test, issue1339)
{
    // empty input
    auto vec = std::string_view{} | seqan3::views::char_to<TypeParam>;
    auto v1 = vec | seqan3::views::translate;

    EXPECT_EQ(v1.size(), 6u);
    for (auto && sequence : v1)
        EXPECT_TRUE(std::ranges::empty(sequence));

    // input of size 1
    auto vec2 = std::string_view{"A"} | seqan3::views::char_to<TypeParam>;
    auto v2 = vec2 | seqan3::views::translate;

    EXPECT_EQ(v2.size(), 6u);
    for (auto && sequence : v2)
        EXPECT_TRUE(std::ranges::empty(sequence));
}
