// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/alphabet/detail/debug_stream_alphabet.hpp>
#include <seqan3/alphabet/exception.hpp>
#include <seqan3/test/pretty_printing.hpp>

template <typename T>
using alphabet = ::testing::Test;

constexpr size_t max_iterations = 65536u;

TYPED_TEST_SUITE_P(alphabet);

TYPED_TEST_P(alphabet, concept_check)
{
    EXPECT_TRUE(seqan3::alphabet<TypeParam>);
    EXPECT_TRUE(seqan3::alphabet<TypeParam &>);
    EXPECT_TRUE(seqan3::alphabet<TypeParam const>);
    EXPECT_TRUE(seqan3::alphabet<TypeParam const &>);

    EXPECT_TRUE(seqan3::writable_alphabet<TypeParam>);
    EXPECT_TRUE(seqan3::writable_alphabet<TypeParam &>);
    EXPECT_FALSE(seqan3::writable_alphabet<TypeParam const>);
    EXPECT_FALSE(seqan3::writable_alphabet<TypeParam const &>);
}

TYPED_TEST_P(alphabet, assign_char_to)
{
    using char_t = seqan3::alphabet_char_t<TypeParam>;
    if constexpr (std::integral<char_t>)
    {
        char_t i = std::numeric_limits<char_t>::min();
        char_t j = std::numeric_limits<char_t>::max();

        TypeParam t0;
        for (size_t k = 0; i < j && k < max_iterations; ++i, ++k)
            seqan3::assign_char_to(i, t0);

        EXPECT_TRUE((std::is_same_v<decltype(seqan3::assign_char_to(0, t0)), TypeParam &>));
        EXPECT_TRUE((std::is_same_v<decltype(seqan3::assign_char_to(0, TypeParam{})), TypeParam>));
    }
}

TYPED_TEST_P(alphabet, char_is_valid_for) // only test negative example for most; more inside specialised tests
{
    if constexpr (seqan3::alphabet_size<TypeParam> < 255) // includes most of our alphabets, but not the adaptations!
    {
        EXPECT_FALSE((seqan3::char_is_valid_for<TypeParam>(0))); // for none of our alphabets char{0} is valid
    }
}

TYPED_TEST_P(alphabet, assign_char_strictly_to)
{
    using char_t = seqan3::alphabet_char_t<TypeParam>;
    if constexpr (std::integral<char_t>)
    {
        char_t i = std::numeric_limits<char_t>::min();
        char_t j = std::numeric_limits<char_t>::max();

        for (size_t k = 0; i < j && k < max_iterations; ++i, ++k)
        {
            if (seqan3::char_is_valid_for<TypeParam>(i))
                EXPECT_NO_THROW(seqan3::assign_char_strictly_to(i, TypeParam{}));
            else
                EXPECT_THROW(seqan3::assign_char_strictly_to(i, TypeParam{}), seqan3::invalid_char_assignment);
        }
    }
}

TYPED_TEST_P(alphabet, to_char)
{
    TypeParam t0;
    EXPECT_TRUE((std::is_same_v<decltype(seqan3::to_char(t0)), seqan3::alphabet_char_t<TypeParam>>));

    // more elaborate tests are done in specific alphabets
}

REGISTER_TYPED_TEST_SUITE_P(alphabet,
                            concept_check,
                            assign_char_to,
                            char_is_valid_for,
                            assign_char_strictly_to,
                            to_char);
