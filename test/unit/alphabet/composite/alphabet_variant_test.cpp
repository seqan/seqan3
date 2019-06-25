// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <variant>

#include <seqan3/alphabet/composite/alphabet_variant.hpp>
#include <seqan3/alphabet/gap/gap.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>

#include "../alphabet_test_template.hpp"
#include "../alphabet_constexpr_test_template.hpp"

using namespace seqan3;

using alphabet_variant_types = ::testing::Types<alphabet_variant<dna4, gap>,
                                                 alphabet_variant<dna4, dna5, gap>,
                                                 alphabet_variant<char, gap>>;

INSTANTIATE_TYPED_TEST_CASE_P(alphabet_variant, alphabet, alphabet_variant_types);
INSTANTIATE_TYPED_TEST_CASE_P(alphabet_variant, alphabet_constexpr, alphabet_variant_types);

template <typename T>
using alphabet_variant_test = ::testing::Test;

TYPED_TEST_CASE(alphabet_variant_test, alphabet_variant_types);

TEST(alphabet_variant_test, initialise_from_component_alphabet)
{
    dna5 l('A'_rna5);

    using alphabet_t = alphabet_variant<dna4, dna5, gap>;
    using variant_t = std::variant<dna4, dna5, gap>;

    constexpr variant_t variant0{'A'_dna4};
    constexpr alphabet_t letter0{'A'_dna4};

    constexpr variant_t variant1 = 'C'_dna4;
    constexpr alphabet_t letter1 = 'C'_dna4;

    constexpr variant_t variant2 = {'G'_dna4};
    constexpr alphabet_t letter2 = {'G'_dna4};

    constexpr variant_t variant3 = static_cast<variant_t>('T'_dna4);
    constexpr alphabet_t letter3 = static_cast<alphabet_t>('T'_dna4);

    constexpr variant_t variant4 = {static_cast<variant_t>('A'_dna5)};
    constexpr alphabet_t letter4 = {static_cast<alphabet_t>('A'_dna5)};

    variant_t variant5{'C'_dna5};
    alphabet_t letter5{'C'_dna5};

    variant_t variant6 = 'G'_dna5;
    alphabet_t letter6 = 'G'_dna5;

    variant_t variant7 = {'N'_dna5};
    alphabet_t letter7 = {'N'_dna5};

    variant_t variant8 = static_cast<variant_t>('T'_dna5);
    alphabet_t letter8 = static_cast<alphabet_t>('T'_dna5);

    variant_t variant9 = {static_cast<variant_t>(gap{})};
    alphabet_t letter9 = {static_cast<alphabet_t>(gap{})};

    EXPECT_EQ(letter0.to_rank(), 0);
    EXPECT_EQ(letter1.to_rank(), 1);
    EXPECT_EQ(letter2.to_rank(), 2);
    EXPECT_EQ(letter3.to_rank(), 3);
    EXPECT_EQ(letter4.to_rank(), 4);
    EXPECT_EQ(letter5.to_rank(), 5);
    EXPECT_EQ(letter6.to_rank(), 6);
    EXPECT_EQ(letter7.to_rank(), 7);
    EXPECT_EQ(letter8.to_rank(), 8);
    EXPECT_EQ(letter9.to_rank(), 9);
}

TEST(alphabet_variant_test, initialise_from_component_alphabet_subtype)
{
    using alphabet_t = alphabet_variant<dna4, dna5, gap>;
    using variant_t = std::variant<dna4, dna5, gap>;

    variant_t variant0{'A'_rna4};
    alphabet_t letter0{'A'_rna4};

    variant_t variant1 = 'C'_rna4;
    alphabet_t letter1 = 'C'_rna4;

    variant_t variant2 = {'G'_rna4};
    alphabet_t letter2 = {'G'_rna4};

    variant_t variant3 = static_cast<variant_t>('T'_rna4);
    alphabet_t letter3 = static_cast<alphabet_t>('T'_rna4);

    variant_t variant4 = {static_cast<variant_t>('A'_rna5)};
    alphabet_t letter4 = {static_cast<alphabet_t>('A'_rna5)};

    variant_t variant5{'C'_rna5};
    alphabet_t letter5{'C'_rna5};

    variant_t variant6 = 'G'_rna5;
    alphabet_t letter6 = 'G'_rna5;

    variant_t variant7 = {'N'_rna5};
    alphabet_t letter7 = {'N'_rna5};

    variant_t variant8 = static_cast<variant_t>('T'_rna5);
    alphabet_t letter8 = static_cast<alphabet_t>('T'_rna5);

    EXPECT_EQ(letter0.to_rank(), 0);
    EXPECT_EQ(letter1.to_rank(), 1);
    EXPECT_EQ(letter2.to_rank(), 2);
    EXPECT_EQ(letter3.to_rank(), 3);
    EXPECT_EQ(letter4.to_rank(), 4);
    EXPECT_EQ(letter5.to_rank(), 5);
    EXPECT_EQ(letter6.to_rank(), 6);
    EXPECT_EQ(letter7.to_rank(), 7);
    EXPECT_EQ(letter8.to_rank(), 8);
}

TEST(alphabet_variant_test, assign_from_component_alphabet)
{
    using alphabet_t = alphabet_variant<dna4, dna5, gap>;
    using variant_t = std::variant<dna4, dna5, gap>;
    alphabet_t letter{};
    variant_t variant{};

    variant = 'A'_dna4;
    letter = 'A'_dna4;
    EXPECT_EQ(variant.index(), 0u);
    EXPECT_EQ(letter.to_rank(), 0);

    variant = {'C'_dna4};
    letter = {'C'_dna4};
    EXPECT_EQ(variant.index(), 0u);
    EXPECT_EQ(letter.to_rank(), 1);

    variant = static_cast<variant_t>('G'_dna4);
    letter = static_cast<alphabet_t>('G'_dna4);
    EXPECT_EQ(variant.index(), 0u);
    EXPECT_EQ(letter.to_rank(), 2);

    variant = {static_cast<variant_t>('T'_dna4)};
    letter = {static_cast<alphabet_t>('T'_dna4)};
    EXPECT_EQ(letter.to_rank(), 3);

    letter = 'A'_dna5;
    EXPECT_EQ(letter.to_rank(), 4);

    letter = 'C'_dna5;
    EXPECT_EQ(letter.to_rank(), 5);

    letter = 'G'_dna5;
    EXPECT_EQ(letter.to_rank(), 6);

    letter = 'N'_dna5;
    EXPECT_EQ(letter.to_rank(), 7);

    letter = 'T'_dna5;
    EXPECT_EQ(letter.to_rank(), 8);

    letter = gap{};
    EXPECT_EQ(letter.to_rank(), 9);
}

TEST(alphabet_variant_test, assign_from_component_alphabet_subtype)
{
    using alphabet_t = alphabet_variant<dna4, dna5, gap>;
    using variant_t = std::variant<dna4, dna5, gap>;
    alphabet_t letter{};
    variant_t variant{};

    variant = 'A'_rna4;
    letter = 'A'_rna4;
    EXPECT_EQ(variant.index(), 0u);
    EXPECT_EQ(letter.to_rank(), 0);

    variant = {'C'_rna4};
    letter = {'C'_rna4};
    EXPECT_EQ(variant.index(), 0u);
    EXPECT_EQ(letter.to_rank(), 1);

    variant = static_cast<variant_t>('G'_rna4);
    letter = static_cast<alphabet_t>('G'_rna4);
    EXPECT_EQ(variant.index(), 0u);
    EXPECT_EQ(letter.to_rank(), 2);

    variant = {static_cast<variant_t>('T'_rna4)};
    letter = {static_cast<alphabet_t>('T'_rna4)};
    EXPECT_EQ(letter.to_rank(), 3);

    letter = 'A'_rna5;
    EXPECT_EQ(letter.to_rank(), 4);

    letter = 'C'_rna5;
    EXPECT_EQ(letter.to_rank(), 5);

    letter = 'G'_rna5;
    EXPECT_EQ(letter.to_rank(), 6);

    letter = 'N'_rna5;
    EXPECT_EQ(letter.to_rank(), 7);

    letter = 'T'_rna5;
    EXPECT_EQ(letter.to_rank(), 8);
}

TEST(alphabet_variant_test, compare_to_component_alphabet)
{
    using alphabet_t = alphabet_variant<dna4, dna5>;

    constexpr alphabet_t letter0{'G'_dna4};

    EXPECT_EQ(letter0, 'G'_dna4);
    EXPECT_NE(letter0, 'A'_dna4);
    EXPECT_NE(letter0, 'A'_dna5);

    EXPECT_EQ('G'_dna4, letter0);
    EXPECT_NE('A'_dna4, letter0);
    EXPECT_NE('A'_dna5, letter0);
}

TEST(alphabet_variant_test, compare_to_component_alphabet_subtype)
{
    using alphabet_t = alphabet_variant<dna4, dna5>;

    constexpr alphabet_t letter0{'G'_dna4};

    EXPECT_EQ(letter0, 'G'_rna4);
    EXPECT_NE(letter0, 'A'_rna4);
    EXPECT_NE(letter0, 'A'_rna5);

    EXPECT_EQ('G'_rna4, letter0);
    EXPECT_NE('A'_rna4, letter0);
    EXPECT_NE('A'_rna5, letter0);
}

TEST(alphabet_variant_test, rank_type)
{
    using alphabet1_t = alphabet_variant<dna4, dna5, gap>;
    using alphabet2_t = alphabet_variant<gap, dna5, dna4>;
    using alphabet3_t = alphabet_variant<char, gap>;

    EXPECT_TRUE((std::is_same_v<alphabet1_t::rank_type, uint8_t>));
    EXPECT_TRUE((std::is_same_v<alphabet2_t::rank_type, uint8_t>));
    EXPECT_TRUE((std::is_same_v<alphabet3_t::rank_type, uint16_t>));
}

TEST(alphabet_variant_test, alphabet_size_)
{
    using alphabet1_t = alphabet_variant<dna4, dna5, gap>;
    using alphabet2_t = alphabet_variant<gap, dna5, dna4>;
    using alphabet3_t = alphabet_variant<char, gap>;

    EXPECT_TRUE((std::is_same_v<decltype(alphabet1_t::alphabet_size), const uint8_t>));
    EXPECT_TRUE((std::is_same_v<decltype(alphabet2_t::alphabet_size), const uint8_t>));
    EXPECT_TRUE((std::is_same_v<decltype(alphabet3_t::alphabet_size), const uint16_t>));

    EXPECT_EQ(alphabet1_t::alphabet_size, 10);
    EXPECT_EQ(alphabet2_t::alphabet_size, 10);
    EXPECT_EQ(alphabet3_t::alphabet_size, 257);
}

TEST(alphabet_variant_test, convert_by_index)
{
    alphabet_variant<dna4, dna5, gap> u;
    u = 'C'_dna5;

    EXPECT_FALSE(u.is_alternative<0>());
    EXPECT_TRUE(u.is_alternative<1>());
    EXPECT_FALSE(u.is_alternative<2>());

    EXPECT_THROW(u.convert_to<0>(), std::bad_variant_access);
    EXPECT_NO_THROW(u.convert_to<1>());
    EXPECT_THROW(u.convert_to<2>(), std::bad_variant_access);

    dna5 out = u.convert_to<1>();
    EXPECT_EQ(out, 'C'_dna5);

    u = gap{};

    gap g = u.convert_unsafely_to<2>();
    EXPECT_EQ(g, gap{});
}

TEST(alphabet_variant_test, convert_by_type)
{
    alphabet_variant<dna4, dna5, gap> u;
    u = 'C'_dna5;

    EXPECT_FALSE(u.is_alternative<dna4>());
    EXPECT_TRUE(u.is_alternative<dna5>());
    EXPECT_FALSE(u.is_alternative<gap>());

    EXPECT_THROW(u.convert_to<dna4>(), std::bad_variant_access);
    EXPECT_NO_THROW(u.convert_to<dna5>());
    EXPECT_THROW(u.convert_to<gap>(), std::bad_variant_access);

    dna5 out = u.convert_to<dna5>();
    EXPECT_EQ(out, 'C'_dna5);

    u = gap{};
    gap g = u.convert_unsafely_to<gap>();
    EXPECT_EQ(g, gap{});
}
