// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <variant>

#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/composite/alphabet_variant.hpp>
#include <seqan3/alphabet/gap/gap.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>

#include "../alphabet_constexpr_test_template.hpp"
#include "../alphabet_test_template.hpp"
#include "../semi_alphabet_constexpr_test_template.hpp"
#include "../semi_alphabet_test_template.hpp"

using seqan3::operator""_dna4;
using seqan3::operator""_dna5;
using seqan3::operator""_rna4;
using seqan3::operator""_rna5;

using alphabet_variant_types = ::testing::Types<seqan3::alphabet_variant<seqan3::dna4, seqan3::gap>,
                                                seqan3::alphabet_variant<seqan3::dna4, seqan3::dna5, seqan3::gap>,
                                                seqan3::alphabet_variant<char, seqan3::gap>>;

INSTANTIATE_TYPED_TEST_SUITE_P(alphabet_variant, alphabet_, alphabet_variant_types, );
INSTANTIATE_TYPED_TEST_SUITE_P(alphabet_variant, semi_alphabet_test, alphabet_variant_types, );
INSTANTIATE_TYPED_TEST_SUITE_P(alphabet_variant, alphabet_constexpr, alphabet_variant_types, );
INSTANTIATE_TYPED_TEST_SUITE_P(alphabet_variant, semi_alphabet_constexpr, alphabet_variant_types, );

template <typename T>
using alphabet_variant_test = ::testing::Test;

TYPED_TEST_SUITE(alphabet_variant_test, alphabet_variant_types, );

TEST(alphabet_variant_test, initialise_from_component_alphabet)
{
    seqan3::dna5 l('A'_rna5);

    using alphabet_t = seqan3::alphabet_variant<seqan3::dna4, seqan3::dna5, seqan3::gap>;
    using variant_t = std::variant<seqan3::dna4, seqan3::dna5, seqan3::gap>;

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

    variant_t variant9 = {static_cast<variant_t>(seqan3::gap{})};
    alphabet_t letter9 = {static_cast<alphabet_t>(seqan3::gap{})};

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
    using alphabet_t = seqan3::alphabet_variant<seqan3::dna4, seqan3::dna5, seqan3::gap>;
    using variant_t = std::variant<seqan3::dna4, seqan3::dna5, seqan3::gap>;

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
    using alphabet_t = seqan3::alphabet_variant<seqan3::dna4, seqan3::dna5, seqan3::gap>;
    using variant_t = std::variant<seqan3::dna4, seqan3::dna5, seqan3::gap>;
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

    letter = seqan3::gap{};
    EXPECT_EQ(letter.to_rank(), 9);
}

TEST(alphabet_variant_test, assign_from_component_alphabet_subtype)
{
    using alphabet_t = seqan3::alphabet_variant<seqan3::dna4, seqan3::dna5, seqan3::gap>;
    using variant_t = std::variant<seqan3::dna4, seqan3::dna5, seqan3::gap>;
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
    using alphabet_t = seqan3::alphabet_variant<seqan3::dna4, seqan3::dna5>;

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
    using alphabet_t = seqan3::alphabet_variant<seqan3::dna4, seqan3::dna5>;

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
    using alphabet1_t = seqan3::alphabet_variant<seqan3::dna4, seqan3::dna5, seqan3::gap>;
    using alphabet2_t = seqan3::alphabet_variant<seqan3::gap, seqan3::dna5, seqan3::dna4>;
    using alphabet3_t = seqan3::alphabet_variant<char, seqan3::gap>;

    EXPECT_TRUE((std::is_same_v<seqan3::alphabet_rank_t<alphabet1_t>, uint8_t>));
    EXPECT_TRUE((std::is_same_v<seqan3::alphabet_rank_t<alphabet2_t>, uint8_t>));
    EXPECT_TRUE((std::is_same_v<seqan3::alphabet_rank_t<alphabet3_t>, uint16_t>));
}

TEST(alphabet_variant_test, alphabet_size_)
{
    using alphabet1_t = seqan3::alphabet_variant<seqan3::dna4, seqan3::dna5, seqan3::gap>;
    using alphabet2_t = seqan3::alphabet_variant<seqan3::gap, seqan3::dna5, seqan3::dna4>;
    using alphabet3_t = seqan3::alphabet_variant<char, seqan3::gap>;

    EXPECT_TRUE((std::is_same_v<decltype(alphabet1_t::alphabet_size), const uint8_t>));
    EXPECT_TRUE((std::is_same_v<decltype(alphabet2_t::alphabet_size), const uint8_t>));
    EXPECT_TRUE((std::is_same_v<decltype(alphabet3_t::alphabet_size), const uint16_t>));

    EXPECT_EQ(alphabet1_t::alphabet_size, 10);
    EXPECT_EQ(alphabet2_t::alphabet_size, 10);
    EXPECT_EQ(alphabet3_t::alphabet_size, 257);
}

TEST(alphabet_variant_test, convert_by_index)
{
    seqan3::alphabet_variant<seqan3::dna4, seqan3::dna5, seqan3::gap> u;
    u = 'C'_dna5;

    EXPECT_FALSE(u.is_alternative<0>());
    EXPECT_TRUE(u.is_alternative<1>());
    EXPECT_FALSE(u.is_alternative<2>());

    EXPECT_THROW(u.convert_to<0>(), std::bad_variant_access);
    EXPECT_NO_THROW(u.convert_to<1>());
    EXPECT_THROW(u.convert_to<2>(), std::bad_variant_access);

    seqan3::dna5 out = u.convert_to<1>();
    EXPECT_EQ(out, 'C'_dna5);

    u = seqan3::gap{};

    seqan3::gap g = u.convert_unsafely_to<2>();
    EXPECT_EQ(g, seqan3::gap{});
}

TEST(alphabet_variant_test, convert_by_type)
{
    seqan3::alphabet_variant<seqan3::dna4, seqan3::dna5, seqan3::gap> u;
    u = 'C'_dna5;

    EXPECT_FALSE(u.is_alternative<seqan3::dna4>());
    EXPECT_TRUE(u.is_alternative<seqan3::dna5>());
    EXPECT_FALSE(u.is_alternative<seqan3::gap>());

    EXPECT_THROW(u.convert_to<seqan3::dna4>(), std::bad_variant_access);
    EXPECT_NO_THROW(u.convert_to<seqan3::dna5>());
    EXPECT_THROW(u.convert_to<seqan3::gap>(), std::bad_variant_access);

    seqan3::dna5 out = u.convert_to<seqan3::dna5>();
    EXPECT_EQ(out, 'C'_dna5);

    u = seqan3::gap{};
    seqan3::gap g = u.convert_unsafely_to<seqan3::gap>();
    EXPECT_EQ(g, seqan3::gap{});
}

TEST(alphabet_variant_test, two_different_variants)
{
    seqan3::alphabet_variant<seqan3::dna4, seqan3::gap> l{seqan3::gap{}};
    seqan3::alphabet_variant<seqan3::dna5, seqan3::gap> r{seqan3::gap{}};

    EXPECT_EQ(l, r);

    l = 'A'_dna4;
    r = 'C'_dna5;
    // r = 'A'_dna5 would also not be equal, because dna4 and dna5 are not comparable, only gap and both aren't gap

    EXPECT_NE(l, r);

    seqan3::alphabet_variant<seqan3::rna4, seqan3::gap> r2{'A'_rna4};
    EXPECT_EQ(l, r2); // this works because rna4 and dna4 are implicitly convertible to each other
}
