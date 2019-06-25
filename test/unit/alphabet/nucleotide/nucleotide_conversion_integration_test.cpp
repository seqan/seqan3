// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <meta/meta.hpp>

#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/core/type_list.hpp>

using namespace seqan3;

template <typename T>
using nucleotide_conversion = ::testing::Test;

using nucleotide_types = type_list<dna4, dna5, dna15, rna4, rna5, rna15>; // needed for some tests
using nucleotide_gtest_types = detail::transfer_template_args_onto_t<nucleotide_types, ::testing::Types>;


TYPED_TEST_CASE(nucleotide_conversion, nucleotide_gtest_types);

// conversion to any other nucleotide type
TYPED_TEST(nucleotide_conversion, explicit_conversion)
{
    meta::for_each(nucleotide_types{}, [&] (auto && nucl) constexpr
    {
        using out_type = std::decay_t<decltype(nucl)>;
        EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('A')), out_type{}.assign_char('A'));
        EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('C')), out_type{}.assign_char('C'));
        EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('G')), out_type{}.assign_char('G'));
        EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('T')), out_type{}.assign_char('T'));
        EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('U')), out_type{}.assign_char('U'));
        EXPECT_EQ(static_cast<out_type>(TypeParam{}.assign_char('T')), out_type{}.assign_char('U'));
    });
}

// conversion to rna/dna of same size
TYPED_TEST(nucleotide_conversion, implicit_conversion)
{
    using other_type = std::conditional_t<std::is_same_v<TypeParam, rna4>,  dna4,
                       std::conditional_t<std::is_same_v<TypeParam, dna4>,  rna4,
                       std::conditional_t<std::is_same_v<TypeParam, rna5>,  dna5,
                       std::conditional_t<std::is_same_v<TypeParam, dna5>,  rna5,
                       std::conditional_t<std::is_same_v<TypeParam, dna15>, rna15,
                       /* must be rna15 */                                  dna15>>>>>;

    // construct
    EXPECT_EQ(other_type{TypeParam{}.assign_char('C')}, other_type{}.assign_char('C'));
    // assign
    other_type l{};
    l = TypeParam{}.assign_char('C');
    EXPECT_EQ(l, other_type{}.assign_char('C'));
}
