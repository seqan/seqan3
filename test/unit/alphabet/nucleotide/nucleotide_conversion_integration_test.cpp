// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <meta/meta.hpp>

#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/core/metafunction/template_inspection.hpp>

using namespace seqan3;

template <typename T>
class nucleotide_conversion : public ::testing::Test
{};

using nucleotide_types2 = meta::list<dna4, dna5, dna15, rna4, rna5, rna15>; // needed for some tests
using nucleotide_types = detail::transfer_template_args_onto_t<nucleotide_types2, ::testing::Types>;

TYPED_TEST_CASE(nucleotide_conversion, nucleotide_types);

// conversion to any other nucleotide type
TYPED_TEST(nucleotide_conversion, explicit_conversion)
{
    meta::for_each(nucleotide_types2{}, [&] (auto && nucl) constexpr
    {
        using out_type = std::decay_t<decltype(nucl)>;
        EXPECT_EQ(static_cast<out_type>(TypeParam::A), out_type::A);
        EXPECT_EQ(static_cast<out_type>(TypeParam::C), out_type::C);
        EXPECT_EQ(static_cast<out_type>(TypeParam::G), out_type::G);
        EXPECT_EQ(static_cast<out_type>(TypeParam::T), out_type::T);
        EXPECT_EQ(static_cast<out_type>(TypeParam::U), out_type::U);
        EXPECT_EQ(static_cast<out_type>(TypeParam::T), out_type::U);
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
    EXPECT_EQ(other_type{TypeParam::C}, other_type::C);
    // assign
    other_type l{};
    l = TypeParam::C;
    EXPECT_EQ(l, other_type::C);
}
