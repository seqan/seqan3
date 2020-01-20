// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/range/hash.hpp>

template <typename T>
using alphabet_hashing = ::testing::Test;

using test_types = ::testing::Types<seqan3::dna4, seqan3::qualified<seqan3::dna4, seqan3::phred42>, seqan3::gapped<seqan3::dna4>>;

TYPED_TEST_SUITE(alphabet_hashing, test_types, );

TYPED_TEST(alphabet_hashing, hash)
{
    {
        TypeParam t0{};
        std::hash<TypeParam> h{};
        if constexpr (std::same_as<TypeParam, char>)
        {
            for (size_t i = 0; i < seqan3::alphabet_size<TypeParam>/2; ++i)
            {
                seqan3::assign_rank_to(i, t0);
                ASSERT_EQ(h(t0), i);
            }
        }
        else
        {
            for (size_t i = 0; i < seqan3::alphabet_size<TypeParam>; ++i)
            {
                seqan3::assign_rank_to(i, t0);
                ASSERT_EQ(h(t0), i);
            }
        }
    }
    {
        std::vector<TypeParam> text;
        text.reserve(4);
        for (size_t i = 0; i < 4; ++i)
        {
            text.push_back(seqan3::assign_rank_to(0, TypeParam{}));
        }
        std::hash<decltype(text)> h{};
        ASSERT_EQ(h(text), 0u);
    }
    {
        std::hash<TypeParam const> h{};
        if constexpr (std::same_as<TypeParam, char>)
        {
            for (size_t i = 0; i < seqan3::alphabet_size<TypeParam>/2; ++i)
            {
                TypeParam const t0 = seqan3::assign_rank_to(i, TypeParam{});
                ASSERT_EQ(h(t0), i);
            }
        }
        else
        {
            for (size_t i = 0; i < seqan3::alphabet_size<TypeParam>; ++i)
            {
                TypeParam const t0 = seqan3::assign_rank_to(i, TypeParam{});
                ASSERT_EQ(h(t0), i);
            }
        }
    }
    {
        std::vector<TypeParam> const text(4, seqan3::assign_rank_to(0, TypeParam{}));
        std::hash<decltype(text)> h{};
        ASSERT_EQ(h(text), 0u);
    }
}
