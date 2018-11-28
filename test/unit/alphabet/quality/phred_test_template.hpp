// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <numeric>
#include <sstream>
#include <tuple>
#include <utility>

#include <gtest/gtest.h>

#include <range/v3/algorithm/for_each.hpp>
#include <range/v3/view/zip.hpp>

#include <seqan3/alphabet/quality/all.hpp>

using namespace seqan3;

template <typename T>
class phred : public ::testing::Test
{};

TYPED_TEST_CASE_P(phred);

// more elaborate test of assign_char and to_char, basic test is in alphabet_test.cpp
TYPED_TEST_P(phred, conversion_char)
{
    using char_type = typename TypeParam::char_type;
    typename TypeParam::rank_type value_size = (std::is_same_v<TypeParam, phred42>) ? phred42::max_value_size : TypeParam::value_size;

    // phred42: ['!' .. '^'], phred63: ['!' .. '_'], phred68legacy: [';' .. '~']
    std::vector<char_type> input_char(value_size);
    std::iota(input_char.begin(), input_char.end(), TypeParam::offset_char);

    // phred42: ['!' .. 'I', 'J' .. 'J'], phred63: =input_char, phred68legacy: =intput_char
    std::vector<char_type> output_char = input_char;
    if (std::is_same_v<TypeParam, phred42>)
        std::replace_if(output_char.begin(), output_char.end(), [](char_type c){if (c > 'J') return true; return false;}, 'J');
    for (auto [ ch, cm ] : ranges::view::zip(input_char, output_char))
        EXPECT_EQ((assign_char(TypeParam{}, ch)).to_char(), cm);
}

// test assign_phred and to_phred
TYPED_TEST_P(phred, conversion_phred)
{
    using phred_type = typename TypeParam::phred_type;
    phred_type phred_offset = (std::is_same_v<TypeParam, phred68legacy>) ? phred68legacy::offset_phred : 0;
    typename TypeParam::rank_type value_size = (std::is_same_v<TypeParam, phred42>) ? phred42::max_value_size : TypeParam::value_size;

    //in:  phred42: ['!' .. '^'], phred63: ['!' .. '_'], phred68legacy: [';' .. '~']
    std::vector<phred_type> input_phred(value_size);
    std::iota(input_phred.begin(), input_phred.end(), phred_offset);

    //out: phred42: ['!' .. 'I', 'J' .. 'J'], phred63: =input_char, phred68legacy: =intput_char
    std::vector<phred_type> output_phred = input_phred;
    if (std::is_same_v<TypeParam, phred42>)
        std::replace_if(output_phred.begin(), output_phred.end(), [](phred_type pt){if (pt > 41) return true; return false;}, 41);

    for (auto [ ph, pm ] : ranges::view::zip(input_phred, output_phred))
    {
        TypeParam p; p.assign_phred(ph);
        EXPECT_EQ(p.to_phred(), pm);
    }
}

// more elaborate test of assign_rank and to_rank, basic test is in alphabet_test.cpp
TYPED_TEST_P(phred, conversion_rank)
{
    using rank_type = typename TypeParam::rank_type;
    typename TypeParam::rank_type value_size = (std::is_same_v<TypeParam, phred42>) ? phred42::max_value_size : TypeParam::value_size;

    // phred42: [0 .. 61], phred63: [0 .. 62], phred68legacy: [0 .. 67]
    std::vector<rank_type> input_rank(value_size);
    std::iota(input_rank.begin(), input_rank.end(), 0);

    // phred42: [0 .. 40, 41 .. 41], phred63: =input_rank, phred68legacy: =intput_rank
    std::vector<rank_type> output_rank = input_rank;
    if (std::is_same_v<TypeParam, phred42>)
        std::replace_if(output_rank.begin(), output_rank.end(),
            [](rank_type r){if (r >= phred42::value_size) return true; return false;},
            phred42::value_size-1);

    for (auto [ ch, cm ] : ranges::view::zip(input_rank, output_rank))
    {
        TypeParam t;
        t.assign_rank(ch);
        EXPECT_EQ(t.to_rank(), cm);
    }
}

// test provision of data type `phred_type` and phred converter.
TYPED_TEST_P(phred, quality_concept)
{
    TypeParam p{0};
    [[maybe_unused]] typename TypeParam::phred_type p_type = p.to_phred();

    for (typename TypeParam::phred_type i = 0; i < TypeParam::value_size; ++i)
    {
        // rank and phred scores are identical for phred42 and phred63, else the offset is -5
        typename TypeParam::phred_type offset = (!std::is_same_v<TypeParam, phred68legacy>) ? 0 : phred68legacy::offset_phred;
        p.assign_phred(i + offset);
        EXPECT_EQ(i + offset, p.to_phred());
    }
}

REGISTER_TYPED_TEST_CASE_P(phred, conversion_char, conversion_phred, conversion_rank, quality_concept);
