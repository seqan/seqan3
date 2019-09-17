// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/decorator/gap_decorator.hpp>
#include <seqan3/range/views/pseudo_random_access.hpp>
#include <seqan3/range/views/to_char.hpp>
#include <seqan3/std/ranges>

#include "../iterator_test_template.hpp"

using namespace seqan3;

template <typename t>
class pseudo_random_access_test : public ::testing::Test
{};

using testing_types = ::testing::Types<std::vector<dna4>, gap_decorator<std::vector<dna4> const &>>;

TYPED_TEST_CASE(pseudo_random_access_test, testing_types);

TYPED_TEST(pseudo_random_access_test, concepts)
{
    using test_type = decltype(std::declval<TypeParam &>() | views::pseudo_random_access);

    // guaranteed concepts
    EXPECT_TRUE(std::ranges::random_access_range<test_type>);
    EXPECT_TRUE(std::ranges::view<test_type>);
    EXPECT_TRUE(std::ranges::viewable_range<test_type>);

    // preserved concepts
    EXPECT_EQ(std::ranges::sized_range<TypeParam>, std::ranges::sized_range<test_type>);
    EXPECT_EQ(std::ranges::common_range<TypeParam>, std::ranges::common_range<test_type>);
    EXPECT_EQ(std::ranges::contiguous_range<TypeParam>, std::ranges::contiguous_range<test_type>);
    EXPECT_EQ(const_iterable_range<TypeParam>, const_iterable_range<test_type>);
    EXPECT_EQ((std::ranges::output_range<TypeParam, dna4>), (std::ranges::output_range<test_type, dna4>));
}

TYPED_TEST(pseudo_random_access_test, adaptor)
{
    std::vector<dna4> source{'A'_dna4, 'C'_dna4, 'G'_dna4};
    TypeParam test_range{source};

    // pipe notation
    auto v = test_range | views::pseudo_random_access;
    EXPECT_EQ(v | views::to_char | std::ranges::to<std::string>, "ACG");

    // function notation
    auto v2 = views::pseudo_random_access(test_range);
    EXPECT_EQ(v2 | views::to_char | std::ranges::to<std::string>, "ACG");

    // combinability
    auto v3 = test_range | views::pseudo_random_access | std::views::drop(1);
    EXPECT_EQ(v3 | views::to_char | std::ranges::to<std::string>, "CG");
}

// ----------------------------------------------------------------------------
// iterator test
// ----------------------------------------------------------------------------

using view_t = decltype(std::declval<gap_decorator<std::vector<dna4> const &> &>() | views::pseudo_random_access);
using pseudo_random_iterator = std::ranges::iterator_t<view_t>;

template <>
class iterator_fixture<pseudo_random_iterator> : public ::testing::Test
{
    using iterator_tag = std::random_access_iterator_tag;
    static constexpr bool const_iterable = true;

    std::vector<dna4> expected_range{'A'_dna4, 'C'_dna4, 'G'_dna4};

    gap_decorator<std::vector<dna4> const &> gap_deco{expected_range};
    view_t test_range{gap_deco | views::pseudo_random_access};
};

using test_type = ::testing::Types<pseudo_random_iterator>;

INSTANTIATE_TYPED_TEST_CASE_P(pseudo_random_access_view_iterator, iterator_fixture, test_type);
