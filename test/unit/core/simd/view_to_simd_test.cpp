// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <deque>
#include <tuple>
#include <vector>

#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/bit_manipulation.hpp>
#include <seqan3/core/simd/simd_traits.hpp>
#include <seqan3/core/simd/simd.hpp>
#include <seqan3/core/simd/view_to_simd.hpp>
#include <seqan3/range/container/aligned_allocator.hpp>
#include <seqan3/range/view/view_all.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/concepts>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>

#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/pretty_printing.hpp>
#include <seqan3/test/simd_utility.hpp>

using namespace seqan3;

template <typename t>
class view_to_simd_test : public ::testing::Test
{
public:
    using container_t = std::tuple_element_t<0, t>;
    using simd_t      = std::tuple_element_t<1, t>;
    using allocator_t = std::conditional_t<simd_traits<simd_t>::length == 1,
                                           std::allocator<simd_t>,
                                           aligned_allocator<simd_t, sizeof(simd_t)>>;

    static constexpr size_t padding_value_dna4 = alphabet_size<value_type_t<container_t>>;
    static constexpr size_t padding_value_custom = 8;

    void SetUp()
    {
        sequences.resize(simd_traits<simd_t>::length);
        for (size_t i = 0; i < simd_traits<simd_t>::length; ++i)
        {
            // Generate sequences that end on different boundaries
            size_t l = simd_traits<simd_t>::length * 64 - (i * simd_traits<simd_t>::length) - i;
            std::ranges::copy(test::generate_sequence<value_type_t<container_t>>(l), std::back_inserter(sequences[i]));
        }

        transformed_simd_vec.resize(simd_traits<simd_t>::length * 64, simd::fill<simd_t>(padding_value_dna4));  // longest sequence in set.
        transformed_simd_vec_padded.resize(simd_traits<simd_t>::length * 64, simd::fill<simd_t>(padding_value_custom));

        for (size_t i = 0; i < simd_traits<simd_t>::length; ++i)
        {
            for (size_t j = 0; j < sequences[i].size(); ++j)
            {
                transformed_simd_vec[j][i] = seqan3::to_rank(sequences[i][j]);
                transformed_simd_vec_padded[j][i] = seqan3::to_rank(sequences[i][j]);
            }
        }
    }

    void compare(auto rng, auto cmp)
    {
        auto it_cmp = cmp.begin();
        for (auto & chunk : rng)
        {
            for (auto & vec : chunk)
            {
                SIMD_EQ(vec, *it_cmp);
                ++it_cmp;
            }
        }
    }

    std::vector<container_t> sequences;
    std::vector<simd_t, allocator_t> transformed_simd_vec{};
    std::vector<simd_t, allocator_t> transformed_simd_vec_padded{};

    using view_to_simd_type = detail::view_to_simd<all_view<std::vector<container_t> &>, simd_t>;
};

using test_types = ::testing::Types<std::tuple<std::vector<dna4>, simd_type_t<int8_t>>,
                                    std::tuple<std::vector<dna4>, simd_type_t<int16_t>>,
                                    std::tuple<std::vector<dna4>, simd_type_t<int32_t>>,
                                    std::tuple<std::vector<dna4>, simd_type_t<int64_t>>,
                                    std::tuple<std::vector<dna4>, simd_type_t<uint8_t>>,
                                    std::tuple<std::vector<dna4>, simd_type_t<uint16_t>>,
                                    std::tuple<std::vector<dna4>, simd_type_t<uint32_t>>,
                                    std::tuple<std::vector<dna4>, simd_type_t<uint64_t>>,
                                    std::tuple<std::deque<dna4>, simd_type_t<int8_t>>,
                                    std::tuple<std::deque<dna4>, simd_type_t<int16_t>>,
                                    std::tuple<std::deque<dna4>, simd_type_t<int32_t>>,
                                    std::tuple<std::deque<dna4>, simd_type_t<int64_t>>,
                                    std::tuple<std::deque<dna4>, simd_type_t<uint8_t>>,
                                    std::tuple<std::deque<dna4>, simd_type_t<uint16_t>>,
                                    std::tuple<std::deque<dna4>, simd_type_t<uint32_t>>,
                                    std::tuple<std::deque<dna4>, simd_type_t<uint64_t>>
                                   >;
TYPED_TEST_CASE(view_to_simd_test, test_types);

TEST(view_to_simd, concept_check)
{
    using cmp_type = std::vector<dna4_vector>;
    using test_type = detail::view_to_simd<all_view<cmp_type &>, simd_type_t<int8_t>>;

    using iter_t = decltype(std::ranges::begin(std::declval<test_type &>()));
    EXPECT_TRUE(std::input_iterator<iter_t>);

    EXPECT_EQ(std::ranges::input_range<cmp_type>, std::ranges::input_range<test_type>);
    EXPECT_NE(std::ranges::forward_range<cmp_type>, std::ranges::forward_range<test_type>);
    EXPECT_NE(std::ranges::bidirectional_range<cmp_type>, std::ranges::bidirectional_range<test_type>);
    EXPECT_NE(std::ranges::random_access_range<cmp_type>, std::ranges::random_access_range<test_type>);
    EXPECT_NE(std::ranges::random_access_range<cmp_type>, std::ranges::random_access_range<test_type>);

    EXPECT_EQ(std::ranges::range<cmp_type>, std::ranges::range<test_type>);
    EXPECT_NE(std::ranges::view<cmp_type>, std::ranges::view<test_type>);
    EXPECT_EQ(std::ranges::sized_range<cmp_type>, std::ranges::sized_range<test_type>);
    EXPECT_NE(std::ranges::common_range<cmp_type>, std::ranges::common_range<test_type>);
    EXPECT_NE(const_iterable_range<cmp_type>, const_iterable_range<test_type>);
    EXPECT_NE((std::ranges::output_range<cmp_type, dna4_vector>), (std::ranges::output_range<test_type, dna4_vector>));
}

TEST(view_to_simd, iter_concept)
{
    using cmp_type = std::vector<dna4_vector>;
    using test_type = detail::view_to_simd<all_view<cmp_type &>, simd_type_t<int8_t>>;
    using iter_t = std::ranges::iterator_t<test_type>;
    using sent_t = std::ranges::sentinel_t<test_type>;

    EXPECT_TRUE(std::input_or_output_iterator<iter_t>);
    EXPECT_TRUE(std::input_iterator<iter_t>);
    EXPECT_FALSE(std::forward_iterator<iter_t>);
    EXPECT_FALSE(std::bidirectional_iterator<iter_t>);
    EXPECT_FALSE(std::random_access_iterator<iter_t>);
    EXPECT_FALSE((std::output_iterator<iter_t, decltype(*std::declval<iter_t &>())>));
    EXPECT_TRUE((std::sentinel_for<sent_t, iter_t>));
}

TYPED_TEST(view_to_simd_test, size)
{
    typename TestFixture::view_to_simd_type to_simd_view{this->sequences};

    EXPECT_EQ(to_simd_view.size(), 64u);
}

TYPED_TEST(view_to_simd_test, empty)
{
    typename TestFixture::view_to_simd_type to_simd_view{this->sequences};

    EXPECT_EQ(to_simd_view.empty(), false);
}

TYPED_TEST(view_to_simd_test, iterate_without_padding)
{
    typename TestFixture::view_to_simd_type to_simd_view{this->sequences};
    this->compare(to_simd_view, this->transformed_simd_vec);
}

TYPED_TEST(view_to_simd_test, iterate_with_padding)
{
    typename TestFixture::view_to_simd_type to_simd_view{this->sequences, TestFixture::padding_value_custom};
    this->compare(to_simd_view, this->transformed_simd_vec_padded);
}

TYPED_TEST(view_to_simd_test, adaptor_pipe)
{
    { // without padding
        auto v = this->sequences | view::to_simd<typename TestFixture::simd_t>;
        this->compare(v, this->transformed_simd_vec);
    }

    { // w padding
        auto v = this->sequences | view::to_simd<typename TestFixture::simd_t>(TestFixture::padding_value_custom);
        this->compare(v, this->transformed_simd_vec_padded);
    }

    { // w padding and calling range
        auto v = view::to_simd<typename TestFixture::simd_t>(TestFixture::padding_value_custom)(this->sequences);
        this->compare(v, this->transformed_simd_vec_padded);
    }
}

TYPED_TEST(view_to_simd_test, adaptor_function)
{
    { // without padding
        auto v = view::to_simd<typename TestFixture::simd_t>(this->sequences);
        this->compare(v, this->transformed_simd_vec);
    }

    { // w padding
        auto v = view::to_simd<typename TestFixture::simd_t>(this->sequences, TestFixture::padding_value_custom);
        this->compare(v, this->transformed_simd_vec_padded);
    }
}

TYPED_TEST(view_to_simd_test, empty_sequences)
{
    using simd_t = typename TestFixture::simd_t;
    std::vector<typename TestFixture::container_t> sequences;
    sequences.resize(simd_traits<simd_t>::length);

    auto v = sequences | view::to_simd<simd_t>;
    this->compare(v, std::vector<simd_t, typename TestFixture::allocator_t>{});

    EXPECT_EQ(v.empty(), true);
    EXPECT_EQ(v.size(), 0u);
}

TYPED_TEST(view_to_simd_test, fewer_sequences)
{
    using simd_t = typename TestFixture::simd_t;
    this->sequences.pop_back();

    for (simd_t & vec : this->transformed_simd_vec)
        vec[simd_traits<simd_t>::length - 1] = TestFixture::padding_value_dna4;

    auto v = this->sequences | view::to_simd<simd_t>;
    this->compare(v, this->transformed_simd_vec);

    if constexpr (simd_traits<simd_t>::length > 1)
    {
        EXPECT_EQ(v.empty(), false);
        EXPECT_EQ(v.size(), 64u);
    }
}

TYPED_TEST(view_to_simd_test, fewer_sequences_w_padding)
{
    using simd_t = typename TestFixture::simd_t;
    this->sequences.pop_back();

    for (simd_t & vec : this->transformed_simd_vec_padded)
        vec[simd_traits<simd_t>::length - 1] = TestFixture::padding_value_custom;

    auto v = this->sequences | view::to_simd<simd_t>(TestFixture::padding_value_custom);
    this->compare(v, this->transformed_simd_vec_padded);

    if constexpr (simd_traits<simd_t>::length > 1)
    {
        EXPECT_EQ(v.empty(), false);
        EXPECT_EQ(v.size(), 64u);
    }
}

TYPED_TEST(view_to_simd_test, empty_underlying_range)
{
    using simd_t = typename TestFixture::simd_t;
    std::vector<typename TestFixture::container_t> sequences{};

    auto v = sequences | view::to_simd<simd_t>;
    this->compare(v, std::vector<simd_t, typename TestFixture::allocator_t>{});

    EXPECT_EQ(v.empty(), true);
    EXPECT_EQ(v.size(), 0u);
}

TYPED_TEST(view_to_simd_test, too_many_sequences)
{
    typename TestFixture::container_t cont;
    std::ranges::copy("ACGTACGACT"_dna4, std::back_inserter(cont));
    this->sequences.push_back(cont);

    EXPECT_THROW(typename TestFixture::view_to_simd_type{this->sequences}, std::invalid_argument);
}
