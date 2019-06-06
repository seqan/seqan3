// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <memory>

#include <range/v3/view/bounded.hpp>

#include <seqan3/alignment/pairwise/policy/unbanded_score_dp_matrix_policy.hpp>
#include <seqan3/core/type_traits/basic.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>

template <typename score_alloc_t>
struct policy_mock :
    public seqan3::detail::unbanded_score_dp_matrix_policy<policy_mock<score_alloc_t>, score_alloc_t>
{
public:

    using base_t = seqan3::detail::unbanded_score_dp_matrix_policy<policy_mock<score_alloc_t>, score_alloc_t>;

    using base_t::base_t;

    // Expose member function for testing
    using base_t::allocate_matrix;
    using base_t::current_column;
    using base_t::go_next_column;

    using base_t::dimension_first_range;
    using base_t::dimension_second_range;
    using base_t::score_matrix;
    using base_t::current_column_index;

};

template <typename alloc_type>
struct unbanded_score_matrix_test : public ::testing::Test
{
    using score_alloc_t = std::allocator<alloc_type>;

    auto fixture()
    {
        return policy_mock<score_alloc_t>{};
    }
};

using test_types = std::tuple<int32_t, int32_t>;
TYPED_TEST_CASE(unbanded_score_matrix_test, ::testing::Types<test_types>);

TYPED_TEST(unbanded_score_matrix_test, constructor)
{
    using mock_t = decltype(this->fixture());

    EXPECT_TRUE(std::is_default_constructible_v<mock_t>);
    EXPECT_TRUE(std::is_copy_constructible_v<mock_t>);
    EXPECT_TRUE(std::is_move_constructible_v<mock_t>);
    EXPECT_TRUE(std::is_copy_assignable_v<mock_t>);
    EXPECT_TRUE(std::is_move_assignable_v<mock_t>);
    EXPECT_TRUE(std::is_destructible_v<mock_t>);
}

TYPED_TEST(unbanded_score_matrix_test, allocate_matrix)
{
    std::string seq1{"garfieldthecat"};
    std::string seq2{"garfieldthefatcat"};

    auto mock = this->fixture();

    mock.allocate_matrix(seq1, seq2);

    EXPECT_EQ(mock.dimension_first_range, seq1.size() + 1);
    EXPECT_EQ(mock.dimension_second_range, seq2.size() + 1);
    EXPECT_EQ(mock.score_matrix.size(), seq2.size() + 1);
}

TYPED_TEST(unbanded_score_matrix_test, current_column)
{
    std::string seq1{"garfieldthecat"};
    std::string seq2{"garfieldthefatcat"};

    auto mock = this->fixture();

    mock.allocate_matrix(seq1, seq2);

    auto zip_view = mock.current_column();

    EXPECT_EQ(std::ranges::size(zip_view), seq2.size() + 1);
    EXPECT_TRUE((std::is_same_v<std::remove_reference_t<decltype(std::get<0>(zip_view[0]))>, TypeParam>));
    EXPECT_TRUE(seqan3::detail::decays_to_ignore_v<std::remove_reference_t<decltype(std::get<2>(zip_view[0]))>>);
}

TYPED_TEST(unbanded_score_matrix_test, go_next_column)
{
    std::string seq1{"garfieldthecat"};
    std::string seq2{"garfieldthefatcat"};

    auto mock = this->fixture();

    mock.allocate_matrix(seq1, seq2);

    // Get the active column
    auto zip_view = mock.current_column();
    for (auto && entry : zip_view)
    {
        std::get<0>(entry) = std::tuple{10, -10};
        std::get<2>(entry) = 1;
    }

    EXPECT_EQ(mock.current_column_index, 0u);
    // Go to next active column and check if active column has moved.
    mock.go_next_column();
    auto zip_view_2 = mock.current_column();
    size_t row_index = 0;
    for (auto const entry : zip_view_2)
    {
        EXPECT_TRUE((std::get<0>(entry) == std::tuple{10, -10}));
        EXPECT_EQ(std::get<1>(entry).first, 1u);
        EXPECT_EQ(std::get<1>(entry).second, row_index++);
        EXPECT_TRUE(seqan3::detail::decays_to_ignore_v<std::remove_reference_t<decltype(std::get<2>(entry))>>);
    }
}

TYPED_TEST(unbanded_score_matrix_test, test_concepts)
{
    auto mock = this->fixture();

    // Get the active column
    auto zip_view = mock.current_column();

    EXPECT_TRUE(std::ranges::BidirectionalRange<decltype(zip_view)>);
    EXPECT_TRUE(std::ranges::SizedRange<decltype(zip_view)>);

    auto it = std::ranges::begin(zip_view);
    EXPECT_TRUE(std::BidirectionalIterator<decltype(it)>);

    auto it_e = std::ranges::end(zip_view);
    EXPECT_TRUE(std::BidirectionalIterator<decltype(it_e)>);
}
