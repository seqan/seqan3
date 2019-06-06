// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <memory>

#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/alignment/pairwise/policy/unbanded_score_trace_dp_matrix_policy.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>

template <typename score_alloc_t, typename trace_alloc_t>
struct policy_mock :
    public seqan3::detail::unbanded_score_trace_dp_matrix_policy<policy_mock<score_alloc_t, trace_alloc_t>,
                                                                 score_alloc_t,
                                                                 trace_alloc_t>
{
public:

    using base_t = seqan3::detail::unbanded_score_trace_dp_matrix_policy<policy_mock<score_alloc_t, trace_alloc_t>,
                                                                         score_alloc_t,
                                                                         trace_alloc_t>;

    using base_t::base_t;

    // Make member functions available for testing
    using base_t::allocate_matrix;
    using base_t::current_column;
    using base_t::go_next_column;

    using base_t::dimension_first_range;
    using base_t::dimension_second_range;
    using base_t::current_column_index;
    using base_t::score_matrix;
    using base_t::trace_matrix;
};

template <typename alloc_types>
struct unbanded_score_trace_test : public ::testing::Test
{
    using score_alloc_t = std::allocator<std::tuple_element_t<0, alloc_types>>;
    using trace_alloc_t = std::allocator<std::tuple_element_t<1, alloc_types>>;

    auto fixture()
    {
        return policy_mock<score_alloc_t, trace_alloc_t>{};
    }
};

using test_type = ::testing::Types<std::tuple<std::tuple<int32_t, int32_t>, seqan3::detail::trace_directions>>;

TYPED_TEST_CASE(unbanded_score_trace_test, test_type);

TYPED_TEST(unbanded_score_trace_test, constructor)
{
    using mock_t = decltype(this->fixture());

    EXPECT_TRUE(std::is_default_constructible_v<mock_t>);
    EXPECT_TRUE(std::is_copy_constructible_v<mock_t>);
    EXPECT_TRUE(std::is_move_constructible_v<mock_t>);
    EXPECT_TRUE(std::is_copy_assignable_v<mock_t>);
    EXPECT_TRUE(std::is_move_assignable_v<mock_t>);
    EXPECT_TRUE(std::is_destructible_v<mock_t>);
}

TYPED_TEST(unbanded_score_trace_test, allocate_matrix)
{
    std::string seq1{"garfieldthecat"};
    std::string seq2{"garfieldthefatcat"};

    auto mock = this->fixture();

    mock.allocate_matrix(seq1, seq2);

    EXPECT_EQ(mock.dimension_first_range, seq1.size() + 1);
    EXPECT_EQ(mock.dimension_second_range, seq2.size() + 1);
    EXPECT_EQ(mock.score_matrix.size(), seq2.size() + 1);
    EXPECT_EQ(mock.trace_matrix.size(), (seq1.size() + 1) * (seq2.size() + 1));
}

TYPED_TEST(unbanded_score_trace_test, current_column)
{
    std::string seq1{"garfieldthecat"};
    std::string seq2{"garfieldthefatcat"};

    auto mock = this->fixture();

    mock.allocate_matrix(seq1, seq2);

    auto zip_view = mock.current_column();

    EXPECT_EQ(std::ranges::size(zip_view), seq2.size() + 1);

    using value_t = seqan3::value_type_t<decltype(zip_view)>;
    EXPECT_EQ(std::tuple_size_v<value_t>, 3u);
    EXPECT_TRUE(std::ranges::SizedRange<decltype(zip_view)>);
    EXPECT_TRUE(std::ranges::BidirectionalRange<decltype(zip_view)>);
}

TYPED_TEST(unbanded_score_trace_test, go_next_column)
{
    std::string seq1{"garfieldthecat"};
    std::string seq2{"garfieldthefatcat"};

    auto mock = this->fixture();

    mock.allocate_matrix(seq1, seq2);

    // Assign to the active column
    auto zip_view = mock.current_column();
    for (auto && entry : zip_view)
    {
        std::get<0>(entry) = std::tuple{10, -10};
        std::get<2>(entry) = seqan3::detail::trace_directions::diagonal;
    }

    // Get the same active column again
    auto zip_view_2 = mock.current_column();
    size_t row_index = 0u;
    for (auto const entry : zip_view_2)
    {
        EXPECT_TRUE((std::get<0>(entry) == std::tuple{10, -10}));
        EXPECT_EQ(std::get<1>(entry).first,  0u);
        EXPECT_EQ(std::get<1>(entry).second,  row_index++);
        EXPECT_TRUE((std::get<2>(entry) == seqan3::detail::trace_directions::diagonal));
    }

    EXPECT_EQ(mock.current_column_index, 0u);

    // Go to next active column and check if active column has moved.
    mock.go_next_column();
    auto zip_view_3 = mock.current_column();

    row_index = 0;
    for (auto const entry : zip_view_3)
    {
        EXPECT_TRUE((std::get<0>(entry) == std::tuple{10, -10}));
        EXPECT_EQ(std::get<1>(entry).first,  1u);
        EXPECT_EQ(std::get<1>(entry).second,  row_index++);
        EXPECT_TRUE((std::get<2>(entry) == seqan3::detail::trace_directions::none));
    }

    EXPECT_EQ(std::ranges::size(zip_view_3), seq2.size() + 1);
}
