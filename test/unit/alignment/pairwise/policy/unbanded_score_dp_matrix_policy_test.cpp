// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2018, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <memory>

#include <seqan3/alignment/pairwise/policy/unbanded_dp_matrix_policy.hpp>
#include <seqan3/core/metafunction/template_inspection.hpp>

template <typename score_alloc_t>
struct policy_mock :
    public seqan3::detail::unbanded_dp_matrix_policy<policy_mock<score_alloc_t>, score_alloc_t>
{
public:

    using base_t = seqan3::detail::unbanded_dp_matrix_policy<policy_mock<score_alloc_t>, score_alloc_t>;

    using base_t::base_t;

    using base_t::allocate_matrix;
    using base_t::active_column;
    using base_t::next_column;

    using base_t::dimension_first_batch;
    using base_t::dimension_second_batch;
    using base_t::score_matrix;

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

TYPED_TEST_CASE_P(unbanded_score_matrix_test);

TYPED_TEST_P(unbanded_score_matrix_test, constructor)
{
    using mock_t = decltype(this->fixture());

    EXPECT_TRUE(std::is_default_constructible_v<mock_t>);
    EXPECT_TRUE(std::is_copy_constructible_v<mock_t>);
    EXPECT_TRUE(std::is_move_constructible_v<mock_t>);
    EXPECT_TRUE(std::is_copy_assignable_v<mock_t>);
    EXPECT_TRUE(std::is_move_assignable_v<mock_t>);
    EXPECT_TRUE(std::is_destructible_v<mock_t>);
}

TYPED_TEST_P(unbanded_score_matrix_test, allocate_matrix)
{
    std::string seq1{"garfieldthecat"};
    std::string seq2{"garfieldthefatcat"};

    auto mock = this->fixture();

    mock.allocate_matrix(seq1, seq2);

    EXPECT_EQ(mock.dimension_first_batch, seq1.size() + 1);
    EXPECT_EQ(mock.dimension_second_batch, seq2.size() + 1);
    EXPECT_EQ(mock.score_matrix.size(), seq2.size() + 1);
}

TYPED_TEST_P(unbanded_score_matrix_test, active_column)
{
    std::string seq1{"garfieldthecat"};
    std::string seq2{"garfieldthefatcat"};

    auto mock = this->fixture();

    mock.allocate_matrix(seq1, seq2);

    auto [sc_col, trace_col] = mock.active_column();

    EXPECT_EQ(seqan3::size(sc_col), seq2.size() + 1);
    EXPECT_TRUE(seqan3::detail::decays_to_ignore_v<decltype(trace_col)>);
}

TYPED_TEST_P(unbanded_score_matrix_test, next_column)
{
    std::string seq1{"garfieldthecat"};
    std::string seq2{"garfieldthefatcat"};

    auto mock = this->fixture();

    mock.allocate_matrix(seq1, seq2);

    // Get the active column
    auto [sc_col, trace_col] = mock.active_column();
    for (auto & sc : sc_col)
        sc = std::tuple{10, -10};

    // Go to next active column and check if active column has moved.
    mock.next_column();
    auto [sc_col2, trace_col2] = mock.active_column();
    for (auto const sc : sc_col2)
        EXPECT_TRUE((sc == std::tuple{10, -10}));
    EXPECT_TRUE(seqan3::detail::decays_to_ignore_v<decltype(trace_col2)>);
}

// ----------------------------------------------------------------------------
// Register and initiate typed test.
// ----------------------------------------------------------------------------

REGISTER_TYPED_TEST_CASE_P(unbanded_score_matrix_test,
                           constructor, allocate_matrix, active_column, next_column);
using test_type = std::tuple<int32_t, int32_t>;
INSTANTIATE_TYPED_TEST_CASE_P(test, unbanded_score_matrix_test, test_type);
