// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <gtest/gtest.h>

#include <seqan3/std/algorithm>

template <typename t>
struct simulated_alignment_test : public ::testing::Test
{
    using score_type = typename t::score_type;

    score_type match = 0;
    score_type mismatch = -1;
    score_type gap = -1;

    std::string first{"abba"};
    std::string second{"baba"};

    t test_data{first, second};

    auto & matrix()
    {
        return test_data.matrix;
    };

    auto & gold_matrix()
    {
        return test_data.gold_matrix;
    }

    auto last_init_column()
    {
        return test_data.last_init_column;
    }
};

TYPED_TEST_SUITE_P(simulated_alignment_test);

TYPED_TEST_P(simulated_alignment_test, linear_alignment)
{
    std::vector<typename TypeParam::score_type> cmp_matrix{};
    cmp_matrix.resize(this->gold_matrix().size());
    auto cmp_it = cmp_matrix.begin();

    // ----------------------------------------------------------------------------
    // Initialise the first column
    // ----------------------------------------------------------------------------

    auto mat_it = this->matrix().begin();
    auto col = *mat_it;
    auto col_it = col.begin();;
    (*col_it).current = 0;
    (*col_it).up = this->gap;
    (*col_it).w_left = this->gap;
    ++col_it;
    size_t step = 0;

    for (; col_it != col.end(); ++col_it, ++step)
    {
        (*col_it).current = (*col_it).up;
        (*col_it).up += this->gap;
        (*col_it).w_left = (*col_it).current + this->gap;
    }

    // dump the computed column.
    cmp_it = std::ranges::transform(col, cmp_it, [] (auto proxy) { return proxy.current; }).out;

    // ----------------------------------------------------------------------------
    // Compute all columns where the first cell is at the top of the matrix
    // ----------------------------------------------------------------------------

    // The inner kernel as lambda to reduce duplicate code.
    auto inner_kernel = [this] (auto col_it, auto const first_idx, auto const second_idx)
    {
        auto proxy = *col_it;
        auto & [current, diagonal, r_left, w_left, up] = proxy;
        current = ((this->first[first_idx] == this->second[second_idx]) ? this->match : this->mismatch) + diagonal;
        current = std::max(current, std::max(r_left, up));
        up = current + this->gap;
        w_left = current + this->gap;
    };

    unsigned first_idx = 0;
    // Go to next column.
    ++mat_it;
    size_t col_id = 0;
    for (; col_id < this->last_init_column(); ++mat_it, ++first_idx, ++col_id)
    {
        unsigned second_idx = 0;
        col = *mat_it;
        auto col_it = col.begin();
        (*col_it).current = (*col_it).r_left;
        (*col_it).up = (*col_it).current + this->gap;
        (*col_it).w_left = (*col_it).r_left + this->gap;
        ++col_it;

        for (; col_it != col.end(); ++col_it, ++second_idx)
            inner_kernel(col_it, first_idx, second_idx);

        // dump the computed column.
        cmp_it = std::ranges::transform(col, cmp_it, [] (auto proxy) { return proxy.current; }).out;
    }

    // ----------------------------------------------------------------------------
    // In banded case compute all columns starting in the middle of the matrix.
    // ----------------------------------------------------------------------------

    for (; col_id < this->first.size(); ++mat_it, ++first_idx, ++col_id)
    {
        unsigned second_idx = col_id - 2; // get correct index of second sequence.
        col = *mat_it;
        auto col_it = col.begin();
        auto proxy = *col_it;
        auto & [current, diagonal, r_left, w_left, up] = proxy;
        current = ((this->first[first_idx] == this->second[second_idx]) ? this->match : this->mismatch) + diagonal;
        current = std::max(current, r_left);
        up = current + this->gap;
        w_left = current + this->gap;
        ++col_it;
        ++second_idx;

        for (; col_it != col.end(); ++col_it, ++second_idx)
            inner_kernel(col_it, first_idx, second_idx);

        // dump the computed column.
        cmp_it = std::ranges::transform(col, cmp_it, [] (auto proxy) { return proxy.current; }).out;
    }

    EXPECT_TRUE(std::equal(cmp_matrix.begin(), cmp_matrix.end(), this->gold_matrix().begin()));
}

REGISTER_TYPED_TEST_SUITE_P(simulated_alignment_test, linear_alignment);
