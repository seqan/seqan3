// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/view/view_all.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>

using namespace seqan3;

template <typename T>
class execution_handler : public ::testing::Test
{};

TYPED_TEST_CASE_P(execution_handler);

TYPED_TEST_P(execution_handler, execute_w_lvalue)
{
    constexpr size_t SIZE = 10000;
    std::vector<std::pair<size_t, size_t>> buffer;
    buffer.resize(SIZE);

    TypeParam exec_handler{};

    auto callable = [](size_t const idx, auto && rng1, auto && rng2)
    {
        return std::pair{idx, rng1.size() + rng2.size()};
    };

    std::vector<dna4_vector> set1;
    std::vector<dna4_vector> set2;

    for (unsigned i = 0; i < SIZE; ++i)
    {
        set1.push_back(test::generate_sequence<dna4>(100, 20, i));
        set2.push_back(test::generate_sequence<dna4>(100, 20, i + SIZE));
    }

    size_t pos = 0;

    for (unsigned i = 0; i < SIZE; ++i, ++pos)
    {
        auto v1 = set1[i] | view::all;
        auto v2 = set2[i] | view::all;
        exec_handler.execute(callable, i, v1, v2, [pos, &buffer] (auto && res) { buffer[pos] = std::move(res); });
    }

    exec_handler.wait();

    for (unsigned i = 0; i < SIZE; ++i)
    {
        EXPECT_EQ(buffer[i].first, i);
        EXPECT_EQ(buffer[i].second, set1[i].size() + set2[i].size());
    }
}

TYPED_TEST_P(execution_handler, execute_w_rvalue)
{
    constexpr size_t SIZE = 10000;
    std::vector<std::pair<size_t, size_t>> buffer;
    buffer.resize(SIZE);

    TypeParam exec_handler{};

    auto callable = [](size_t const idx, auto && rng1, auto && rng2)
    {
        return std::pair{idx, rng1.size() + rng2.size()};
    };

    std::vector<dna4_vector> set1;
    std::vector<dna4_vector> set2;

    for (unsigned i = 0; i < SIZE; ++i)
    {
        set1.push_back(test::generate_sequence<dna4>(100, 20, i));
        set2.push_back(test::generate_sequence<dna4>(100, 20, i + SIZE));
    }

    size_t pos = 0;

    for (unsigned i = 0; i < SIZE; ++i, ++pos)
    {
        exec_handler.execute(callable, i, set1[i] | view::all, set2[i] | view::all,
                             [&buffer, pos] (auto && res) { buffer[pos] = std::move(res); });
    }

    exec_handler.wait();

    for (unsigned i = 0; i < SIZE; ++i)
    {
        EXPECT_EQ(buffer[i].first, i);
        EXPECT_EQ(buffer[i].second, set1[i].size() + set2[i].size());
    }
}

REGISTER_TYPED_TEST_CASE_P(execution_handler, execute_w_lvalue, execute_w_rvalue);
