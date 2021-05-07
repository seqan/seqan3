// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <seqan3/std/algorithm>
#include <seqan3/std/iterator>
#include <seqan3/std/ranges>
#include <thread>
#include <utility>
#include <vector>

#include <gtest/gtest.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/utility/views/chunk.hpp>
#include <seqan3/utility/views/zip.hpp>

//Forward declaration
namespace seqan3::detail
{

struct execution_handler_sequential;

}

template <typename execution_handler_t>
struct execution_handler : public ::testing::Test
{
    static constexpr size_t total_size = 10000;

    void SetUp()
    {
        for (unsigned i = 0; i < total_size; ++i)
        {
            sequence_collection1.push_back(seqan3::test::generate_sequence<seqan3::dna4>(100, 20, i));
            sequence_collection2.push_back(seqan3::test::generate_sequence<seqan3::dna4>(100, 20, i + total_size));
        }
    }

    template <typename buffer_t>
    void check_result(buffer_t const & buffer) const
    {
        for (unsigned i = 0; i < total_size; ++i)
        {
            EXPECT_EQ(buffer[i].first, i) << "Position: " << i;
            EXPECT_EQ(buffer[i].second,
                      sequence_collection1[i].size() + sequence_collection2[i].size()) << "Position: " << i;
        }
    }

    std::vector<seqan3::dna4_vector> sequence_collection1{};
    std::vector<seqan3::dna4_vector> sequence_collection2{};

    // Do not use more than 4 threads if running in parallel
    execution_handler_t execution_helper()
    {
        if constexpr(std::same_as<execution_handler_t, seqan3::detail::execution_handler_sequential>)
            return execution_handler_t{};
        else
            return execution_handler_t{std::min<uint32_t>(4, std::thread::hardware_concurrency())};
    }
};

auto simulate_alignment_with_range = [] (auto indexed_sequence_pairs,
                                         std::function<void(std::pair<size_t, size_t>)> && callback)
{
    for (auto && [sequence_pair, idx] : indexed_sequence_pairs)
        callback(std::pair{idx, std::get<0>(sequence_pair).size() + std::get<1>(sequence_pair).size()});
};

TYPED_TEST_SUITE_P(execution_handler);

TYPED_TEST_P(execution_handler, execute_as_indexed_sequence_pairs)
{
    std::vector<std::pair<size_t, size_t>> buffer;
    buffer.resize(this->total_size);

    TypeParam exec_handler{this->execution_helper()};

    size_t pos = 0;
    size_t chunk_size = 4; // total_size is a multiple of chunk size.

    auto indexed_sequence_pairs = seqan3::views::zip(seqan3::views::zip(this->sequence_collection1,
                                                                        this->sequence_collection2),
                                                     std::views::iota(0));
    using range_iterator_t = std::ranges::iterator_t<decltype(indexed_sequence_pairs)>;

    for (range_iterator_t it = indexed_sequence_pairs.begin();
         it != indexed_sequence_pairs.end();
         it += chunk_size, pos += chunk_size)
    {
        std::ranges::subrange<range_iterator_t, range_iterator_t> chunk{it, std::next(it, chunk_size)};
        exec_handler.execute(simulate_alignment_with_range, std::move(chunk), [pos, &buffer] (auto && res) mutable
        {
            *(buffer.begin() + pos) = std::forward<decltype(res)>(res);
            ++pos;
        });
    }

    exec_handler.wait();

    this->check_result(buffer);
}

TYPED_TEST_P(execution_handler, bulk_execute)
{
    std::vector<std::pair<size_t, size_t>> buffer{};
    buffer.resize(this->total_size);

    TypeParam exec_handler{this->execution_helper()};

    size_t chunk_size = 4;
    auto indexed_sequence_pairs = seqan3::views::zip(seqan3::views::zip(this->sequence_collection1,
                                                                        this->sequence_collection2),
                                                     std::views::iota(0))
                                | seqan3::views::chunk(chunk_size);

    exec_handler.bulk_execute(simulate_alignment_with_range, indexed_sequence_pairs, [&] (auto && result_pair)
    {
        buffer[result_pair.first] = result_pair;
    });
    this->check_result(buffer);
}

REGISTER_TYPED_TEST_SUITE_P(execution_handler, execute_as_indexed_sequence_pairs, bulk_execute);
