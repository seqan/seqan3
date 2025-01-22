// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <algorithm>
#include <ranges>
#include <string>

#include <seqan3/core/algorithm/detail/algorithm_executor_blocking.hpp>
#include <seqan3/test/pretty_printing.hpp>
#include <seqan3/utility/views/zip.hpp>

// A dummy algorithm that just counts the number of equal characters in two sequences.
struct dummy_algorithm
{
    template <typename sequences_t, typename callback_t>
    constexpr void operator()(sequences_t && sequence_pairs, callback_t && callback) const
    {
        auto && [first_seq, second_seq] = sequence_pairs;

        size_t count = 0;
        for (auto && [lhs, rhs] : seqan3::views::zip(first_seq, second_seq))
            if (lhs == rhs)
                ++count;

        if (count != 0) // Simulating not to call the callback without a result.
            callback(count);
    }
};

template <typename resource_t>
struct algorithm_type_for_input
{
    using algorithm_input_t = std::ranges::range_value_t<resource_t>;
    using callback_t = std::function<void(size_t)>;
    using type = std::function<void(algorithm_input_t, callback_t)>;
};

template <typename execution_handler_t>
struct algorithm_executor_blocking_test : public ::testing::Test
{
    // Some globally defined test types
    using sequence_pair_t = std::pair<std::string, std::string>;
    using sequence_pairs_t = std::vector<sequence_pair_t>;

    sequence_pair_t sequence_pair{"AACGTACGT", "ATCGTCCGT"}; // Hamming distance is 7
    sequence_pairs_t sequence_pairs{5, sequence_pair};

    // Do not use more than 4 threads if running in parallel
    execution_handler_t execution_handler()
    {
        if constexpr (std::same_as<execution_handler_t, seqan3::detail::execution_handler_sequential>)
            return execution_handler_t{};
        else
            return execution_handler_t{std::min<uint32_t>(4, std::thread::hardware_concurrency())};
    }
};

using testing_types =
    testing::Types<seqan3::detail::execution_handler_sequential, seqan3::detail::execution_handler_parallel>;
TYPED_TEST_SUITE(algorithm_executor_blocking_test, testing_types, );

TYPED_TEST(algorithm_executor_blocking_test, construction)
{
    using algorithm_t = typename algorithm_type_for_input<typename TestFixture::sequence_pairs_t &>::type;
    using executor_t = seqan3::detail::
        algorithm_executor_blocking<typename TestFixture::sequence_pairs_t &, algorithm_t, size_t, TypeParam>;

    EXPECT_FALSE(std::is_default_constructible_v<executor_t>);
    EXPECT_FALSE(std::is_copy_constructible_v<executor_t>);
    EXPECT_TRUE(std::is_move_constructible_v<executor_t>);
    EXPECT_FALSE(std::is_copy_assignable_v<executor_t>);
    EXPECT_TRUE(std::is_move_assignable_v<executor_t>);
}

TYPED_TEST(algorithm_executor_blocking_test, is_eof)
{
    using algorithm_t = typename algorithm_type_for_input<typename TestFixture::sequence_pairs_t &>::type;
    using executor_t = seqan3::detail::
        algorithm_executor_blocking<typename TestFixture::sequence_pairs_t &, algorithm_t, size_t, TypeParam>;

    executor_t exec{this->sequence_pairs, algorithm_t{dummy_algorithm{}}, 0u, this->execution_handler()};
    EXPECT_FALSE(exec.is_eof());
}

TYPED_TEST(algorithm_executor_blocking_test, type_deduction)
{
    using algorithm_t = typename algorithm_type_for_input<typename TestFixture::sequence_pairs_t &>::type;
    seqan3::detail::algorithm_executor_blocking exec{this->sequence_pairs, algorithm_t{dummy_algorithm{}}, 0u};
    EXPECT_FALSE(exec.is_eof());
}

TYPED_TEST(algorithm_executor_blocking_test, next_result)
{
    using algorithm_t = typename algorithm_type_for_input<typename TestFixture::sequence_pairs_t &>::type;
    using executor_t = seqan3::detail::
        algorithm_executor_blocking<typename TestFixture::sequence_pairs_t &, algorithm_t, size_t, TypeParam>;

    executor_t exec{this->sequence_pairs, algorithm_t{dummy_algorithm{}}, 0u, this->execution_handler()};

    EXPECT_EQ(exec.next_result().value(), 7u);
    EXPECT_EQ(exec.next_result().value(), 7u);
    EXPECT_EQ(exec.next_result().value(), 7u);
    EXPECT_EQ(exec.next_result().value(), 7u);
    EXPECT_EQ(exec.next_result().value(), 7u);
    EXPECT_FALSE(static_cast<bool>(exec.next_result()));
}

TYPED_TEST(algorithm_executor_blocking_test, move_assignment)
{
    using algorithm_t = typename algorithm_type_for_input<typename TestFixture::sequence_pairs_t &>::type;
    using executor_t = seqan3::detail::
        algorithm_executor_blocking<typename TestFixture::sequence_pairs_t &, algorithm_t, size_t, TypeParam>;

    executor_t exec{this->sequence_pairs, algorithm_t{dummy_algorithm{}}, 0u, this->execution_handler()};
    executor_t exec_move_assigned{this->sequence_pairs, algorithm_t{dummy_algorithm{}}, 0u, this->execution_handler()};

    exec_move_assigned = std::move(exec);

    EXPECT_EQ(exec_move_assigned.next_result().value(), 7u);
    EXPECT_EQ(exec_move_assigned.next_result().value(), 7u);
    EXPECT_EQ(exec_move_assigned.next_result().value(), 7u);

    executor_t exec_move_constructed{std::move(exec_move_assigned)};
    EXPECT_EQ(exec_move_constructed.next_result().value(), 7u);
    EXPECT_EQ(exec_move_constructed.next_result().value(), 7u);
    EXPECT_FALSE(static_cast<bool>(exec_move_constructed.next_result()));
}

TYPED_TEST(algorithm_executor_blocking_test, lvalue_sequence_pair_view)
{
    auto v = std::views::single(this->sequence_pair);
    using algorithm_t = typename algorithm_type_for_input<decltype(v)>::type;
    using executor_t = seqan3::detail::algorithm_executor_blocking<decltype(v), algorithm_t, size_t, TypeParam>;

    executor_t exec{v, algorithm_t{dummy_algorithm{}}, 0u, this->execution_handler()};
    EXPECT_EQ(exec.next_result().value(), 7u);
    EXPECT_FALSE(static_cast<bool>(exec.next_result()));
}

TYPED_TEST(algorithm_executor_blocking_test, rvalue_sequence_pair_view)
{
    using single_pair_t = decltype(std::views::single(this->sequence_pair));
    using algorithm_t = typename algorithm_type_for_input<single_pair_t>::type;
    using executor_t = seqan3::detail::algorithm_executor_blocking<single_pair_t, algorithm_t, size_t, TypeParam>;

    executor_t exec{std::views::single(this->sequence_pair),
                    algorithm_t{dummy_algorithm{}},
                    0u,
                    this->execution_handler()};
    EXPECT_EQ(exec.next_result().value(), 7u);
    EXPECT_FALSE(static_cast<bool>(exec.next_result()));
}

TYPED_TEST(algorithm_executor_blocking_test, lvalue_sequence_pairs)
{
    using algorithm_t = typename algorithm_type_for_input<typename TestFixture::sequence_pairs_t &>::type;
    using executor_t = seqan3::detail::
        algorithm_executor_blocking<typename TestFixture::sequence_pairs_t &, algorithm_t, size_t, TypeParam>;

    executor_t exec{this->sequence_pairs, algorithm_t{dummy_algorithm{}}, 0u, this->execution_handler()};
    EXPECT_EQ(exec.next_result().value(), 7u);
    EXPECT_EQ(exec.next_result().value(), 7u);
    EXPECT_EQ(exec.next_result().value(), 7u);
    EXPECT_EQ(exec.next_result().value(), 7u);
    EXPECT_EQ(exec.next_result().value(), 7u);
    EXPECT_FALSE(static_cast<bool>(exec.next_result()));
}

TYPED_TEST(algorithm_executor_blocking_test, rvalue_sequence_pairs_view)
{
    using persist_pairs_t = decltype(this->sequence_pairs | std::views::all);
    using algorithm_t = typename algorithm_type_for_input<persist_pairs_t>::type;
    using executor_t = seqan3::detail::algorithm_executor_blocking<persist_pairs_t, algorithm_t, size_t, TypeParam>;

    executor_t exec{this->sequence_pairs | std::views::all,
                    algorithm_t{dummy_algorithm{}},
                    0u,
                    this->execution_handler()};
    EXPECT_EQ(exec.next_result().value(), 7u);
    EXPECT_EQ(exec.next_result().value(), 7u);
    EXPECT_EQ(exec.next_result().value(), 7u);
    EXPECT_EQ(exec.next_result().value(), 7u);
    EXPECT_EQ(exec.next_result().value(), 7u);
    EXPECT_FALSE(static_cast<bool>(exec.next_result()));
}

TYPED_TEST(algorithm_executor_blocking_test, empty_result_bucket)
{
    using algorithm_t = typename algorithm_type_for_input<typename TestFixture::sequence_pairs_t &>::type;
    using executor_t = seqan3::detail::
        algorithm_executor_blocking<typename TestFixture::sequence_pairs_t &, algorithm_t, size_t, TypeParam>;
    this->sequence_pairs[3].first = "";
    executor_t exec{this->sequence_pairs, algorithm_t{dummy_algorithm{}}, 0u, this->execution_handler()};

    EXPECT_EQ(exec.next_result().value(), 7u);
    EXPECT_EQ(exec.next_result().value(), 7u);
    EXPECT_EQ(exec.next_result().value(), 7u);
    EXPECT_EQ(exec.next_result().value(), 7u);
    EXPECT_FALSE(static_cast<bool>(exec.next_result()));
}

//See issue: https://github.com/seqan/seqan3/issues/1801
TEST(algorithm_executor_blocking_test, issue_1801)
{
    std::vector<std::thread::id> thread_ids{}; // Stores the thread ids.
    std::mutex push_mutex{}; // Used to synchronise concurrent push back operation on the thread_ids vector.

    using callback_t = std::function<void(size_t)>;
    std::function algorithm = [&](std::string const & seq, callback_t && callback)
    {
        { // Tell which thread is working.
            std::unique_lock push_lock{push_mutex};
            thread_ids.push_back(std::this_thread::get_id());
        }

        callback(seq.size());
    };

    // The sequence vector.
    std::vector<std::string> sequences{10000, std::string{"sequence"}};

    // Define an parallel algorithm executor.
    using algorithm_t = decltype(algorithm);
    using executor_t = seqan3::detail::algorithm_executor_blocking<std::vector<std::string> &,
                                                                   algorithm_t,
                                                                   size_t,
                                                                   seqan3::detail::execution_handler_parallel>;

    // Allow at most 2 threads and then execute until no results are available anymore.
    static constexpr size_t thread_count = 2u;
    executor_t executor{sequences, algorithm, 0ull, seqan3::detail::execution_handler_parallel{thread_count}};
    auto result = executor.next_result();

    while (result.has_value())
        result = executor.next_result();

    // Expect exactly many ids as sequences were procesed.
    EXPECT_EQ(thread_ids.size(), sequences.size());

    // Sort and remove unique ids.
    std::sort(thread_ids.begin(), thread_ids.end());
    thread_ids.erase(std::unique(thread_ids.begin(), thread_ids.end()), thread_ids.end());

    // Expect at most thread count many ids. Note it can also be fewer threads since it is not guaranteed, that
    // all threads will get a piece of the cake.
    EXPECT_LE(thread_ids.size(), thread_count);
}
