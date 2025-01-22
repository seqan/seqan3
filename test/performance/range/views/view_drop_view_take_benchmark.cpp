// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h>

#include <deque>
#include <forward_list>
#include <list>
#include <random>
#include <string>
#include <vector>

#include <seqan3/utility/views/single_pass_input.hpp>

// THIS FILE IMPLICITLY TESTS seqan3::views::slice, because that is just drop piped into take

template <typename container_t, typename drop_t, typename take_t, bool single_pass = false>
void sequential_read(benchmark::State & state)
{
    container_t c;
    c.resize(1'003'000);
    uint8_t i = 0;
    for (auto & e : c)
        e = ++i; // dummy values

    uint8_t dummy = 0;

    // if single_pass, add seqan3::views::single_pass_input, otherwise just &
    using single_t = std::conditional_t<single_pass, decltype(c | seqan3::views::single_pass_input), container_t &>;

    if constexpr (std::same_as<drop_t, void>)
    {
        for (auto _ : state)
        {
            single_t s{c};
            for (auto e : s)
                dummy += e;
        }
    }
    else
    {
        auto dr = drop_t{}(1000);
        auto ta = take_t{}(1'000'000);

        for (auto _ : state)
        {
            single_t s{c};
            auto v = s | dr | ta | dr | ta | dr | ta;
            for (auto e : v)
                dummy += e;
        }
    }

    [[maybe_unused]] volatile uint8_t dummy2 = dummy;
}

BENCHMARK_TEMPLATE(sequential_read, std::string, void, void);
BENCHMARK_TEMPLATE(sequential_read, std::string, decltype(std::views::drop), decltype(std::views::take));

BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, void, void);
BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, decltype(std::views::drop), decltype(std::views::take));

BENCHMARK_TEMPLATE(sequential_read, std::deque<uint8_t>, void, void);
BENCHMARK_TEMPLATE(sequential_read, std::deque<uint8_t>, decltype(std::views::drop), decltype(std::views::take));

BENCHMARK_TEMPLATE(sequential_read, std::list<uint8_t>, void, void);
BENCHMARK_TEMPLATE(sequential_read, std::list<uint8_t>, decltype(std::views::drop), decltype(std::views::take));

BENCHMARK_TEMPLATE(sequential_read, std::forward_list<uint8_t>, void, void);
BENCHMARK_TEMPLATE(sequential_read, std::forward_list<uint8_t>, decltype(std::views::drop), decltype(std::views::take));

BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, void, void, true);
BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, decltype(std::views::drop), decltype(std::views::take), true);

BENCHMARK_TEMPLATE(sequential_read, std::forward_list<uint8_t>, void, void, true);
BENCHMARK_TEMPLATE(sequential_read,
                   std::forward_list<uint8_t>,
                   decltype(std::views::drop),
                   decltype(std::views::take),
                   true);

// ============================================================================
//  random access
// ============================================================================

template <typename container_t, typename drop_t, typename take_t>
void random_access(benchmark::State & state)
{
    container_t c;
    c.resize(1'003'000);
    uint8_t i = 0;
    for (auto & e : c)
        e = ++i; // dummy values

    std::vector<size_t> access_positions;
    access_positions.resize(1'000'000);
    std::mt19937 gen(42);
    std::uniform_int_distribution<size_t> dis(0, 998000 - 1);

    for (size_t i = 0; i < 1'000'000; ++i)
        access_positions[i] = dis(gen);

    uint8_t dummy = 0;

    if constexpr (std::same_as<drop_t, void>)
    {
        for (auto _ : state)
        {
            for (auto i : access_positions)
                dummy += c[i];
        }
    }
    else
    {
        auto dr = drop_t{}(1000);
        auto ta = take_t{}(1'000'000);
        auto v = c | dr | ta | dr | ta | dr | ta;

        for (auto _ : state)
        {
            for (auto i : access_positions)
                dummy += v[i];
        }
    }

    [[maybe_unused]] volatile uint8_t dummy2 = dummy;
}

BENCHMARK_TEMPLATE(random_access, std::string, void, void);
BENCHMARK_TEMPLATE(random_access, std::string, decltype(std::views::drop), decltype(std::views::take));

BENCHMARK_TEMPLATE(random_access, std::vector<uint8_t>, void, void);
BENCHMARK_TEMPLATE(random_access, std::vector<uint8_t>, decltype(std::views::drop), decltype(std::views::take));

BENCHMARK_TEMPLATE(random_access, std::deque<uint8_t>, void, void);
BENCHMARK_TEMPLATE(random_access, std::deque<uint8_t>, decltype(std::views::drop), decltype(std::views::take));

// ============================================================================
//  run
// ============================================================================

BENCHMARK_MAIN();
