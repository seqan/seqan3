// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h>

#include <deque>
#include <list>
#include <ranges>
#include <string>
#include <vector>

#include <seqan3/utility/views/type_reduce.hpp>

// ============================================================================
//  sequential_read
// ============================================================================

template <typename container_t, typename adaptor_t>
void sequential_read(benchmark::State & state)
{
    container_t c;
    c.resize(1'000'000);
    uint8_t i = 0;
    for (auto & e : c)
        e = ++i; // dummy values

    uint8_t dummy = 0;

    if constexpr (std::same_as<adaptor_t, void>)
    {
        for (auto _ : state)
        {
            for (auto e : c)
                dummy += e;
        }
    }
    else
    {
        adaptor_t adaptor;
        auto v = c | adaptor;

        for (auto _ : state)
        {
            for (auto e : v)
                dummy += e;
        }
    }

    [[maybe_unused]] volatile uint8_t dummy2 = dummy;
}

BENCHMARK_TEMPLATE(sequential_read, std::string, void);
BENCHMARK_TEMPLATE(sequential_read, std::string, decltype(std::views::all));
BENCHMARK_TEMPLATE(sequential_read, std::string, decltype(seqan3::views::type_reduce));

BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, void);
BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, decltype(std::views::all));
BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, decltype(seqan3::views::type_reduce));

BENCHMARK_TEMPLATE(sequential_read, std::deque<uint8_t>, void);
BENCHMARK_TEMPLATE(sequential_read, std::deque<uint8_t>, decltype(std::views::all));
BENCHMARK_TEMPLATE(sequential_read, std::deque<uint8_t>, decltype(seqan3::views::type_reduce));

BENCHMARK_TEMPLATE(sequential_read, std::list<uint8_t>, void);
BENCHMARK_TEMPLATE(sequential_read, std::list<uint8_t>, decltype(std::views::all));
BENCHMARK_TEMPLATE(sequential_read, std::list<uint8_t>, decltype(seqan3::views::type_reduce));

// ============================================================================
//  run
// ============================================================================

BENCHMARK_MAIN();
