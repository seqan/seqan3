// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <cstdlib>

#include <benchmark/benchmark.h>

#include <seqan3/core/bit_manipulation.hpp>

template <typename size_type>
static void is_power_of_two_popcount(benchmark::State & state) {
    std::srand(0);
    size_type n = 0;
    for (auto _ : state)
    {
        n = std::rand();
        if constexpr(std::is_same_v<size_type, unsigned long long>)
            benchmark::DoNotOptimize(__builtin_popcountll(n) == 1);
        else if constexpr(std::is_same_v<size_type, unsigned long>)
            benchmark::DoNotOptimize(__builtin_popcountl(n) == 1);
        else if constexpr(std::is_same_v<size_type, unsigned>)
            benchmark::DoNotOptimize(__builtin_popcount(n) == 1);
        else
            static_assert(std::is_same_v<size_type, void>, "FAILED");
    }
}
BENCHMARK_TEMPLATE(is_power_of_two_popcount, unsigned);
BENCHMARK_TEMPLATE(is_power_of_two_popcount, unsigned long);
BENCHMARK_TEMPLATE(is_power_of_two_popcount, unsigned long long);

static void is_power_of_two_arithmetic(benchmark::State & state) {
    std::srand(0);
    size_t n = 0;
    for (auto _ : state)
    {
        n = std::rand();
        benchmark::DoNotOptimize(n > 0 && (n & (n-1)) == 0);
    }
}
BENCHMARK(is_power_of_two_arithmetic);

static void is_power_of_two_seqan3(benchmark::State & state) {
    std::srand(0);
    size_t n = 0;
    for (auto _ : state)
    {
        n = std::rand();
        benchmark::DoNotOptimize(seqan3::detail::is_power_of_two(n));
    }
}
BENCHMARK(is_power_of_two_seqan3);

static void next_power_of_two_seqan3(benchmark::State & state) {
    std::srand(0);
    size_t n = 0;
    for (auto _ : state)
    {
        n = std::rand();
        benchmark::DoNotOptimize(seqan3::detail::next_power_of_two(n));
    }
}
BENCHMARK(next_power_of_two_seqan3);

BENCHMARK_MAIN();
