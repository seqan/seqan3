// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <cstdlib>

#include <benchmark/benchmark.h>

#include <seqan3/core/platform.hpp> // pre-define SEQAN3_DEPRECATED_HEADER
#pragma push_macro("SEQAN3_DEPRECATED_HEADER")
#undef SEQAN3_DEPRECATED_HEADER
#define SEQAN3_DEPRECATED_HEADER(...)
#include <seqan3/core/bit_manipulation.hpp> // include header without warning
#pragma pop_macro("SEQAN3_DEPRECATED_HEADER")

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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
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
#pragma GCC diagnostic pop

static void is_power_of_two_std(benchmark::State & state) {
    std::srand(0);
    size_t n = 0;
    for (auto _ : state)
    {
        n = std::rand();
        benchmark::DoNotOptimize(std::has_single_bit(n));
    }
}
BENCHMARK(is_power_of_two_std);

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
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
#pragma GCC diagnostic pop

static void next_power_of_two_std(benchmark::State & state) {
    std::srand(0);
    size_t n = 0;
    for (auto _ : state)
    {
        n = std::rand();
        benchmark::DoNotOptimize(std::bit_ceil(n));
    }
}
BENCHMARK(next_power_of_two_std);

BENCHMARK_MAIN();
