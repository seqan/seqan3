// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h>

#include <bit>
#include <cstdlib>

template <typename size_type>
static void is_power_of_two_popcount(benchmark::State & state)
{
    std::srand(0);
    size_type n = 0;
    for (auto _ : state)
    {
        n = std::rand();
        if constexpr (std::is_same_v<size_type, unsigned long long>)
        {
            bool result = __builtin_popcountll(n) == 1;
            benchmark::DoNotOptimize(result);
        }
        else if constexpr (std::is_same_v<size_type, unsigned long>)
        {
            bool result = __builtin_popcountl(n) == 1;
            benchmark::DoNotOptimize(result);
        }
        else if constexpr (std::is_same_v<size_type, unsigned>)
        {
            bool result = __builtin_popcount(n) == 1;
            benchmark::DoNotOptimize(result);
        }
        else
            static_assert(std::is_same_v<size_type, void>, "FAILED");
    }
}
BENCHMARK_TEMPLATE(is_power_of_two_popcount, unsigned);
BENCHMARK_TEMPLATE(is_power_of_two_popcount, unsigned long);
BENCHMARK_TEMPLATE(is_power_of_two_popcount, unsigned long long);

static void is_power_of_two_arithmetic(benchmark::State & state)
{
    std::srand(0);
    size_t n = 0;
    for (auto _ : state)
    {
        n = std::rand();
        bool result = n > 0 && (n & (n - 1)) == 0;
        benchmark::DoNotOptimize(result);
    }
}
BENCHMARK(is_power_of_two_arithmetic);

static void is_power_of_two_std(benchmark::State & state)
{
    std::srand(0);
    size_t n = 0;
    for (auto _ : state)
    {
        n = std::rand();
        bool result = std::has_single_bit(n);
        benchmark::DoNotOptimize(result);
    }
}
BENCHMARK(is_power_of_two_std);

static void next_power_of_two_std(benchmark::State & state)
{
    std::srand(0);
    size_t n = 0;
    for (auto _ : state)
    {
        n = std::rand();
        bool result = std::bit_ceil(n);
        benchmark::DoNotOptimize(result);
    }
}
BENCHMARK(next_power_of_two_std);

BENCHMARK_MAIN();
