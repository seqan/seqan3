// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h>

#include <algorithm>
#include <vector>

#include <seqan3/test/literal/bytes.hpp>
#include <seqan3/test/performance/units.hpp>

using namespace seqan3::test::literals;

// This test allocates 1Gb data with a string and copies it to another memory location.
static void memcpy_benchmark(benchmark::State & state)
{
    size_t const size = state.range(0);
    std::vector<char> const src(size, '-');
    std::vector<char> dst(size);

    for (auto _ : state)
        std::copy(src.begin(), src.end(), dst.begin());

    state.counters["bytes_per_second"] = seqan3::test::bytes_per_second(size);
}

// Register the function as a benchmark
BENCHMARK(memcpy_benchmark)->RangeMultiplier(2)->Range(16_MiB, 1_GiB); // 16 MiB to 1 GiB data

BENCHMARK_MAIN();
