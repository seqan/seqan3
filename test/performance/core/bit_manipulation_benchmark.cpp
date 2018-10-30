// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

#include <cstdlib>

#include <benchmark/benchmark.h>

#include <seqan3/core/bit_manipulation.hpp>

template <typename size_type>
static void is_power_of_two_popcount(benchmark::State& state) {
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

static void is_power_of_two_arithmetic(benchmark::State& state) {
    std::srand(0);
    size_t n = 0;
    for (auto _ : state)
    {
        n = std::rand();
        benchmark::DoNotOptimize(n > 0 && (n & (n-1)) == 0);
    }
}
BENCHMARK(is_power_of_two_arithmetic);

static void is_power_of_two_seqan3(benchmark::State& state) {
    std::srand(0);
    size_t n = 0;
    for (auto _ : state)
    {
        n = std::rand();
        benchmark::DoNotOptimize(seqan3::detail::is_power_of_two(n));
    }
}
BENCHMARK(is_power_of_two_seqan3);

static void next_power_of_two_seqan3(benchmark::State& state) {
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
