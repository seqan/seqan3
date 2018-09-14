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

#include <algorithm>
#include <cctype>
#include <cstring>

#include <benchmark/benchmark.h>

#include <seqan3/io/stream/parse_condition.hpp>

using namespace seqan3;

constexpr std::array<char, 1 << 20> arr{};

template <bool stl>
static void simple(benchmark::State& state)
{
    size_t sum = 0;
    size_t i = 0;

    for (auto _ : state)
    {
        i = (i + 1) % (1 << 20);
        if constexpr (stl)
            sum += std::isupper(arr[i]);
        else
            sum += is_upper(arr[i]);
    }

    // prevent complete optimisation
    [[maybe_unused]] volatile size_t fin = sum;
}

BENCHMARK_TEMPLATE(simple, true);
BENCHMARK_TEMPLATE(simple, false);

template <bool stl>
static void combined(benchmark::State& state)
{
    size_t sum = 0;
    size_t i = 0;

    [[maybe_unused]] auto constexpr mycheck = is_punct || is_upper || is_digit;

    for (auto _ : state)
    {
        i = (i + 1) % (1 << 20);
        if constexpr (stl)
            sum += std::ispunct(arr[i]) || std::isupper(arr[i]) || std::isdigit(arr[i]);
        else
            sum += mycheck(arr[i]);
    }

    // prevent complete optimisation
    [[maybe_unused]] volatile size_t fin = sum;
}

BENCHMARK_TEMPLATE(combined, true);
BENCHMARK_TEMPLATE(combined, false);

BENCHMARK_MAIN();
