// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <algorithm>
#include <cctype>
#include <cstring>

#include <benchmark/benchmark.h>

#include <seqan3/core/char_operations/predicate.hpp>

using namespace seqan3;

constexpr std::array<char, 1 << 20> arr{};

template <bool stl>
static void simple(benchmark::State & state)
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
static void combined(benchmark::State & state)
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
