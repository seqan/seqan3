// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <algorithm>
#include <cctype>
#include <cstring>

#include <benchmark/benchmark.h>

#include <seqan3/core/char_operations/predicate.hpp>
#include <seqan3/test/seqan2.hpp>

#if SEQAN3_HAS_SEQAN2
#include <seqan/stream.h>
#endif

enum class tag
{
    std,
    seqan2,
    seqan2_serial,
    seqan3,
    seqan3_serial
};

constexpr std::array<char, 1 << 20> arr{};

// ============================================================================
//  simple
// ============================================================================

template <tag id>
static void simple(benchmark::State & state)
{
    size_t sum = 0;
    size_t i = 0;

    for (auto _ : state)
    {
        i = (i + 1) % (1 << 20);
        if constexpr (id == tag::std)
            sum += std::isalpha(arr[i]);
        else if constexpr (id == tag::seqan3)
            sum += seqan3::is_alpha(arr[i]);
#if SEQAN3_HAS_SEQAN2
        else if constexpr (id == tag::seqan2)
            sum += seqan::IsAlpha{}(arr[i]);
#endif
    }

    // prevent complete optimisation
    [[maybe_unused]] volatile size_t fin = sum;
}

BENCHMARK_TEMPLATE(simple, tag::std);
BENCHMARK_TEMPLATE(simple, tag::seqan3);
#if SEQAN3_HAS_SEQAN2
BENCHMARK_TEMPLATE(simple, tag::seqan2);
#endif

// ============================================================================
//  combined
// ============================================================================

template <tag id>
static void combined(benchmark::State & state)
{
    size_t sum = 0;
    size_t i = 0;

    for (auto _ : state)
    {
        i = (i + 1) % (1 << 20);
        if constexpr (id == tag::std)
            sum += std::isalpha(arr[i]) || std::isblank(arr[i]) || std::isdigit(arr[i]);
        else if constexpr (id == tag::seqan3)
            sum += (seqan3::is_alpha || seqan3::is_blank || seqan3::is_digit)(arr[i]);
        else if constexpr (id == tag::seqan3_serial)
            sum += seqan3::is_alpha(arr[i]) || seqan3::is_blank(arr[i]) || seqan3::is_digit(arr[i]);
#if SEQAN3_HAS_SEQAN2
        else if constexpr (id == tag::seqan2)
            sum += seqan::OrFunctor<seqan::OrFunctor<seqan::IsAlpha, seqan::IsBlank>, seqan::IsDigit>{}(arr[i]);
        else if constexpr (id == tag::seqan2_serial)
            sum += seqan::IsAlpha{}(arr[i]) || seqan::IsBlank{}(arr[i]) || seqan::IsDigit{}(arr[i]);
#endif
    }

    // prevent complete optimisation
    [[maybe_unused]] volatile size_t fin = sum;
}

BENCHMARK_TEMPLATE(combined, tag::std);
BENCHMARK_TEMPLATE(combined, tag::seqan3);
BENCHMARK_TEMPLATE(combined, tag::seqan3_serial);
#if SEQAN3_HAS_SEQAN2
BENCHMARK_TEMPLATE(combined, tag::seqan2);
BENCHMARK_TEMPLATE(combined, tag::seqan2_serial);
#endif

BENCHMARK_MAIN();
