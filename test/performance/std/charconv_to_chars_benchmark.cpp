// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h>

#include <cstring>
#include <seqan3/std/charconv>
#include <sstream>
#include <string_view>

// -----------------------------------------------------------------------------
// to_chars for integral types
// -----------------------------------------------------------------------------

template <typename arithmetic_type>
static void to_char(benchmark::State & state)
{
    arithmetic_type val{120};
    std::stringstream ss;
    char buffer[10];

    for (auto _ : state)
    {
        auto result = std::to_chars(std::begin(buffer), std::end(buffer), val);
        *result.ptr = 0;
        ss << buffer;
    }

    // prevent complete optimisation
    [[maybe_unused]] volatile auto fin = ss.rdbuf();
}

BENCHMARK_TEMPLATE(to_char, int8_t);
BENCHMARK_TEMPLATE(to_char, uint8_t);
BENCHMARK_TEMPLATE(to_char, int16_t);
BENCHMARK_TEMPLATE(to_char, uint16_t);
BENCHMARK_TEMPLATE(to_char, int32_t);
BENCHMARK_TEMPLATE(to_char, uint32_t);
BENCHMARK_TEMPLATE(to_char, int64_t);
BENCHMARK_TEMPLATE(to_char, uint64_t);

template <typename arithmetic_type>
static void to_stream(benchmark::State & state)
{
    arithmetic_type val{120};
    std::stringstream ss;

    for (auto _ : state)
    {
        if constexpr (sizeof(arithmetic_type) == 1)
            ss << (int)val;
        else
            ss << val;
    }

    // prevent complete optimisation
    [[maybe_unused]] volatile auto fin = ss.rdbuf();
}

BENCHMARK_TEMPLATE(to_stream, int8_t);
BENCHMARK_TEMPLATE(to_stream, uint8_t);
BENCHMARK_TEMPLATE(to_stream, int16_t);
BENCHMARK_TEMPLATE(to_stream, uint16_t);
BENCHMARK_TEMPLATE(to_stream, int32_t);
BENCHMARK_TEMPLATE(to_stream, uint32_t);
BENCHMARK_TEMPLATE(to_stream, int64_t);
BENCHMARK_TEMPLATE(to_stream, uint64_t);

BENCHMARK_MAIN();
