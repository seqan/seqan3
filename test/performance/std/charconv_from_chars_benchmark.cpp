// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <cstring>
#include <sstream>
#include <string_view>

#include <seqan3/std/charconv>

// -----------------------------------------------------------------------------
// from_char for integral types
// -----------------------------------------------------------------------------
char const * str = "122";

template <typename arithmetic_type>
static void from_char(benchmark::State & state)
{
    arithmetic_type val{};
    size_t sum{};

    for (auto _ : state)
    {
        std::from_chars(&str[0], &str[0] + sizeof(str), val);
        sum += val;
    }

    // prevent complete optimisation
    [[maybe_unused]] volatile auto fin = sum;
}

BENCHMARK_TEMPLATE(from_char, int8_t);
BENCHMARK_TEMPLATE(from_char, uint8_t);
BENCHMARK_TEMPLATE(from_char, int16_t);
BENCHMARK_TEMPLATE(from_char, uint16_t);
BENCHMARK_TEMPLATE(from_char, int32_t);
BENCHMARK_TEMPLATE(from_char, uint32_t);
BENCHMARK_TEMPLATE(from_char, int64_t);
BENCHMARK_TEMPLATE(from_char, uint64_t);

template <typename arithmetic_type>
static void from_stream(benchmark::State & state)
{
    arithmetic_type val{};
    size_t sum{};

    for (auto _ : state)
    {
        std::stringstream ss;
        ss << str;
        ss >> val;
        sum += val;
    }

    // prevent complete optimisation
    [[maybe_unused]] volatile auto fin = sum;
}

BENCHMARK_TEMPLATE(from_stream, int8_t);
BENCHMARK_TEMPLATE(from_stream, uint8_t);
BENCHMARK_TEMPLATE(from_stream, int16_t);
BENCHMARK_TEMPLATE(from_stream, uint16_t);
BENCHMARK_TEMPLATE(from_stream, int32_t);
BENCHMARK_TEMPLATE(from_stream, uint32_t);
BENCHMARK_TEMPLATE(from_stream, int64_t);
BENCHMARK_TEMPLATE(from_stream, uint64_t);

template <typename arithmetic_type>
static void from_atol(benchmark::State & state)
{
    arithmetic_type val{};
    size_t sum{};

    for (auto _ : state)
    {
        val = atol(str);
        sum += val;
    }

    // prevent complete optimisation
    [[maybe_unused]] volatile auto fin = sum;
}

BENCHMARK_TEMPLATE(from_atol, int8_t);
BENCHMARK_TEMPLATE(from_atol, uint8_t);
BENCHMARK_TEMPLATE(from_atol, int16_t);
BENCHMARK_TEMPLATE(from_atol, uint16_t);
BENCHMARK_TEMPLATE(from_atol, int32_t);
BENCHMARK_TEMPLATE(from_atol, uint32_t);
BENCHMARK_TEMPLATE(from_atol, int64_t);
BENCHMARK_TEMPLATE(from_atol, uint64_t);

#if __has_include(<boost/spirit/include/qi.hpp>)

#include <boost/spirit/include/qi.hpp>

template <typename arithmetic_type>
static void from_boost(benchmark::State & state)
{
    arithmetic_type val{};
    size_t sum{};

    for (auto _ : state)
    {
        auto it = &str[0];
        boost::spirit::qi::phrase_parse(it, &str[0] + sizeof(str),
                                        boost::spirit::qi::int_, boost::spirit::ascii::space, val);
        sum += val;
    }

    // prevent complete optimisation
    [[maybe_unused]] volatile auto fin = sum;
}

BENCHMARK_TEMPLATE(from_boost, int8_t);
BENCHMARK_TEMPLATE(from_boost, uint8_t);
BENCHMARK_TEMPLATE(from_boost, int16_t);
BENCHMARK_TEMPLATE(from_boost, uint16_t);
BENCHMARK_TEMPLATE(from_boost, int32_t);
BENCHMARK_TEMPLATE(from_boost, uint32_t);
BENCHMARK_TEMPLATE(from_boost, int64_t);
BENCHMARK_TEMPLATE(from_boost, uint64_t);

#endif // __has_include(<boost/spirit/include/qi.hpp>)

// -----------------------------------------------------------------------------
// from_char for floating point types
// -----------------------------------------------------------------------------
char const * str_float = "122.45e-2";

template <typename arithmetic_type>
static void from_chars_to_float(benchmark::State & state)
{
    arithmetic_type val{};
    size_t sum{};

    for (auto _ : state)
    {
        std::from_chars(&str_float[0], &str_float[0] + sizeof(str_float), val);
        sum += val;
    }

    // prevent complete optimisation
    [[maybe_unused]] volatile auto fin = sum;
}

template <typename arithmetic_type>
static void from_stream_to_float(benchmark::State & state)
{
    arithmetic_type val{};
    size_t sum{};

    for (auto _ : state)
    {
        std::stringstream ss;
        ss << str_float;
        ss >> val;
        sum += val;
    }

    // prevent complete optimisation
    [[maybe_unused]] volatile auto fin = sum;
}

BENCHMARK_TEMPLATE(from_chars_to_float, float);
BENCHMARK_TEMPLATE(from_chars_to_float, double);
BENCHMARK_TEMPLATE(from_stream_to_float, float);
BENCHMARK_TEMPLATE(from_stream_to_float, double);

BENCHMARK_MAIN();
