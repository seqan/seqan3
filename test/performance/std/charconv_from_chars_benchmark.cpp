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
