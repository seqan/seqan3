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
