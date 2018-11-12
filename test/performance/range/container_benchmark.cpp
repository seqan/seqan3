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
#include <cmath>
#include <cstring>

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/range/container/all.hpp>

using namespace seqan3;

template <typename t>
using sdsl_int_vec = sdsl::int_vector<sizeof(t) * 8>;

// ============================================================================
//  push_back
// ============================================================================

template <template <typename> typename container_t, typename alphabet_t, bool reserve>
static void push_back(benchmark::State& state)
{
    container_t<alphabet_t> c;

    if constexpr(reserve)
        c.reserve(1 << 30);

    for (auto _ : state)
    {
        c.push_back(alphabet_t{});
    }

    state.counters["sizeof"] = sizeof(alphabet_t);
    if constexpr (Alphabet<alphabet_t>)
        state.counters["alph_size"] = seqan3::alphabet_size_v<alphabet_t>;
#if 0
    state.counters["preallocated_mem"] = reserve;
#endif
}

BENCHMARK_TEMPLATE(push_back, std::vector, uint8_t, false);
BENCHMARK_TEMPLATE(push_back, std::vector, uint16_t, false);
BENCHMARK_TEMPLATE(push_back, std::vector, uint32_t, false);
BENCHMARK_TEMPLATE(push_back, std::vector, uint64_t, false);

BENCHMARK_TEMPLATE(push_back, sdsl_int_vec, uint8_t, false);
BENCHMARK_TEMPLATE(push_back, sdsl_int_vec, uint16_t, false);
BENCHMARK_TEMPLATE(push_back, sdsl_int_vec, uint32_t, false);
BENCHMARK_TEMPLATE(push_back, sdsl_int_vec, uint64_t, false);

BENCHMARK_TEMPLATE(push_back, std::vector, gap, false);
BENCHMARK_TEMPLATE(push_back, std::vector, dna4, false);
BENCHMARK_TEMPLATE(push_back, std::vector, gapped<dna4>, false);
BENCHMARK_TEMPLATE(push_back, std::vector, dna15, false);
BENCHMARK_TEMPLATE(push_back, std::vector, aa27, false);
BENCHMARK_TEMPLATE(push_back, std::vector, char, false);
BENCHMARK_TEMPLATE(push_back, std::vector, union_composition<char, dna4>, false);

BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, gap, false);
BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, dna4, false);
BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, gapped<dna4>, false);
BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, dna15, false);
BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, aa27, false);
BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, char, false);
BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, union_composition<char, dna4>, false);

// ============================================================================
//  push_back with prealloc
// ============================================================================

#if 0 // no difference right now, disabled by default
BENCHMARK_TEMPLATE(push_back, std::vector, uint8_t, true);
BENCHMARK_TEMPLATE(push_back, std::vector, uint16_t, true);
BENCHMARK_TEMPLATE(push_back, std::vector, uint32_t, true);
BENCHMARK_TEMPLATE(push_back, std::vector, uint64_t, true);

BENCHMARK_TEMPLATE(push_back, sdsl_int_vec, uint8_t, true);
BENCHMARK_TEMPLATE(push_back, sdsl_int_vec, uint16_t, true);
BENCHMARK_TEMPLATE(push_back, sdsl_int_vec, uint32_t, true);
BENCHMARK_TEMPLATE(push_back, sdsl_int_vec, uint64_t, true);

BENCHMARK_TEMPLATE(push_back, std::vector, gap, true);
BENCHMARK_TEMPLATE(push_back, std::vector, dna4, true);
BENCHMARK_TEMPLATE(push_back, std::vector, gapped<dna4>, true);
BENCHMARK_TEMPLATE(push_back, std::vector, dna15, true);
BENCHMARK_TEMPLATE(push_back, std::vector, aa27, true);
BENCHMARK_TEMPLATE(push_back, std::vector, char, true);
BENCHMARK_TEMPLATE(push_back, std::vector, union_composition<char, dna4>, true);

BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, gap, true);
BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, dna4, true);
BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, gapped<dna4>, true);
BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, dna15, true);
BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, aa27, true);
BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, char, true);
BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, union_composition<char, dna4>, true);
#endif

// ============================================================================
//  sequential_read
// ============================================================================

template <template <typename> typename container_t, typename alphabet_t, bool const_ = false>
static void sequential_read(benchmark::State& state)
{
    [[maybe_unused]] container_t<alphabet_t> source;
    [[maybe_unused]] auto const & c_source = source;

    [[maybe_unused]] std::vector<alphabet_t> target;

    source.resize(1 << 20);
    target.resize(1 << 20);

    size_t pos = 0;
    for (auto _ : state)
    {
        for (size_t i = 0; i < 8; ++i)
        {
            pos = (pos + 1) % (1 << 20);
            if constexpr (const_)
                target[pos] = c_source[pos];
            else
                target[pos] = source[pos];
        }
    }

    // prevent complete optimisation
    [[maybe_unused]] volatile alphabet_t target2 = target.back();

    state.counters["sizeof"] = sizeof(alphabet_t);
    if constexpr (Alphabet<alphabet_t>)
        state.counters["alph_size"] = seqan3::alphabet_size_v<alphabet_t>;
    state.counters["const"] = const_;
}

BENCHMARK_TEMPLATE(sequential_read, std::vector, uint8_t);
BENCHMARK_TEMPLATE(sequential_read, std::vector, uint16_t);
BENCHMARK_TEMPLATE(sequential_read, std::vector, uint32_t);
BENCHMARK_TEMPLATE(sequential_read, std::vector, uint64_t);

BENCHMARK_TEMPLATE(sequential_read, sdsl_int_vec, uint8_t);
BENCHMARK_TEMPLATE(sequential_read, sdsl_int_vec, uint16_t);
BENCHMARK_TEMPLATE(sequential_read, sdsl_int_vec, uint32_t);
BENCHMARK_TEMPLATE(sequential_read, sdsl_int_vec, uint64_t);

BENCHMARK_TEMPLATE(sequential_read, std::vector, gap);
BENCHMARK_TEMPLATE(sequential_read, std::vector, dna4);
BENCHMARK_TEMPLATE(sequential_read, std::vector, gapped<dna4>);
BENCHMARK_TEMPLATE(sequential_read, std::vector, dna15);
BENCHMARK_TEMPLATE(sequential_read, std::vector, aa27);
BENCHMARK_TEMPLATE(sequential_read, std::vector, char);
BENCHMARK_TEMPLATE(sequential_read, std::vector, union_composition<char, dna4>);

BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, gap);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, dna4);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, gapped<dna4>);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, dna15);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, aa27);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, char);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, union_composition<char, dna4, dna15>);

// ============================================================================
//  sequential_read (const)
// ============================================================================

BENCHMARK_TEMPLATE(sequential_read, std::vector, uint8_t, true);
BENCHMARK_TEMPLATE(sequential_read, std::vector, uint16_t, true);
BENCHMARK_TEMPLATE(sequential_read, std::vector, uint32_t, true);
BENCHMARK_TEMPLATE(sequential_read, std::vector, uint64_t, true);

BENCHMARK_TEMPLATE(sequential_read, sdsl_int_vec, uint8_t, true);
BENCHMARK_TEMPLATE(sequential_read, sdsl_int_vec, uint16_t, true);
BENCHMARK_TEMPLATE(sequential_read, sdsl_int_vec, uint32_t, true);
BENCHMARK_TEMPLATE(sequential_read, sdsl_int_vec, uint64_t, true);

BENCHMARK_TEMPLATE(sequential_read, std::vector, gap, true);
BENCHMARK_TEMPLATE(sequential_read, std::vector, dna4, true);
BENCHMARK_TEMPLATE(sequential_read, std::vector, gapped<dna4>, true);
BENCHMARK_TEMPLATE(sequential_read, std::vector, dna15, true);
BENCHMARK_TEMPLATE(sequential_read, std::vector, aa27, true);
BENCHMARK_TEMPLATE(sequential_read, std::vector, char, true);
BENCHMARK_TEMPLATE(sequential_read, std::vector, union_composition<char, dna4>, true);

BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, gap, true);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, dna4, true);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, gapped<dna4>, true);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, dna15, true);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, aa27, true);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, char, true);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, union_composition<char, dna4, dna15>, true);

// ============================================================================
//  sequential_write
// ============================================================================

template <template <typename> typename container_t, typename alphabet_t>
static void sequential_write(benchmark::State& state)
{
    [[maybe_unused]] std::vector<alphabet_t> source;
    [[maybe_unused]] container_t<alphabet_t> target;

    source.resize(1 << 20);
    target.resize(1 << 20);

    size_t pos = 0;
    for (auto _ : state)
    {
        for (size_t i = 0; i < 8; ++i)
        {
            pos = (pos + 1) % (1 << 20);
            target[pos] = source[pos];

        }
    }

    // prevent complete optimisation
    [[maybe_unused]] volatile alphabet_t target2 = target.back();

    state.counters["sizeof"] = sizeof(alphabet_t);
    if constexpr (Alphabet<alphabet_t>)
        state.counters["alph_size"] = seqan3::alphabet_size_v<alphabet_t>;
}

BENCHMARK_TEMPLATE(sequential_write, std::vector, uint8_t);
BENCHMARK_TEMPLATE(sequential_write, std::vector, uint16_t);
BENCHMARK_TEMPLATE(sequential_write, std::vector, uint32_t);
BENCHMARK_TEMPLATE(sequential_write, std::vector, uint64_t);

BENCHMARK_TEMPLATE(sequential_write, sdsl_int_vec, uint8_t);
BENCHMARK_TEMPLATE(sequential_write, sdsl_int_vec, uint16_t);
BENCHMARK_TEMPLATE(sequential_write, sdsl_int_vec, uint32_t);
BENCHMARK_TEMPLATE(sequential_write, sdsl_int_vec, uint64_t);

BENCHMARK_TEMPLATE(sequential_write, std::vector, gap);
BENCHMARK_TEMPLATE(sequential_write, std::vector, dna4);
BENCHMARK_TEMPLATE(sequential_write, std::vector, gapped<dna4>);
BENCHMARK_TEMPLATE(sequential_write, std::vector, dna15);
BENCHMARK_TEMPLATE(sequential_write, std::vector, aa27);
BENCHMARK_TEMPLATE(sequential_write, std::vector, char);
BENCHMARK_TEMPLATE(sequential_write, std::vector, union_composition<char, dna4>);

BENCHMARK_TEMPLATE(sequential_write, seqan3::bitcompressed_vector, gap);
BENCHMARK_TEMPLATE(sequential_write, seqan3::bitcompressed_vector, dna4);
BENCHMARK_TEMPLATE(sequential_write, seqan3::bitcompressed_vector, gapped<dna4>);
BENCHMARK_TEMPLATE(sequential_write, seqan3::bitcompressed_vector, dna15);
BENCHMARK_TEMPLATE(sequential_write, seqan3::bitcompressed_vector, aa27);
BENCHMARK_TEMPLATE(sequential_write, seqan3::bitcompressed_vector, char);
BENCHMARK_TEMPLATE(sequential_write, seqan3::bitcompressed_vector, union_composition<char, dna4>);

// ============================================================================
//  run
// ============================================================================

BENCHMARK_MAIN();
