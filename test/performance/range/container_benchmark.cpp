// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

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
        state.counters["alph_size"] = seqan3::alphabet_size<alphabet_t>;
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
BENCHMARK_TEMPLATE(push_back, std::vector, alphabet_variant<char, dna4>, false);

BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, gap, false);
BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, dna4, false);
BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, gapped<dna4>, false);
BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, dna15, false);
BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, aa27, false);
BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, char, false);
BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, alphabet_variant<char, dna4>, false);

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
BENCHMARK_TEMPLATE(push_back, std::vector, alphabet_variant<char, dna4>, true);

BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, gap, true);
BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, dna4, true);
BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, gapped<dna4>, true);
BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, dna15, true);
BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, aa27, true);
BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, char, true);
BENCHMARK_TEMPLATE(push_back, seqan3::bitcompressed_vector, alphabet_variant<char, dna4>, true);
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
        state.counters["alph_size"] = seqan3::alphabet_size<alphabet_t>;
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
BENCHMARK_TEMPLATE(sequential_read, std::vector, alphabet_variant<char, dna4>);

BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, gap);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, dna4);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, gapped<dna4>);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, dna15);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, aa27);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, char);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, alphabet_variant<char, dna4, dna15>);

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
BENCHMARK_TEMPLATE(sequential_read, std::vector, alphabet_variant<char, dna4>, true);

BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, gap, true);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, dna4, true);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, gapped<dna4>, true);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, dna15, true);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, aa27, true);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, char, true);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, alphabet_variant<char, dna4, dna15>, true);

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
        state.counters["alph_size"] = seqan3::alphabet_size<alphabet_t>;
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
BENCHMARK_TEMPLATE(sequential_write, std::vector, alphabet_variant<char, dna4>);

BENCHMARK_TEMPLATE(sequential_write, seqan3::bitcompressed_vector, gap);
BENCHMARK_TEMPLATE(sequential_write, seqan3::bitcompressed_vector, dna4);
BENCHMARK_TEMPLATE(sequential_write, seqan3::bitcompressed_vector, gapped<dna4>);
BENCHMARK_TEMPLATE(sequential_write, seqan3::bitcompressed_vector, dna15);
BENCHMARK_TEMPLATE(sequential_write, seqan3::bitcompressed_vector, aa27);
BENCHMARK_TEMPLATE(sequential_write, seqan3::bitcompressed_vector, char);
BENCHMARK_TEMPLATE(sequential_write, seqan3::bitcompressed_vector, alphabet_variant<char, dna4>);

// ============================================================================
//  run
// ============================================================================

BENCHMARK_MAIN();
