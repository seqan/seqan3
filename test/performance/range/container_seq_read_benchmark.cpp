// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <deque>
#include <list>
#include <vector>

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/range/container/all.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>

template <typename t>
using sdsl_int_vec = sdsl::int_vector<sizeof(t) * 8>;

template <typename t>
using small_vec = seqan3::small_vector<t, 10'000>;

// ============================================================================
//  sequential_read
// ============================================================================

template <template <typename> typename container_t, typename alphabet_t, bool const_ = false>
void sequential_read(benchmark::State & state)
{
    auto cont_rando = seqan3::test::generate_sequence<alphabet_t>(10'000, 0, 0);
    container_t<alphabet_t> source(cont_rando.begin(), cont_rando.end());

    using source_ref_t = std::conditional_t<const_,
                                            container_t<alphabet_t> const &,
                                            container_t<alphabet_t> &>;

    source_ref_t source_ref{source};

    alphabet_t a;
    for (auto _ : state)
        for (auto && c : source_ref)
            benchmark::DoNotOptimize(a = c);

    state.counters["sizeof"] = sizeof(alphabet_t);
    if constexpr (seqan3::alphabet<alphabet_t>)
        state.counters["alph_size"] = seqan3::alphabet_size<alphabet_t>;
    state.counters["const"] = const_;
}

BENCHMARK_TEMPLATE(sequential_read, std::vector, char);
BENCHMARK_TEMPLATE(sequential_read, std::vector, uint8_t);
BENCHMARK_TEMPLATE(sequential_read, std::vector, uint16_t);
BENCHMARK_TEMPLATE(sequential_read, std::vector, uint32_t);
BENCHMARK_TEMPLATE(sequential_read, std::vector, uint64_t);
BENCHMARK_TEMPLATE(sequential_read, std::vector, seqan3::gap);
BENCHMARK_TEMPLATE(sequential_read, std::vector, seqan3::dna4);
BENCHMARK_TEMPLATE(sequential_read, std::vector, seqan3::gapped<seqan3::dna4>);
BENCHMARK_TEMPLATE(sequential_read, std::vector, seqan3::dna15);
BENCHMARK_TEMPLATE(sequential_read, std::vector, seqan3::aa27);
BENCHMARK_TEMPLATE(sequential_read, std::vector, seqan3::alphabet_variant<char, seqan3::dna4>);

BENCHMARK_TEMPLATE(sequential_read, std::deque, char);
BENCHMARK_TEMPLATE(sequential_read, std::deque, uint8_t);
BENCHMARK_TEMPLATE(sequential_read, std::deque, uint16_t);
BENCHMARK_TEMPLATE(sequential_read, std::deque, uint32_t);
BENCHMARK_TEMPLATE(sequential_read, std::deque, uint64_t);
BENCHMARK_TEMPLATE(sequential_read, std::deque, seqan3::gap);
BENCHMARK_TEMPLATE(sequential_read, std::deque, seqan3::dna4);
BENCHMARK_TEMPLATE(sequential_read, std::deque, seqan3::gapped<seqan3::dna4>);
BENCHMARK_TEMPLATE(sequential_read, std::deque, seqan3::dna15);
BENCHMARK_TEMPLATE(sequential_read, std::deque, seqan3::aa27);
BENCHMARK_TEMPLATE(sequential_read, std::deque, seqan3::alphabet_variant<char, seqan3::dna4>);

BENCHMARK_TEMPLATE(sequential_read, std::list, char);
BENCHMARK_TEMPLATE(sequential_read, std::list, uint8_t);
BENCHMARK_TEMPLATE(sequential_read, std::list, uint16_t);
BENCHMARK_TEMPLATE(sequential_read, std::list, uint32_t);
BENCHMARK_TEMPLATE(sequential_read, std::list, uint64_t);
BENCHMARK_TEMPLATE(sequential_read, std::list, seqan3::gap);
BENCHMARK_TEMPLATE(sequential_read, std::list, seqan3::dna4);
BENCHMARK_TEMPLATE(sequential_read, std::list, seqan3::gapped<seqan3::dna4>);
BENCHMARK_TEMPLATE(sequential_read, std::list, seqan3::dna15);
BENCHMARK_TEMPLATE(sequential_read, std::list, seqan3::aa27);
BENCHMARK_TEMPLATE(sequential_read, std::list, seqan3::alphabet_variant<char, seqan3::dna4>);

BENCHMARK_TEMPLATE(sequential_read, sdsl_int_vec, uint8_t);
BENCHMARK_TEMPLATE(sequential_read, sdsl_int_vec, uint16_t);
BENCHMARK_TEMPLATE(sequential_read, sdsl_int_vec, uint32_t);
BENCHMARK_TEMPLATE(sequential_read, sdsl_int_vec, uint64_t);

BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, char);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, seqan3::gap);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, seqan3::dna4);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, seqan3::gapped<seqan3::dna4>);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, seqan3::dna15);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, seqan3::aa27);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, seqan3::alphabet_variant<char, seqan3::dna4>);

BENCHMARK_TEMPLATE(sequential_read, small_vec, char);
BENCHMARK_TEMPLATE(sequential_read, small_vec, seqan3::gap);
BENCHMARK_TEMPLATE(sequential_read, small_vec, seqan3::dna4);
BENCHMARK_TEMPLATE(sequential_read, small_vec, seqan3::gapped<seqan3::dna4>);
BENCHMARK_TEMPLATE(sequential_read, small_vec, seqan3::dna15);
BENCHMARK_TEMPLATE(sequential_read, small_vec, seqan3::aa27);
BENCHMARK_TEMPLATE(sequential_read, small_vec, seqan3::alphabet_variant<char, seqan3::dna4>);

// ============================================================================
//  sequential_read (const)
// ============================================================================

BENCHMARK_TEMPLATE(sequential_read, std::vector, char, true);
BENCHMARK_TEMPLATE(sequential_read, std::vector, uint8_t, true);
BENCHMARK_TEMPLATE(sequential_read, std::vector, uint16_t, true);
BENCHMARK_TEMPLATE(sequential_read, std::vector, uint32_t, true);
BENCHMARK_TEMPLATE(sequential_read, std::vector, uint64_t, true);
BENCHMARK_TEMPLATE(sequential_read, std::vector, seqan3::gap, true);
BENCHMARK_TEMPLATE(sequential_read, std::vector, seqan3::dna4, true);
BENCHMARK_TEMPLATE(sequential_read, std::vector, seqan3::gapped<seqan3::dna4>, true);
BENCHMARK_TEMPLATE(sequential_read, std::vector, seqan3::dna15, true);
BENCHMARK_TEMPLATE(sequential_read, std::vector, seqan3::aa27, true);
BENCHMARK_TEMPLATE(sequential_read, std::vector, seqan3::alphabet_variant<char, seqan3::dna4>, true);

BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, char, true);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, seqan3::gap, true);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, seqan3::dna4, true);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, seqan3::gapped<seqan3::dna4>, true);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, seqan3::dna15, true);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, seqan3::aa27, true);
BENCHMARK_TEMPLATE(sequential_read, seqan3::bitcompressed_vector, seqan3::alphabet_variant<char, seqan3::dna4>, true);

BENCHMARK_TEMPLATE(sequential_read, small_vec, char, true);
BENCHMARK_TEMPLATE(sequential_read, small_vec, seqan3::gap, true);
BENCHMARK_TEMPLATE(sequential_read, small_vec, seqan3::dna4, true);
BENCHMARK_TEMPLATE(sequential_read, small_vec, seqan3::gapped<seqan3::dna4>, true);
BENCHMARK_TEMPLATE(sequential_read, small_vec, seqan3::dna15, true);
BENCHMARK_TEMPLATE(sequential_read, small_vec, seqan3::aa27, true);
BENCHMARK_TEMPLATE(sequential_read, small_vec, seqan3::alphabet_variant<char, seqan3::dna4>, true);

// ============================================================================
//  run
// ============================================================================

BENCHMARK_MAIN();
