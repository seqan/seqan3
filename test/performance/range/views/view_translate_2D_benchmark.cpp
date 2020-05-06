// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <random>
#include <vector>

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/views/join.hpp>
#include <seqan3/range/views/to.hpp>
#include <seqan3/range/views/translate.hpp>
#include <seqan3/range/views/translate_join.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/seqan2.hpp>

#ifdef SEQAN3_HAS_SEQAN2
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/translation.h>
#endif

// Tags used to define the benchmark type
struct baseline_tag{}; // Baseline where view is applied and only iterating the output range is benchmarked
struct translate_tag{}; // Benchmark view_translate followed by seqan3::views::join
struct translate_join_tag{}; // Benchmark seqan3::views::translate_join

// ============================================================================
//  sequential_read
// ============================================================================

template <typename rng_t>
void sequential_read_impl(benchmark::State & state, rng_t && rng)
{
    for (auto _ : state)
        for (auto && outer : rng)
            benchmark::DoNotOptimize(outer[0].to_rank());
}

template <typename tag_t>
void sequential_read(benchmark::State & state)
{
    std::vector<std::vector<seqan3::dna4>> dna_sequence_collection{};
    dna_sequence_collection.resize(1000);

    for (size_t i = 0; i < dna_sequence_collection.size(); ++i)
        dna_sequence_collection[i] = seqan3::test::generate_sequence<seqan3::dna4>(100, 0, 0);

    if constexpr (std::is_same_v<tag_t, baseline_tag>)
    {
        std::vector<seqan3::aa27_vector> translated_aa_sequences = dna_sequence_collection
                                                                 | seqan3::views::translate_join
                                                                 | seqan3::views::to<std::vector<seqan3::aa27_vector>>;
        sequential_read_impl(state, translated_aa_sequences);
    }
    else if constexpr (std::is_same_v<tag_t, translate_tag>)
    {
        auto translated_aa_view = dna_sequence_collection | seqan3::views::translate | seqan3::views::join;
        sequential_read_impl(state, translated_aa_view);
    }
    else
    {
        auto translated_aa_view = dna_sequence_collection | seqan3::views::translate_join;
        sequential_read_impl(state, translated_aa_view);
    }
}

BENCHMARK_TEMPLATE(sequential_read, baseline_tag);
BENCHMARK_TEMPLATE(sequential_read, translate_tag);
BENCHMARK_TEMPLATE(sequential_read, translate_join_tag);

// ============================================================================
//  random_access
// ============================================================================

template <typename rng_t>
void random_access_impl(benchmark::State & state, rng_t && rng, std::vector<size_t> const & access_positions)
{
    for (auto _ : state)
        for (auto pos : access_positions)
            benchmark::DoNotOptimize(rng[pos][0].to_rank());
}

template <typename tag_t>
void random_access(benchmark::State & state)
{
    std::vector<std::vector<seqan3::dna4>> dna_sequence_collection{};
    dna_sequence_collection.resize(1000);

    for (size_t i = 0; i < dna_sequence_collection.size(); ++i)
        dna_sequence_collection[i] = seqan3::test::generate_sequence<seqan3::dna4>(100, 0, 0);

    std::vector<size_t> access_positions{};
    access_positions.resize(200);
    std::mt19937 gen(42);
    std::uniform_int_distribution<size_t> position_generator(0, 1000 - 1);

    for (size_t i = 0; i < access_positions.size(); ++i)
        access_positions[i] = position_generator(gen);

    if constexpr (std::is_same_v<tag_t, baseline_tag>)
    {
        std::vector<seqan3::aa27_vector> translated_aa_sequences = dna_sequence_collection
                                                                 | seqan3::views::translate_join
                                                                 | seqan3::views::to<std::vector<seqan3::aa27_vector>>;
        random_access_impl(state, translated_aa_sequences, access_positions);
    }
    else
    {
        auto translated_aa_view = dna_sequence_collection | seqan3::views::translate_join;
        random_access_impl(state, translated_aa_view, access_positions);
    }
}

BENCHMARK_TEMPLATE(random_access, baseline_tag);
BENCHMARK_TEMPLATE(random_access, translate_join_tag);

// ============================================================================
//  run
// ============================================================================

BENCHMARK_MAIN();
