// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h>

#include <random>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/views/translate.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/seqan2.hpp>
#include <seqan3/utility/range/to.hpp>

#ifdef SEQAN3_HAS_SEQAN2
#    include <seqan/seq_io.h>
#    include <seqan/sequence.h>
#    include <seqan/translation.h>
#endif

// Tags used to define the benchmark type
struct baseline_tag
{}; // Baseline where view is applied and only iterating the output range is benchmarked
struct translate_tag
{}; // Benchmark seqan3::views::translate_single

// ============================================================================
//  sequential_read
// ============================================================================

template <typename rng_t>
void sequential_read_impl(benchmark::State & state, rng_t && rng)
{
    for (auto _ : state)
    {
        for (seqan3::aa27 c : rng)
        {
            auto rank = c.to_rank();
            benchmark::DoNotOptimize(rank);
        }
    }
}

template <typename tag_t>
void sequential_read(benchmark::State & state)
{
    std::vector<seqan3::dna4> dna_sequence{seqan3::test::generate_sequence<seqan3::dna4>(1000, 0, 0)};

    if constexpr (std::is_same_v<tag_t, baseline_tag>)
    {
        seqan3::aa27_vector translated_aa_sequence =
            dna_sequence | seqan3::views::translate_single | seqan3::ranges::to<seqan3::aa27_vector>();
        sequential_read_impl(state, translated_aa_sequence);
    }
    else if constexpr (std::is_same_v<tag_t, translate_tag>)
    {
        auto translated_aa_view = dna_sequence | seqan3::views::translate_single;
        sequential_read_impl(state, translated_aa_view);
    }
}

BENCHMARK_TEMPLATE(sequential_read, baseline_tag);
BENCHMARK_TEMPLATE(sequential_read, translate_tag);

// ============================================================================
//  random_access
// ============================================================================

template <typename rng_t>
void random_access_impl(benchmark::State & state, rng_t && rng, std::vector<size_t> const & access_positions)
{
    for (auto _ : state)
    {
        for (auto pos : access_positions)
        {
            auto access = rng[pos].to_rank();
            benchmark::DoNotOptimize(access);
        }
    }
}

template <typename tag_t>
void random_access(benchmark::State & state)
{
    std::vector<seqan3::dna4> dna_sequence{seqan3::test::generate_sequence<seqan3::dna4>(10000, 0, 0)};

    std::vector<size_t> access_positions{};
    access_positions.resize(200);
    std::mt19937 gen(42);
    std::uniform_int_distribution<size_t> position_generator(0, 1000 - 1);

    for (size_t i = 0; i < access_positions.size(); ++i)
        access_positions[i] = position_generator(gen);

    if constexpr (std::is_same_v<tag_t, baseline_tag>)
    {
        seqan3::aa27_vector translated_aa_sequence =
            dna_sequence | seqan3::views::translate_single | seqan3::ranges::to<seqan3::aa27_vector>();
        random_access_impl(state, translated_aa_sequence, access_positions);
    }
    else
    {
        auto translated_aa_view = dna_sequence | seqan3::views::translate_single;
        random_access_impl(state, translated_aa_view, access_positions);
    }
}

BENCHMARK_TEMPLATE(random_access, baseline_tag);
BENCHMARK_TEMPLATE(random_access, translate_tag);

// ============================================================================
//  copy_vector
// ============================================================================

template <typename adaptor_t>
void copy_impl(benchmark::State & state, std::vector<seqan3::dna4> const & dna_sequence, adaptor_t & adaptor)
{
    for (auto _ : state)
    {
        seqan3::aa27_vector translated_aa_sequence{};
        benchmark::DoNotOptimize(translated_aa_sequence =
                                     dna_sequence | adaptor | seqan3::ranges::to<seqan3::aa27_vector>());
    }
}

#ifdef SEQAN3_HAS_SEQAN2
template <typename tag_t>
void copy_impl_seqan2(benchmark::State & state, seqan2::DnaString const & dna_sequence)
{
    for (auto _ : state)
    {
        seqan2::String<seqan2::AminoAcid> out{};
        seqan2::translate(out, dna_sequence, seqan2::SINGLE_FRAME, seqan2::CANONICAL, tag_t{});
    }
}
#endif // SEQAN3_HAS_SEQAN2

template <typename tag_t>
void copy(benchmark::State & state)
{
    std::vector<seqan3::dna4> seqan3_dna_sequence{seqan3::test::generate_sequence<seqan3::dna4>(1000, 0, 0)};

#ifdef SEQAN3_HAS_SEQAN2
    seqan2::DnaString seqan2_dna_sequence{seqan3::test::generate_sequence_seqan2<seqan2::Dna>(1000, 0, 0)};
#endif // SEQAN3_HAS_SEQAN2

    if constexpr (std::is_same_v<tag_t, translate_tag>)
    {
        auto adaptor = seqan3::views::translate_single;
        copy_impl(state, seqan3_dna_sequence, adaptor);
    }
#ifdef SEQAN3_HAS_SEQAN2
    else
    {
        copy_impl_seqan2<tag_t>(state, seqan2_dna_sequence);
    }
#endif // SEQAN3_HAS_SEQAN2
}

BENCHMARK_TEMPLATE(copy, translate_tag);

#ifdef SEQAN3_HAS_SEQAN2
BENCHMARK_TEMPLATE(copy, seqan2::Serial);
BENCHMARK_TEMPLATE(copy, seqan2::Parallel);
#endif // SEQAN3_HAS_SEQAN2

// ============================================================================
//  run
// ============================================================================

BENCHMARK_MAIN();
