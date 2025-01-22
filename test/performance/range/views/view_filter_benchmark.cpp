// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>

auto const always_true = [](seqan3::dna4 const letter)
{
    return letter.to_rank() < 42u;
};
auto const randomly_true = [](seqan3::dna4 const letter)
{
    return letter.to_rank() > 1u;
};
auto const never_true = [](seqan3::dna4 const letter)
{
    return letter.to_rank() > 42u;
};

void baseline(benchmark::State & state)
{
    std::vector<seqan3::dna4> input = seqan3::test::generate_sequence<seqan3::dna4>(10000, 0, 0);

    for (auto _ : state)
    {
        std::vector<seqan3::dna4> output{};
        for (seqan3::dna4 const letter : input)
            output.push_back(letter);
        output.clear();
    }
}

template <typename predicate_t>
void loop_if(benchmark::State & state, predicate_t predicate)
{
    std::vector<seqan3::dna4> input = seqan3::test::generate_sequence<seqan3::dna4>(10000, 0, 0);

    for (auto _ : state)
    {
        std::vector<seqan3::dna4> output{};
        for (seqan3::dna4 const letter : input)
            if (predicate(letter))
                output.push_back(letter);
        output.clear();
    }
}

template <typename predicate_t>
void loop_view(benchmark::State & state, predicate_t predicate)
{
    std::vector<seqan3::dna4> input = seqan3::test::generate_sequence<seqan3::dna4>(10000, 0, 0);

    for (auto _ : state)
    {
        std::vector<seqan3::dna4> output{};
        for (seqan3::dna4 const letter : input | std::views::filter(predicate))
            output.push_back(letter);
        output.clear();
    }
}

template <typename predicate_t>
void copy_view(benchmark::State & state, predicate_t predicate)
{
    std::vector<seqan3::dna4> input = seqan3::test::generate_sequence<seqan3::dna4>(10000, 0, 0);

    for (auto _ : state)
    {
        std::vector<seqan3::dna4> output{};
        std::ranges::copy(input | std::views::filter(predicate), std::back_inserter(output));
        output.clear();
    }
}

BENCHMARK(baseline);

BENCHMARK_CAPTURE(loop_if, "always_true", always_true);
BENCHMARK_CAPTURE(loop_view, "always_true", always_true);
BENCHMARK_CAPTURE(copy_view, "always_true", always_true);

BENCHMARK_CAPTURE(loop_if, "randomly_true", randomly_true);
BENCHMARK_CAPTURE(loop_view, "randomly_true", randomly_true);
BENCHMARK_CAPTURE(copy_view, "randomly_true", randomly_true);

BENCHMARK_CAPTURE(loop_if, "never_true", never_true);
BENCHMARK_CAPTURE(loop_view, "never_true", never_true);
BENCHMARK_CAPTURE(copy_view, "never_true", never_true);

// ============================================================================
//  run
// ============================================================================

BENCHMARK_MAIN();
