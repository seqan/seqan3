// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <vector>

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/std/ranges>
#include <seqan3/test/performance/sequence_generator.hpp>

/** This benchmarks evaluates views::filter against other forms of filtering **/

auto const happy  = [] (seqan3::dna4 const c) { return c.to_rank() < 42; }; // always true
auto const rando  = [] (seqan3::dna4 const c) { return c.to_rank() > 1;  }; // true in 50% of the cases
auto const sad    = [] (seqan3::dna4 const c) { return c.to_rank() > 42; }; // never true

void baseline(benchmark::State & state)
{
    std::vector<seqan3::dna4> input = seqan3::test::generate_sequence<seqan3::dna4>(10000, 0, 0);

    for (auto _ : state)
    {
        std::vector<seqan3::dna4> output;
        for (seqan3::dna4 const c : input)
            output.push_back(c);
        output.clear();
    }
}

template <typename fun_t>
void loop_if(benchmark::State & state, fun_t fun)
{
    std::vector<seqan3::dna4> input = seqan3::test::generate_sequence<seqan3::dna4>(10000, 0, 0);

    for (auto _ : state)
    {
        std::vector<seqan3::dna4> output;
        for (seqan3::dna4 const c : input)
            if (fun(c))
                output.push_back(c);
        output.clear();
    }
}

template <typename fun_t>
void loop_view(benchmark::State & state, fun_t fun)
{
    std::vector<seqan3::dna4> input = seqan3::test::generate_sequence<seqan3::dna4>(10000, 0, 0);

    for (auto _ : state)
    {
        std::vector<seqan3::dna4> output;
        for (seqan3::dna4 const c : input | std::views::filter(fun))
            output.push_back(c);
        output.clear();
    }
}

template <typename fun_t>
void copy_view(benchmark::State & state, fun_t fun)
{
    std::vector<seqan3::dna4> input = seqan3::test::generate_sequence<seqan3::dna4>(10000, 0, 0);

    for (auto _ : state)
    {
        std::vector<seqan3::dna4> output;
        std::ranges::copy(input | std::views::filter(fun),
                          std::ranges::back_inserter(output));
        output.clear();
    }
}

BENCHMARK(baseline);

BENCHMARK_CAPTURE(loop_if, "happy", happy);
BENCHMARK_CAPTURE(loop_if, "random", rando);
BENCHMARK_CAPTURE(loop_if, "sad", sad);

BENCHMARK_CAPTURE(loop_view, "happy", happy);
BENCHMARK_CAPTURE(loop_view, "random", rando);
BENCHMARK_CAPTURE(loop_view, "sad", sad);

BENCHMARK_CAPTURE(copy_view, "happy", happy);
BENCHMARK_CAPTURE(copy_view, "random", rando);
BENCHMARK_CAPTURE(copy_view, "sad", sad);

// ============================================================================
//  run
// ============================================================================

BENCHMARK_MAIN();
