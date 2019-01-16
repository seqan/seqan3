// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <algorithm>
#include <cstring>

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/all.hpp>

using namespace seqan3;

template <alphabet_concept alphabet_t>
static void assign_char(benchmark::State& state)
{
    using char_t = underlying_char_t<alphabet_t>;
    for (auto _ : state)
    {
        alphabet_t chr{};
        for (char_t c = std::numeric_limits<char_t>::min(); c < std::numeric_limits<char_t>::max(); ++c)
        {
            benchmark::DoNotOptimize
            (
                assign_char(chr, c)
            );
        }
    }
    state.counters["loop_iterations"] = std::numeric_limits<char_t>::max() - std::numeric_limits<char_t>::min();
}

template <alphabet_concept alphabet_t>
static void to_char(benchmark::State& state)
{
    using char_t = underlying_char_t<alphabet_t>;
    for (auto _ : state)
    {
        alphabet_t chr{};
        for (char_t c = std::numeric_limits<char_t>::min(); c < std::numeric_limits<char_t>::max(); ++c)
        {
            benchmark::DoNotOptimize
            (
                to_char(assign_char(chr, c))
            );
        }
    }
    state.counters["loop_iterations"] = std::numeric_limits<char_t>::max() - std::numeric_limits<char_t>::min();
}

template <semi_alphabet_concept alphabet_t>
static void assign_rank(benchmark::State& state)
{
    using rank_t_ = underlying_rank_t<alphabet_t>;
    using rank_t = std::conditional_t<std::is_same_v<rank_t_, bool>, uint8_t, rank_t_>;
    for (auto _ : state)
    {
        alphabet_t chr{};
        for (rank_t r = std::numeric_limits<rank_t>::min(); r < std::numeric_limits<rank_t>::max(); ++r)
        {
            benchmark::DoNotOptimize
            (
                assign_rank(chr, r % alphabet_size_v<alphabet_t>)
            );
        }
    }
    state.counters["loop_iterations"] = std::numeric_limits<rank_t>::max() - std::numeric_limits<rank_t>::min();
}

template <semi_alphabet_concept alphabet_t>
static void to_rank(benchmark::State& state)
{
    using rank_t_ = underlying_rank_t<alphabet_t>;
    using rank_t = std::conditional_t<std::is_same_v<rank_t_, bool>, uint8_t, rank_t_>;

    for (auto _ : state)
    {
        alphabet_t chr{};
        for (rank_t r = std::numeric_limits<rank_t>::min(); r < std::numeric_limits<rank_t>::max(); ++r)
        {
            benchmark::DoNotOptimize
            (
                to_rank(assign_rank(chr, r % alphabet_size_v<alphabet_t>))
            );
        }
    }
    state.counters["loop_iterations"] = std::numeric_limits<rank_t>::max() - std::numeric_limits<rank_t>::min();
}

BENCHMARK_TEMPLATE(assign_char, gap);
BENCHMARK_TEMPLATE(assign_char, dna4);
BENCHMARK_TEMPLATE(assign_char, dna5);
BENCHMARK_TEMPLATE(assign_char, dna15);
BENCHMARK_TEMPLATE(assign_char, rna15);
BENCHMARK_TEMPLATE(assign_char, rna4);
BENCHMARK_TEMPLATE(assign_char, rna5);
BENCHMARK_TEMPLATE(assign_char, char);
BENCHMARK_TEMPLATE(assign_char, gapped<dna4>);
BENCHMARK_TEMPLATE(assign_char, union_composition<gap,dna4,dna5,dna15,rna15,rna4,rna5>);
BENCHMARK_TEMPLATE(assign_char, union_composition<char, dna4, dna5, dna15>);

BENCHMARK_TEMPLATE(to_char, gap);
BENCHMARK_TEMPLATE(to_char, dna4);
BENCHMARK_TEMPLATE(to_char, dna5);
BENCHMARK_TEMPLATE(to_char, dna15);
BENCHMARK_TEMPLATE(to_char, rna15);
BENCHMARK_TEMPLATE(to_char, rna4);
BENCHMARK_TEMPLATE(to_char, rna5);
BENCHMARK_TEMPLATE(to_char, char);
BENCHMARK_TEMPLATE(to_char, gapped<dna4>);
BENCHMARK_TEMPLATE(to_char, union_composition<gap,dna4,dna5,dna15,rna15,rna4,rna5>);
BENCHMARK_TEMPLATE(to_char, union_composition<char, dna4, dna5, dna15>);

BENCHMARK_TEMPLATE(assign_rank, gap);
BENCHMARK_TEMPLATE(assign_rank, dna4);
BENCHMARK_TEMPLATE(assign_rank, dna5);
BENCHMARK_TEMPLATE(assign_rank, dna15);
BENCHMARK_TEMPLATE(assign_rank, rna15);
BENCHMARK_TEMPLATE(assign_rank, rna4);
BENCHMARK_TEMPLATE(assign_rank, rna5);
BENCHMARK_TEMPLATE(assign_rank, char);
BENCHMARK_TEMPLATE(assign_rank, gapped<dna4>);
BENCHMARK_TEMPLATE(assign_rank, union_composition<gap,dna4,dna5,dna15,rna15,rna4,rna5>);
BENCHMARK_TEMPLATE(assign_rank, union_composition<char, dna4, dna5, dna15>);

BENCHMARK_TEMPLATE(to_rank, gap);
BENCHMARK_TEMPLATE(to_rank, dna4);
BENCHMARK_TEMPLATE(to_rank, dna5);
BENCHMARK_TEMPLATE(to_rank, dna15);
BENCHMARK_TEMPLATE(to_rank, rna15);
BENCHMARK_TEMPLATE(to_rank, rna4);
BENCHMARK_TEMPLATE(to_rank, rna5);
BENCHMARK_TEMPLATE(to_rank, char);
BENCHMARK_TEMPLATE(to_rank, gapped<dna4>);
BENCHMARK_TEMPLATE(to_rank, union_composition<gap,dna4,dna5,dna15,rna15,rna4,rna5>);
BENCHMARK_TEMPLATE(to_rank, union_composition<char, dna4, dna5, dna15>);

BENCHMARK_MAIN();
