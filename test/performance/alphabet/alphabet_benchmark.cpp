// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <algorithm>
#include <cstring>

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/all.hpp>

using namespace seqan3;

template <Alphabet alphabet_t>
static void assign_char_(benchmark::State& state)
{
    using char_t = alphabet_char_t<alphabet_t>;
    for (auto _ : state)
    {
        alphabet_t chr{};
        for (char_t c = std::numeric_limits<char_t>::min(); c < std::numeric_limits<char_t>::max(); ++c)
        {
            benchmark::DoNotOptimize
            (
                assign_char_to(c, chr)
            );
        }
    }
    state.counters["loop_iterations"] = std::numeric_limits<char_t>::max() - std::numeric_limits<char_t>::min();
}

template <Alphabet alphabet_t>
static void to_char_(benchmark::State& state)
{
    using char_t = alphabet_char_t<alphabet_t>;
    for (auto _ : state)
    {
        alphabet_t chr{};
        for (char_t c = std::numeric_limits<char_t>::min(); c < std::numeric_limits<char_t>::max(); ++c)
        {
            benchmark::DoNotOptimize
            (
                to_char(assign_char_to(c, chr))
            );
        }
    }
    state.counters["loop_iterations"] = std::numeric_limits<char_t>::max() - std::numeric_limits<char_t>::min();
}

template <Semialphabet alphabet_t>
static void assign_rank_(benchmark::State& state)
{
    using rank_t_ = alphabet_rank_t<alphabet_t>;
    using rank_t = std::conditional_t<std::is_same_v<rank_t_, bool>, uint8_t, rank_t_>;
    for (auto _ : state)
    {
        alphabet_t chr{};
        for (rank_t r = std::numeric_limits<rank_t>::min(); r < std::numeric_limits<rank_t>::max(); ++r)
        {
            benchmark::DoNotOptimize
            (
                assign_rank_to(r % alphabet_size<alphabet_t>, chr)
            );
        }
    }
    state.counters["loop_iterations"] = std::numeric_limits<rank_t>::max() - std::numeric_limits<rank_t>::min();
}

template <Semialphabet alphabet_t>
static void to_rank_(benchmark::State& state)
{
    using rank_t_ = alphabet_rank_t<alphabet_t>;
    using rank_t = std::conditional_t<std::is_same_v<rank_t_, bool>, uint8_t, rank_t_>;

    for (auto _ : state)
    {
        alphabet_t chr{};
        for (rank_t r = std::numeric_limits<rank_t>::min(); r < std::numeric_limits<rank_t>::max(); ++r)
        {
            benchmark::DoNotOptimize
            (
                to_rank(assign_rank_to(r % alphabet_size<alphabet_t>, chr))
            );
        }
    }
    state.counters["loop_iterations"] = std::numeric_limits<rank_t>::max() - std::numeric_limits<rank_t>::min();
}

BENCHMARK_TEMPLATE(assign_char_, gap);
BENCHMARK_TEMPLATE(assign_char_, dna4);
BENCHMARK_TEMPLATE(assign_char_, dna5);
BENCHMARK_TEMPLATE(assign_char_, dna15);
BENCHMARK_TEMPLATE(assign_char_, rna15);
BENCHMARK_TEMPLATE(assign_char_, rna4);
BENCHMARK_TEMPLATE(assign_char_, rna5);
BENCHMARK_TEMPLATE(assign_char_, char);
BENCHMARK_TEMPLATE(assign_char_, gapped<dna4>);
BENCHMARK_TEMPLATE(assign_char_, alphabet_variant<gap,dna4,dna5,dna15,rna15,rna4,rna5>);
BENCHMARK_TEMPLATE(assign_char_, alphabet_variant<char, dna4, dna5, dna15>);

BENCHMARK_TEMPLATE(to_char_, gap);
BENCHMARK_TEMPLATE(to_char_, dna4);
BENCHMARK_TEMPLATE(to_char_, dna5);
BENCHMARK_TEMPLATE(to_char_, dna15);
BENCHMARK_TEMPLATE(to_char_, rna15);
BENCHMARK_TEMPLATE(to_char_, rna4);
BENCHMARK_TEMPLATE(to_char_, rna5);
BENCHMARK_TEMPLATE(to_char_, char);
BENCHMARK_TEMPLATE(to_char_, gapped<dna4>);
BENCHMARK_TEMPLATE(to_char_, alphabet_variant<gap,dna4,dna5,dna15,rna15,rna4,rna5>);
BENCHMARK_TEMPLATE(to_char_, alphabet_variant<char, dna4, dna5, dna15>);

BENCHMARK_TEMPLATE(assign_rank_, gap);
BENCHMARK_TEMPLATE(assign_rank_, dna4);
BENCHMARK_TEMPLATE(assign_rank_, dna5);
BENCHMARK_TEMPLATE(assign_rank_, dna15);
BENCHMARK_TEMPLATE(assign_rank_, rna15);
BENCHMARK_TEMPLATE(assign_rank_, rna4);
BENCHMARK_TEMPLATE(assign_rank_, rna5);
BENCHMARK_TEMPLATE(assign_rank_, char);
BENCHMARK_TEMPLATE(assign_rank_, gapped<dna4>);
BENCHMARK_TEMPLATE(assign_rank_, alphabet_variant<gap,dna4,dna5,dna15,rna15,rna4,rna5>);
BENCHMARK_TEMPLATE(assign_rank_, alphabet_variant<char, dna4, dna5, dna15>);

BENCHMARK_TEMPLATE(to_rank_, gap);
BENCHMARK_TEMPLATE(to_rank_, dna4);
BENCHMARK_TEMPLATE(to_rank_, dna5);
BENCHMARK_TEMPLATE(to_rank_, dna15);
BENCHMARK_TEMPLATE(to_rank_, rna15);
BENCHMARK_TEMPLATE(to_rank_, rna4);
BENCHMARK_TEMPLATE(to_rank_, rna5);
BENCHMARK_TEMPLATE(to_rank_, char);
BENCHMARK_TEMPLATE(to_rank_, gapped<dna4>);
BENCHMARK_TEMPLATE(to_rank_, alphabet_variant<gap,dna4,dna5,dna15,rna15,rna4,rna5>);
BENCHMARK_TEMPLATE(to_rank_, alphabet_variant<char, dna4, dna5, dna15>);

BENCHMARK_MAIN();
