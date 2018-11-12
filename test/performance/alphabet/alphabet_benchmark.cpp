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
#include <cstring>

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/all.hpp>

using namespace seqan3;

template <Alphabet alphabet_t>
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

template <Alphabet alphabet_t>
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

template <semi_Alphabet alphabet_t>
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

template <semi_Alphabet alphabet_t>
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
