// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <seqan3/alignment/multiple/align_multiple.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>

// Sequence generator functions which can later be moved to
// seqan3/test/include/seqan3/test/performance/sequence_generator.hpp

namespace seqan3::test
{

template <typename alphabet_t>
auto generate_random_sequence_set(size_t const sequence_length,
                                  size_t const set_size,
                                  size_t const sequence_variance = 0)
{
    using sequence_t = decltype(generate_sequence<alphabet_t>());

    std::vector<sequence_t> vec;

    for (unsigned i = 0; i < set_size; ++i)
    {
        //by varying the seed every iteration, we get pseudo random sequences.
        sequence_t seq = generate_sequence<alphabet_t>(sequence_length, sequence_variance, i);
        vec.push_back(seq);
    }

    return vec;
}

template <typename alphabet_t>
auto generate_similar_sequence_set(size_t const sequence_length,
                                   size_t const set_size,
                                   float const mutation_rate)
{
    using sequence_t = decltype(generate_sequence<alphabet_t>());

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distrib_pos(0, sequence_length - 1); // distribution for choosing a random position.
    std::uniform_int_distribution<> distrib_base(0, 3); // distribution for choosing a random base.

    std::vector<sequence_t> vec;
    // generate first sequences
    sequence_t seq = generate_sequence<alphabet_t>(sequence_length, 0, 1);
    vec.push_back(seq);

    size_t mutations = std::round(sequence_length * mutation_rate); // number of mutations per sequence.
    for (size_t i = 0; i < set_size - 1; ++i)
    {
        // mutate the original sequence.
        sequence_t seq_temp = seq;
        for ( size_t i = 0 ; i < mutations ; i++ )
        {
            size_t pos = distrib_pos(gen);
            seq_temp[pos] = seq_temp[pos].assign_rank(distrib_base(gen));
        }
        vec.push_back(seq_temp);
    }

    return vec;
}

} // namespace seqan3::test


static void arguments(benchmark::internal::Benchmark* b)
{
    for (int32_t i : {10, 50})     // i = sequence length.
    {
        for (int32_t n : {20, 50}) // n = number of sequences.
        {
            b->Args({n, i});
        }
    }
}

static void seqan3_msa_similar_sequences(benchmark::State & state)
{
    auto n = static_cast<size_t>(state.range(0));
    assert(n > 1);
    auto i = static_cast<size_t>(state.range(1));
    assert(i > 0);
    // generate set of similar sequences with (0.3 * sequence_length) mutations per sequence.
    auto seqs = seqan3::test::generate_similar_sequence_set<seqan3::dna4>(i, n, 0.3);

    for (auto _ : state)
    {
        auto result = seqan3::align_multiple(seqs);
        benchmark::DoNotOptimize(result);
    }
}

static void seqan3_msa_random_sequences(benchmark::State & state)
{
    auto n = static_cast<size_t>(state.range(0));
    assert(n > 1);
    auto i = static_cast<size_t>(state.range(1));
    assert(i > 0);
    // generate random sequence set.
    auto seqs = seqan3::test::generate_random_sequence_set<seqan3::dna4>(i, n, 0);

    for (auto _ : state)
    {
        auto result = seqan3::align_multiple(seqs);
        benchmark::DoNotOptimize(result);
    }
}

BENCHMARK(seqan3_msa_similar_sequences)->Apply(arguments);
BENCHMARK(seqan3_msa_random_sequences)->Apply(arguments);

BENCHMARK_MAIN();
