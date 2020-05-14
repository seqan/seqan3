// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/views/minimiser_hash.hpp>
#include <seqan3/test/performance/naive_minimiser_hash.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/performance/units.hpp>

inline benchmark::Counter bp_per_second(size_t const basepairs)
{
    return benchmark::Counter(basepairs,
                              benchmark::Counter::kIsIterationInvariantRate,
                              benchmark::Counter::OneK::kIs1000);
}

inline seqan3::shape make_gapped_shape(size_t const k)
{
    seqan3::shape shape_{};

    for (size_t i{0}; i < k - 1; ++i)
        shape_.push_back((i + 1) % 2);

    shape_.push_back(1u);
    shape_.push_back(0u);
    return shape_;
}

static void arguments(benchmark::internal::Benchmark* b)
{
    for (int32_t sequence_length : {1'000, 50'000, /*1'000'000*/})
    {
        for (int32_t k : {8, /*16, 24,*/ 30})
        {
            for (int32_t w : {k + 5, k + 10, k + 20})
            {
                b->Args({sequence_length, k, w});
            }
        }
    }
}

enum class method_tag
{
    seqan3_ungapped,
    seqan3_gapped,
    naive
};

template <method_tag tag>
void compute_minimisers(benchmark::State & state)
{
    auto sequence_length = state.range(0);
    size_t k = static_cast<size_t>(state.range(1));
    uint32_t w = static_cast<size_t>(state.range(2));
    assert(sequence_length > 0);
    assert(k > 0);
    assert(w > k);
    auto seq = seqan3::test::generate_sequence<seqan3::dna4>(sequence_length, 0, 0);

    size_t sum{0};

    for (auto _ : state)
    {
        if constexpr (tag == method_tag::naive)
        {
            for (auto h : seq | seqan3::views::naive_minimiser_hash(seqan3::ungapped{static_cast<uint8_t>(k)}, w))
                benchmark::DoNotOptimize(sum += h);
        }
        else if (tag == method_tag::seqan3_ungapped)
        {
            for (auto h : seq | seqan3::views::minimiser_hash(seqan3::ungapped{static_cast<uint8_t>(k)}, window_size{w}))
                benchmark::DoNotOptimize(sum += h);
        }
        else
        {
            for (auto h : seq | seqan3::views::minimiser_hash(make_gapped_shape(k), window_size{w}))
                benchmark::DoNotOptimize(sum += h);
        }
    }

    state.counters["Throughput[bp/s]"] = bp_per_second(sequence_length - k + 1);
}

BENCHMARK_TEMPLATE(compute_minimisers, method_tag::naive)->Apply(arguments);
BENCHMARK_TEMPLATE(compute_minimisers, method_tag::seqan3_ungapped)->Apply(arguments);
BENCHMARK_TEMPLATE(compute_minimisers, method_tag::seqan3_gapped)->Apply(arguments);

BENCHMARK_MAIN();
