// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <seqan3/range/views/zip.hpp>
#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>

inline benchmark::Counter hashes_per_second(size_t const count)
{
    return benchmark::Counter(count,
                              benchmark::Counter::kIsIterationInvariantRate,
                              benchmark::Counter::OneK::kIs1000);
}

static void arguments(benchmark::internal::Benchmark* b)
{
    for (int32_t bins : {64, 8192})
    {
        for (int32_t bits = 1<<15; bits <= 1<<20/* Increase for more extensive benchmarks*/; bits <<= 5)
        {
            for (int32_t hash_num = 2; hash_num < 3/* Increase for more extensive benchmarks*/; ++hash_num)
            {
                b->Args({bins, bits/bins, hash_num, 1'000/* Increase for more extensive benchmarks*/});
            }
        }
    }
}

template <typename ibf_type>
auto set_up(size_t bins, size_t bits, size_t hash_num, size_t sequence_length)
{
    auto bin_indices = seqan3::test::generate_numeric_sequence<size_t>(sequence_length, 0u, bins - 1);
    auto hash_values = seqan3::test::generate_numeric_sequence<size_t>(sequence_length);
    seqan3::interleaved_bloom_filter tmp_ibf(seqan3::bin_count{bins},
                                             seqan3::bin_size{bits},
                                             seqan3::hash_function_count{hash_num});

    ibf_type ibf{std::move(tmp_ibf)};

    return std::make_tuple(bin_indices, hash_values, ibf);
}

template <typename ibf_type>
void emplace_benchmark(::benchmark::State & state)
{
    auto && [ bin_indices, hash_values, ibf ] = set_up<ibf_type>(state.range(0),
                                                                 state.range(1),
                                                                 state.range(2),
                                                                 state.range(3));

    for (auto _ : state)
    {
        for (auto [hash, bin] : seqan3::views::zip(hash_values, bin_indices))
            ibf.emplace(hash, seqan3::bin_index{bin});
    }

    state.counters["hashes/sec"] = hashes_per_second(std::ranges::size(hash_values));
}

template <typename ibf_type>
void bulk_contains_benchmark(::benchmark::State & state)
{
    auto && [ bin_indices, hash_values, ibf ] = set_up<ibf_type>(state.range(0),
                                                                 state.range(1),
                                                                 state.range(2),
                                                                 state.range(3));
    (void) bin_indices;

    for (auto _ : state)
    {
        for (auto hash : hash_values)
            [[maybe_unused]] auto & res = ibf.bulk_contains(hash);
    }

    state.counters["hashes/sec"] = hashes_per_second(std::ranges::size(hash_values));
}

BENCHMARK_TEMPLATE(emplace_benchmark,
                   seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>)->Apply(arguments);

BENCHMARK_TEMPLATE(bulk_contains_benchmark,
                   seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>)->Apply(arguments);
BENCHMARK_TEMPLATE(bulk_contains_benchmark,
                   seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed>)->Apply(arguments);

BENCHMARK_MAIN();
