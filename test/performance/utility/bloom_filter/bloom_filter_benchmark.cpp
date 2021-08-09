// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <seqan3/range/views/to.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/utility/bloom_filter/bloom_filter.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>

inline benchmark::Counter hashes_per_second(size_t const count)
{
    return benchmark::Counter(count,
                              benchmark::Counter::kIsIterationInvariantRate,
                              benchmark::Counter::OneK::kIs1000);
}

static void arguments(benchmark::internal::Benchmark* b)
{
    // Size of the IBF will be 2^bits bits
    for (int32_t bits = 15; bits <= 20; bits += 5)
    {
        // The bits must fit in an int32_t
        if (bits < 32)
        {
            for (int32_t hash_num = 2; hash_num < 3; ++hash_num)
            {
                b->Args({(1LL << bits), hash_num, 1'000});
            }
        }
    }
}

template <typename bf_type>
auto set_up(size_t bits, size_t hash_num, size_t sequence_length)
{
    auto hash_values = seqan3::test::generate_numeric_sequence<size_t>(sequence_length);
    seqan3::bloom_filter tmp_bf(seqan3::bin_size{bits},
                                seqan3::hash_function_count{hash_num});

    bf_type bf{std::move(tmp_bf)};

    return std::make_tuple(hash_values, bf);
}

template <typename ibf_type>
void emplace_benchmark(::benchmark::State & state)
{
    auto && [ hash_values, bf ] = set_up<ibf_type>(state.range(0),
                                                    state.range(1),
                                                    state.range(2));

    for (auto _ : state)
    {
        for (auto hash : hash_values)
            bf.emplace(hash);
    }

    state.counters["hashes/sec"] = hashes_per_second(std::ranges::size(hash_values));
}

template <typename bf_type>
void clear_benchmark(::benchmark::State & state)
{
    auto && [ hash_values, bf ] = set_up<bf_type>(state.range(0),
                                                  state.range(1),
                                                  state.range(2));
    (void) hash_values;

    for (auto _ : state)
    {
        bf.clear();
    }
}

template <typename bf_type>
void contains_benchmark(::benchmark::State & state)
{
    auto && [ hash_values, bf ] = set_up<bf_type>(state.range(0),
                                                  state.range(1),
                                                  state.range(2));

    for (auto _ : state)
    {
        for (auto hash : hash_values)
            benchmark::DoNotOptimize(bf.contains(hash));
    }

    state.counters["hashes/sec"] = hashes_per_second(std::ranges::size(hash_values));
}

template <typename bf_type>
void count_benchmark(::benchmark::State & state)
{
    auto && [ hash_values, bf ] = set_up<bf_type>(state.range(0),
                                                  state.range(1),
                                                  state.range(2));

    for (auto _ : state)
    {
        benchmark::DoNotOptimize(bf.count(hash_values));
    }

    state.counters["hashes/sec"] = hashes_per_second(std::ranges::size(hash_values));
}

BENCHMARK_TEMPLATE(emplace_benchmark,
                   seqan3::bloom_filter<seqan3::data_layout::uncompressed>)->Apply(arguments);
BENCHMARK_TEMPLATE(clear_benchmark,
                   seqan3::bloom_filter<seqan3::data_layout::uncompressed>)->Apply(arguments);

BENCHMARK_TEMPLATE(contains_benchmark,
                   seqan3::bloom_filter<seqan3::data_layout::uncompressed>)->Apply(arguments);
BENCHMARK_TEMPLATE(contains_benchmark,
                   seqan3::bloom_filter<seqan3::data_layout::compressed>)->Apply(arguments);

BENCHMARK_TEMPLATE(count_benchmark,
                   seqan3::bloom_filter<seqan3::data_layout::uncompressed>)->Apply(arguments);
BENCHMARK_TEMPLATE(count_benchmark,
                   seqan3::bloom_filter<seqan3::data_layout::compressed>)->Apply(arguments);

BENCHMARK_MAIN();
