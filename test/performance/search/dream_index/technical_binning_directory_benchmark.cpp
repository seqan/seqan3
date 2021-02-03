// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <seqan3/range/views/repeat_n.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/search/dream_index/technical_binning_directory.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>

using seqan3::operator""_dna4;

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

template <typename tbd_type>
auto set_up(size_t bins, size_t bits, size_t hash_num, size_t sequence_length)
{
    seqan3::ibf_config cfg{seqan3::bin_count{bins},
                           seqan3::bin_size{bits},
                           seqan3::hash_function_count{hash_num}};

    auto bin_indices = seqan3::test::generate_numeric_sequence<size_t>(sequence_length, 0u, bins - 1);
    auto hash_values = seqan3::test::generate_numeric_sequence<size_t>(sequence_length);
    auto v = seqan3::views::kmer_hash(seqan3::ungapped{20u});
    seqan3::technical_binning_directory tmp_tbd{seqan3::views::repeat_n(std::vector<seqan3::dna4>{}, bins),
                                                std::move(v),
                                                cfg};

    tbd_type tbd{std::move(tmp_tbd)};

    return std::make_tuple(bin_indices, hash_values, tbd);
}

template <typename tbd_type>
void emplace_benchmark(::benchmark::State & state)
{
    auto && [ bin_indices, hash_values, tbd ] = set_up<tbd_type>(state.range(0),
                                                                 state.range(1),
                                                                 state.range(2),
                                                                 state.range(3));

    for (auto _ : state)
    {
        for (auto [hash, bin] : seqan3::views::zip(hash_values, bin_indices))
            tbd.emplace(hash, seqan3::bin_index{bin});
    }

    state.counters["hashes/sec"] = hashes_per_second(std::ranges::size(hash_values));
}

template <typename tbd_type>
void bulk_contains_benchmark(::benchmark::State & state)
{
    auto && [ bin_indices, hash_values, tbd ] = set_up<tbd_type>(state.range(0),
                                                                 state.range(1),
                                                                 state.range(2),
                                                                 state.range(3));
    (void) bin_indices;

    auto agent = tbd.membership_agent();
    for (auto _ : state)
    {
        for (auto hash : hash_values)
            [[maybe_unused]] auto & res = agent.bulk_contains(hash);
    }

    state.counters["hashes/sec"] = hashes_per_second(std::ranges::size(hash_values));
}

template <typename tbd_type>
void count_benchmark(::benchmark::State & state)
{
    auto && [ bin_indices, hash_values, tbd ] = set_up<tbd_type>(state.range(0),
                                                                 state.range(1),
                                                                 state.range(2),
                                                                 state.range(3));
    (void) bin_indices;
    (void) hash_values;

    auto agent = tbd.counting_agent();
    auto const query = std::vector<seqan3::dna4>{100, 'A'_dna4};
    for (auto _ : state)
    {
            [[maybe_unused]] auto & res = agent.count_query(query);
    }

    state.counters["hashes/sec"] = hashes_per_second(81); // 100-20+1 per call to `count_query`
}

BENCHMARK_TEMPLATE(emplace_benchmark,
                   seqan3::technical_binning_directory<seqan3::data_layout::uncompressed>)->Apply(arguments);

BENCHMARK_TEMPLATE(bulk_contains_benchmark,
                   seqan3::technical_binning_directory<seqan3::data_layout::uncompressed>)->Apply(arguments);
BENCHMARK_TEMPLATE(bulk_contains_benchmark,
                   seqan3::technical_binning_directory<seqan3::data_layout::compressed>)->Apply(arguments);

BENCHMARK_TEMPLATE(count_benchmark,
                   seqan3::technical_binning_directory<seqan3::data_layout::uncompressed>)->Apply(arguments);
BENCHMARK_TEMPLATE(count_benchmark,
                   seqan3::technical_binning_directory<seqan3::data_layout::compressed>)->Apply(arguments);

BENCHMARK_MAIN();
