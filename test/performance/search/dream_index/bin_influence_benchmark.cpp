// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <seqan3/range/views/repeat_n.hpp>
#include <seqan3/range/views/slice.hpp> 
#include <seqan3/range/views/zip.hpp>
#include <seqan3/search/dream_index/technical_binning_directory.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>

static constexpr size_t const genome_size{5000}; // 4'300'000'000
static constexpr size_t const read_size{100};
static constexpr size_t const read_count{1'000}; // 1'000'000
static constexpr size_t const ibf_size{8'388'608/*=1MiB*/}; // 68'719'476'736/*=8GiB*/

static std::vector<seqan3::dna4> const genome{seqan3::test::generate_sequence<seqan3::dna4>(genome_size, 0, 0)};
static std::vector<std::vector<seqan3::dna4>> const reads{[] (auto const & genome) { 
    std::vector<std::vector<seqan3::dna4>> result;
    auto const read_start_positions{seqan3::test::generate_numeric_sequence<size_t>(read_count, 0, genome_size - read_size + 1, 0)};
    for (size_t const & read_start : read_start_positions)
        result.emplace_back(genome | seqan3::views::slice(read_start, read_start + read_size) | seqan3::views::to<std::vector<seqan3::dna4>>);
    return result;
    } (genome)};

static void search_benchmark(benchmark::State & state)
{
    size_t const bin_count = static_cast<size_t>(state.range(0));
    size_t const hash_num{2u};
    size_t const bin_size{ibf_size / bin_count};
    size_t const chunk_size{(genome_size + bin_count - 1) / bin_count};

    seqan3::ibf_config const cfg{seqan3::bin_count{bin_count},
                                 seqan3::bin_size{bin_size},
                                 seqan3::hash_function_count{hash_num}};

    seqan3::technical_binning_directory tbd{genome | ranges::views::chunk(chunk_size),
                                            seqan3::views::kmer_hash(seqan3::ungapped{19u}),
                                            cfg};

    auto agent = tbd.counting_agent();
    for (auto _ : state)
        for (auto && query : reads)
            benchmark::DoNotOptimize(agent.count_query(query));
}

BENCHMARK(search_benchmark)->RangeMultiplier(2)->Range(64, 128); // 65536

BENCHMARK_MAIN();
