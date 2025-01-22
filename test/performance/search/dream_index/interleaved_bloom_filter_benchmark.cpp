// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h>

#include <seqan3/search/dream_index/interleaved_bloom_filter.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/utility/range/to.hpp>
#include <seqan3/utility/views/zip.hpp>

inline benchmark::Counter hashes_per_second(size_t const count)
{
    return benchmark::Counter(count, benchmark::Counter::kIsIterationInvariantRate, benchmark::Counter::OneK::kIs1000);
}

static void arguments(benchmark::internal::Benchmark * b)
{
    // Bins must be powers of two
    for (int32_t bins : {64, 8192})
    {
        // Size of the IBF will be 2^bits bits
        for (int32_t bits = 15; bits <= 20; bits += 5)
        {
            // The bits per bin must fit in an int32_t
            if (bits - std::countr_zero(static_cast<uint32_t>(bins)) < 32)
            {
                for (int32_t hash_num = 2; hash_num < 3; ++hash_num)
                {
                    b->Args({bins, (1LL << bits) / bins, hash_num, 1000});
                }
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
    auto && [bin_indices, hash_values, ibf] =
        set_up<ibf_type>(state.range(0), state.range(1), state.range(2), state.range(3));

    for (auto _ : state)
    {
        for (auto [hash, bin] : seqan3::views::zip(hash_values, bin_indices))
            ibf.emplace(hash, seqan3::bin_index{bin});
    }

    state.counters["hashes/sec"] = hashes_per_second(std::ranges::size(hash_values));
}

template <typename ibf_type>
void clear_benchmark(::benchmark::State & state)
{
    auto && [bin_indices, hash_values, ibf] =
        set_up<ibf_type>(state.range(0), state.range(1), state.range(2), state.range(3));
    (void)bin_indices;
    (void)hash_values;

    std::vector<seqan3::bin_index> bin_range = std::views::iota(0u, static_cast<size_t>(state.range(0)))
                                             | std::views::transform(
                                                   [](size_t i)
                                                   {
                                                       return seqan3::bin_index{i};
                                                   })
                                             | seqan3::ranges::to<std::vector>();

    for (auto _ : state)
    {
        for (auto bin : bin_range)
            ibf.clear(bin);
    }

    state.counters["bins/sec"] = hashes_per_second(std::ranges::size(bin_range));
}

template <typename ibf_type>
void clear_range_benchmark(::benchmark::State & state)
{
    auto && [bin_indices, hash_values, ibf] =
        set_up<ibf_type>(state.range(0), state.range(1), state.range(2), state.range(3));
    (void)bin_indices;
    (void)hash_values;

    std::vector<seqan3::bin_index> bin_range = std::views::iota(0u, static_cast<size_t>(state.range(0)))
                                             | std::views::transform(
                                                   [](size_t i)
                                                   {
                                                       return seqan3::bin_index{i};
                                                   })
                                             | seqan3::ranges::to<std::vector>();

    for (auto _ : state)
    {
        ibf.clear(bin_range);
    }

    state.counters["bins/sec"] = hashes_per_second(std::ranges::size(bin_range));
}

template <typename ibf_type>
void bulk_contains_benchmark(::benchmark::State & state)
{
    auto && [bin_indices, hash_values, ibf] =
        set_up<ibf_type>(state.range(0), state.range(1), state.range(2), state.range(3));
    (void)bin_indices;

    auto agent = ibf.membership_agent();
    for (auto _ : state)
    {
        for (auto hash : hash_values) [[maybe_unused]]
            auto & res = agent.bulk_contains(hash);
    }

    state.counters["hashes/sec"] = hashes_per_second(std::ranges::size(hash_values));
}

template <typename ibf_type>
void bulk_count_benchmark(::benchmark::State & state)
{
    auto && [bin_indices, hash_values, ibf] =
        set_up<ibf_type>(state.range(0), state.range(1), state.range(2), state.range(3));
    (void)bin_indices;

    auto agent = ibf.counting_agent();
    for (auto _ : state)
    {
        [[maybe_unused]] auto & res = agent.bulk_count(hash_values);
    }

    state.counters["hashes/sec"] = hashes_per_second(std::ranges::size(hash_values));
}

BENCHMARK_TEMPLATE(emplace_benchmark, seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>)
    ->Apply(arguments);
BENCHMARK_TEMPLATE(clear_benchmark, seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>)
    ->Apply(arguments);
BENCHMARK_TEMPLATE(clear_range_benchmark, seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>)
    ->Apply(arguments);

BENCHMARK_TEMPLATE(bulk_contains_benchmark, seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>)
    ->Apply(arguments);
BENCHMARK_TEMPLATE(bulk_contains_benchmark, seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed>)
    ->Apply(arguments);

BENCHMARK_TEMPLATE(bulk_count_benchmark, seqan3::interleaved_bloom_filter<seqan3::data_layout::uncompressed>)
    ->Apply(arguments);
BENCHMARK_TEMPLATE(bulk_count_benchmark, seqan3::interleaved_bloom_filter<seqan3::data_layout::compressed>)
    ->Apply(arguments);

BENCHMARK_MAIN();
