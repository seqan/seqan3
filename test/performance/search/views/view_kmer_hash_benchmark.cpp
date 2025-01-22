// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/performance/units.hpp>

#ifdef SEQAN3_HAS_SEQAN2
#    include <seqan/index.h>
#endif // SEQAN3_HAS_SEQAN2

inline benchmark::Counter bp_per_second(size_t const basepairs)
{
    return benchmark::Counter(basepairs,
                              benchmark::Counter::kIsIterationInvariantRate,
                              benchmark::Counter::OneK::kIs1000);
}

inline seqan3::shape make_gapped_shape(size_t const k)
{
    seqan3::shape shape{};

    for (size_t i{0}; i < k - 1; ++i)
        shape.push_back((i + 1) % 2);

    shape.push_back(1u);
    shape.push_back(0u);
    return shape;
}

static void arguments(benchmark::internal::Benchmark * b)
{
    for (int32_t sequence_length : {1000, 50000, /*1'000'000*/})
    {
        for (int32_t k : {8, /*16, 24,*/ 30})
        {
            b->Args({sequence_length, k});
        }
    }
}

static void seqan_kmer_hash_ungapped(benchmark::State & state)
{
    auto sequence_length = state.range(0);
    assert(sequence_length > 0);
    size_t k = static_cast<size_t>(state.range(1));
    assert(k > 0);
    auto seq = seqan3::test::generate_sequence<seqan3::dna4>(sequence_length, 0, 0);

    size_t sum{0};

    for (auto _ : state)
    {
        for (auto h : seq | seqan3::views::kmer_hash(seqan3::ungapped{static_cast<uint8_t>(k)}))
            benchmark::DoNotOptimize(sum += h);
    }

    // prevent complete optimisation
    [[maybe_unused]] volatile auto fin = sum;

    state.counters["Throughput[bp/s]"] = bp_per_second(sequence_length - k + 1);
}

static void seqan_kmer_hash_gapped(benchmark::State & state)
{
    auto sequence_length = state.range(0);
    assert(sequence_length > 0);
    size_t k = static_cast<size_t>(state.range(1));
    assert(k > 0);
    auto seq = seqan3::test::generate_sequence<seqan3::dna4>(sequence_length, 0, 0);

    size_t sum{0};

    for (auto _ : state)
    {
        for (auto h : seq | seqan3::views::kmer_hash(make_gapped_shape(k)))
            benchmark::DoNotOptimize(sum += h);
    }

    // prevent complete optimisation
    [[maybe_unused]] volatile auto fin = sum;

    state.counters["Throughput[bp/s]"] = bp_per_second(sequence_length - k + 1);
}

static void naive_kmer_hash(benchmark::State & state)
{
    auto sequence_length = state.range(0);
    assert(sequence_length > 0);
    size_t k = static_cast<size_t>(state.range(1));
    assert(k > 0);
    auto seq = seqan3::test::generate_sequence<seqan3::dna4>(sequence_length, 0, 0);

    size_t sum{0};

    for (auto _ : state)
    {
        // creates a initial subrange covering exactly 'k-1' characters
        auto subrange_begin = begin(seq);
        auto subrange_end = std::ranges::next(subrange_begin, k - 1, end(seq));

        // slide over the range
        while (subrange_end != end(seq))
        {
            ++subrange_end; // extend to 'k' characters
            auto r = std::ranges::subrange{subrange_begin, subrange_end};
            auto h = std::hash<decltype(r)>{}(r);
            ++subrange_begin; // move front forward, back to 'k-1' characters range

            benchmark::DoNotOptimize(sum += h);
        }
    }

    // prevent complete optimisation
    [[maybe_unused]] volatile auto fin = sum;

    state.counters["Throughput[bp/s]"] = bp_per_second(sequence_length - k + 1);
}

#ifdef SEQAN3_HAS_SEQAN2
inline auto make_gapped_shape_seqan2(size_t const k)
{
    seqan2::String<char> bitmap;

    for (size_t i{0}; i < k - 1; ++i)
        seqan2::append(bitmap, seqan2::CharString(std::to_string((i + 1) % 2)));

    seqan2::append(bitmap, seqan2::CharString("1"));

    seqan2::Shape<seqan2::Dna, seqan2::GenericShape> s(bitmap);
    return s;
}

static void seqan2_kmer_hash_ungapped(benchmark::State & state)
{
    auto sequence_length = state.range(0);
    assert(sequence_length > 0);
    size_t k = static_cast<size_t>(state.range(1));
    assert(k > 0);
    auto seq = seqan3::test::generate_sequence_seqan2<seqan2::Dna>(sequence_length, 0, 0);
    seqan2::Shape<seqan2::Dna, seqan2::SimpleShape> s;
    seqan2::resize(s, k);

    size_t sum{0};

    for (auto _ : state)
    {
        auto it = seqan2::begin(seq);
        seqan2::hashInit(s, it);
        for (size_t i{0}; i < seqan2::length(seq) - k + 1; ++i, ++it)
        {
            benchmark::DoNotOptimize(sum += seqan2::hashNext(s, it));
        }
    }

    // prevent complete optimisation
    [[maybe_unused]] volatile auto fin = sum;

    state.counters["Throughput[bp/s]"] = bp_per_second(sequence_length - k + 1);
}

static void seqan2_kmer_hash_gapped(benchmark::State & state)
{
    auto sequence_length = state.range(0);
    assert(sequence_length > 0);
    size_t k = static_cast<size_t>(state.range(1));
    assert(k > 0);
    auto seq = seqan3::test::generate_sequence_seqan2<seqan2::Dna>(sequence_length, 0, 0);
    seqan2::Shape<seqan2::Dna, seqan2::GenericShape> s = make_gapped_shape_seqan2(k);

    size_t sum{0};

    for (auto _ : state)
    {
        auto it = seqan2::begin(seq);
        seqan2::hashInit(s, it);
        for (size_t i{0}; i < seqan2::length(seq) - k + 1; ++i, ++it)
        {
            benchmark::DoNotOptimize(sum += seqan2::hashNext(s, it));
        }
    }

    // prevent complete optimisation
    [[maybe_unused]] volatile auto fin = sum;

    state.counters["Throughput[bp/s]"] = bp_per_second(sequence_length - k + 1);
}

BENCHMARK(seqan2_kmer_hash_ungapped)->Apply(arguments);
BENCHMARK(seqan2_kmer_hash_gapped)->Apply(arguments);
#endif // SEQAN3_HAS_SEQAN2

BENCHMARK(seqan_kmer_hash_ungapped)->Apply(arguments);
BENCHMARK(seqan_kmer_hash_gapped)->Apply(arguments);
BENCHMARK(naive_kmer_hash)->Apply(arguments);

BENCHMARK_MAIN();
