// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/test/performance/naive_minimiser_hash.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/performance/units.hpp>

#ifdef SEQAN3_HAS_SEQAN2
#include <seqan/index.h>
#include <seqan3/test/performance/seqan2_minimiser.h>
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

static void arguments(benchmark::internal::Benchmark* b)
{
    for (int32_t sequence_length : {50'000, /*1'000'000*/})
    {
        for (int32_t k : {8, /*16, 24,*/ 30})
        {
            for (int32_t w : {k + 5, k + 20})
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
    naive,
    seqan2_ungapped,
    seqan2_gapped
};

#ifdef SEQAN3_HAS_SEQAN2
inline auto make_gapped_shape_seqan2(size_t const k)
{
    seqan::String<char> bitmap;

    for (size_t i{0}; i < k - 1; ++i)
        seqan::append(bitmap, seqan::CharString(std::to_string((i + 1) % 2)));

    seqan::append(bitmap, seqan::CharString("1"));

    return seqan::Shape<seqan::Dna, seqan::GenericShape>(bitmap);
}
#endif // SEQAN3_HAS_SEQAN2

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
        else if constexpr (tag == method_tag::seqan3_ungapped)
        {
            for (auto h : seq | seqan3::views::minimiser_hash(seqan3::ungapped{static_cast<uint8_t>(k)}, seqan3::window_size{w}))
                benchmark::DoNotOptimize(sum += h);
        }
        else if constexpr (tag == method_tag::seqan3_gapped)
        {
            for (auto h : seq | seqan3::views::minimiser_hash(make_gapped_shape(k), seqan3::window_size{w}))
                benchmark::DoNotOptimize(sum += h);
        }
        #ifdef SEQAN3_HAS_SEQAN2
        else
        {
            auto seqan2_seq = seqan3::test::generate_sequence_seqan2<seqan::Dna>(sequence_length, 0, 0);
            using shape_t = std::conditional_t<tag == method_tag::seqan2_ungapped,
                                               seqan::Shape<seqan::Dna, seqan::SimpleShape>,
                                               seqan::Shape<seqan::Dna, seqan::GenericShape>>;

            shape_t shape;
            if constexpr (tag == method_tag::seqan2_ungapped)
               seqan::resize(shape, k);
            else
                shape = make_gapped_shape_seqan2(k);

            minimiser<decltype(shape)> seqan_minimiser(window{w}, kmer{k}, shape);
            seqan_minimiser.compute(seqan2_seq);

            for (auto h : seqan_minimiser.minimiser_hash)
                benchmark::DoNotOptimize(sum += h);
        }
        #endif // SEQAN3_HAS_SEQAN2
    }

    state.counters["Throughput[bp/s]"] = bp_per_second(sequence_length - k + 1);
}

template <method_tag tag>
void compute_minimisers_on_poly_A_sequence(benchmark::State & state)
{
    auto sequence_length = state.range(0);
    size_t k = static_cast<size_t>(state.range(1));
    uint32_t w = static_cast<size_t>(state.range(2));
    assert(sequence_length > 0);
    assert(k > 0);
    assert(w > k);
    auto seq = std::vector<seqan3::dna4>(sequence_length);

    size_t sum{0};

    for (auto _ : state)
    {
        if constexpr (tag == method_tag::seqan3_ungapped)
        {
            for (auto h : seq | seqan3::views::minimiser_hash(seqan3::ungapped{static_cast<uint8_t>(k)}, seqan3::window_size{w}))
                benchmark::DoNotOptimize(sum += h);
        }
        else if constexpr (tag == method_tag::seqan3_gapped)
        {
            for (auto h : seq | seqan3::views::minimiser_hash(make_gapped_shape(k), seqan3::window_size{w}))
                benchmark::DoNotOptimize(sum += h);
        }
    }

    state.counters["Throughput[bp/s]"] = bp_per_second(sequence_length - k + 1);
}


#ifdef SEQAN3_HAS_SEQAN2
BENCHMARK_TEMPLATE(compute_minimisers, method_tag::seqan2_ungapped)->Apply(arguments);
BENCHMARK_TEMPLATE(compute_minimisers, method_tag::seqan2_gapped)->Apply(arguments);
#endif // SEQAN3_HAS_SEQAN2

BENCHMARK_TEMPLATE(compute_minimisers, method_tag::naive)->Apply(arguments);
BENCHMARK_TEMPLATE(compute_minimisers, method_tag::seqan3_ungapped)->Apply(arguments);
BENCHMARK_TEMPLATE(compute_minimisers, method_tag::seqan3_gapped)->Apply(arguments);

BENCHMARK_TEMPLATE(compute_minimisers_on_poly_A_sequence, method_tag::seqan3_ungapped)->Apply(arguments);
BENCHMARK_TEMPLATE(compute_minimisers_on_poly_A_sequence, method_tag::seqan3_gapped)->Apply(arguments);

BENCHMARK_MAIN();
