// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/search/views/minimiser_hash.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/performance/units.hpp>
#include <seqan3/utility/views/zip.hpp>

#ifdef SEQAN3_HAS_SEQAN2
#    include <seqan3/test/performance/seqan2_minimiser.h>

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
    for (int32_t sequence_length : {50000, /*1'000'000*/})
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
    seqan2::String<char> bitmap;

    for (size_t i{0}; i < k - 1; ++i)
        seqan2::append(bitmap, seqan2::CharString(std::to_string((i + 1) % 2)));

    seqan2::append(bitmap, seqan2::CharString("1"));

    return seqan2::Shape<seqan2::Dna, seqan2::GenericShape>(bitmap);
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
            seqan3::shape shape = seqan3::ungapped{static_cast<uint8_t>(k)};
            uint64_t const seed = 0x8F'3F'73'B5'CF'1C'9A'DE;

            // Use random seed to randomise order on forward strand.
            auto forward = seq | seqan3::views::kmer_hash(shape)
                         | std::views::transform(
                               [](uint64_t const i)
                               {
                                   return i ^ seed;
                               });
            // Create reverse complement strand and use random seed to randomise order on reverse complement strand.
            auto reverse = seq | seqan3::views::complement // Create complement.
                         | std::views::reverse             // Reverse order.
                         | seqan3::views::kmer_hash(shape) // Get hash values.
                         | std::views::transform(
                               [](uint64_t const i)
                               {
                                   return i ^ seed;
                               })               // Randomise.
                         | std::views::reverse; // Reverse again, so that the first hash value
                                                // is the reverse complement of the first
                                                // hash value in the forward strand.
            // Get minimum between forward and reverse strand for each value.
            auto both = seqan3::views::zip(forward, reverse)
                      | std::views::transform(
                            [](auto && fwd_rev_hash_pair)
                            {
                                return std::min(std::get<0>(fwd_rev_hash_pair), std::get<1>(fwd_rev_hash_pair));
                            });

            // Setup to slide over a range with `w - shape.size()+1` window size.
            // Initialize a sub range of size ` w - shape.size()`
            auto subrange_begin = begin(both);
            auto subrange_end = std::ranges::next(subrange_begin, w - shape.size(), end(both));

            // Slide over the range
            while (subrange_end != end(both))
            {
                ++subrange_end; // Extends the subrange to `w - shape.size()+1`
                auto h = *std::ranges::min_element(subrange_begin, subrange_end);

                ++subrange_begin; // Move the beginning one forward
                benchmark::DoNotOptimize(sum += h);
            }
        }
        else if constexpr (tag == method_tag::seqan3_ungapped)
        {
            for (auto h :
                 seq | seqan3::views::minimiser_hash(seqan3::ungapped{static_cast<uint8_t>(k)}, seqan3::window_size{w}))
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
            auto seqan2_seq = seqan3::test::generate_sequence_seqan2<seqan2::Dna>(sequence_length, 0, 0);
            using shape_t = std::conditional_t<tag == method_tag::seqan2_ungapped,
                                               seqan2::Shape<seqan2::Dna, seqan2::SimpleShape>,
                                               seqan2::Shape<seqan2::Dna, seqan2::GenericShape>>;

            shape_t shape;
            if constexpr (tag == method_tag::seqan2_ungapped)
                seqan2::resize(shape, k);
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
            for (auto h :
                 seq | seqan3::views::minimiser_hash(seqan3::ungapped{static_cast<uint8_t>(k)}, seqan3::window_size{w}))
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
