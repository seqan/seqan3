// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/type_traits/lazy.hpp>
#include <seqan3/range/views/persist.hpp>
#include <seqan3/range/views/rank_to.hpp>
#include <seqan3/search/fm_index/all.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/seqan2.hpp>

#if SEQAN3_HAS_SEQAN2
#include <seqan/index.h>
#endif

static constexpr int32_t max_length{50'000};
static constexpr size_t seed{0x6126f};

static void arguments(benchmark::internal::Benchmark * b)
{
    for (int32_t length : {50, 5000, 50'000})
    {
        if (length > max_length)
            throw std::logic_error{"Increase max_length to at least " + std::to_string(length)};

        b->Args({length, 5});
    }

    b->Args({500, 1'000});
}

enum class tag
{
    fm_index,
    bi_fm_index
};

struct sequence_store_seqan3
{
    std::vector<seqan3::dna4> const dna4_rng{seqan3::test::generate_sequence<seqan3::dna4>(max_length, 0, seed)};
    std::vector<seqan3::aa27> const aa27_rng{seqan3::test::generate_sequence<seqan3::aa27>(max_length, 0, seed)};
    std::string const char_rng{seqan3::test::generate_numeric_sequence<uint8_t>(max_length, 0, 253, seed)
                               | seqan3::views::persist
                               | seqan3::views::rank_to<char>
                               | seqan3::views::to<std::string>};
};

sequence_store_seqan3 store{};

template <tag index_tag, typename rng_t>
void index_benchmark_seqan3(benchmark::State & state)
{
    using alphabet_t = seqan3::innermost_value_type_t<rng_t>;
    using inner_rng_t = std::conditional_t<seqan3::dimension_v<rng_t> == 1, rng_t, std::ranges::range_value_t<rng_t>>;

    rng_t sequence;
    inner_rng_t inner_sequence;
    if constexpr (std::same_as<alphabet_t, seqan3::dna4>)
        inner_sequence = store.dna4_rng | seqan3::views::take(state.range(0)) | seqan3::views::to<inner_rng_t>;
    else if constexpr (std::same_as<alphabet_t, seqan3::aa27>)
        inner_sequence = store.aa27_rng | seqan3::views::take(state.range(0)) | seqan3::views::to<inner_rng_t>;
    else
        inner_sequence = store.char_rng | seqan3::views::take(state.range(0)) | seqan3::views::to<inner_rng_t>;

    if constexpr (seqan3::dimension_v<rng_t> == 1)
    {
        sequence = std::move(inner_sequence);
    }
    else
    {
        for (int32_t i = 0; i < state.range(1); ++i)
            sequence.push_back(inner_sequence);
    }

    for (auto _ : state)
    {
        if constexpr (index_tag == tag::fm_index)
            seqan3::fm_index index{sequence};
        else
            seqan3::bi_fm_index index{sequence};
    }
}

#if SEQAN3_HAS_SEQAN2
struct sequence_store_seqan2
{
    seqan::String<seqan::Dna> const dna4_rng{seqan3::test::generate_sequence_seqan2<seqan::Dna>(max_length, 0, seed)};
    seqan::String<seqan::AminoAcid> const aa27_rng{seqan3::test::generate_sequence_seqan2<seqan::AminoAcid>(max_length,
                                                                                                            0,
                                                                                                            seed)};
    seqan::String<char> const char_rng{seqan3::test::generate_numeric_sequence<uint8_t>(max_length, 0, 253, seed)
                                       | seqan3::views::persist
                                       | seqan3::views::rank_to<char>
                                       | seqan3::views::to<std::string>};
};

sequence_store_seqan2 store2{};

// Since the seqan2 alphabet has a value type, the dimension is actually one less than seqan3::dimension_v reports.
template <typename rng_t>
constexpr size_t seqan_dimension_v()
{
    if (seqan3::detail::is_type_specialisation_of_v<rng_t, seqan::String>)
        return 1;
    if (seqan3::detail::is_type_specialisation_of_v<rng_t, seqan::StringSet>)
        return 2;
   return 0;
}

template <tag index_tag, typename rng_t>
void index_benchmark_seqan2(benchmark::State & state)
{
    constexpr size_t dimension = seqan_dimension_v<rng_t>();
    static_assert(dimension != 0, "Use seqan::String or seqan::StringSet for SeqAn2 index benchmarks!");

    // Calling std::ranges::range_value_t twice on seqan::String<char> is not valid.
    using alphabet_t = seqan3::detail::lazy_conditional_t<dimension == 1,
                                                         std::ranges::range_value_t<rng_t>,
                                                         seqan3::detail::lazy<std::ranges::range_value_t,
                                                                              std::ranges::range_value_t<rng_t>>>;

    using inner_rng_t = std::conditional_t<dimension == 1, rng_t, std::ranges::range_value_t<rng_t>>;
    using index_cfg = seqan::FastFMIndexConfig<void, uint64_t>;
    using index_t = std::conditional_t<index_tag == tag::fm_index,
                                       seqan::Index<rng_t, seqan::FMIndex<void, index_cfg>>,
                                       seqan::Index<rng_t, seqan::BidirectionalIndex<seqan::FMIndex<void, index_cfg>>>>;

    rng_t sequence;
    inner_rng_t inner_sequence;

    if constexpr (std::same_as<alphabet_t, seqan::Dna>)
        inner_sequence = seqan::prefix(store2.dna4_rng, state.range(0));
    else if constexpr (std::same_as<alphabet_t, seqan::AminoAcid>)
        inner_sequence = seqan::prefix(store2.aa27_rng, state.range(0));
    else
        inner_sequence = seqan::prefix(store2.char_rng, state.range(0));

    if constexpr (dimension == 1)
    {
        sequence = std::move(inner_sequence);
    }
    else
    {
        for (int32_t i = 0; i < state.range(1); ++i)
            seqan::appendValue(sequence, inner_sequence);
    }

    for (auto _ : state)
    {
        index_t index{sequence};
        seqan::indexCreate(index, seqan::FibreSALF());
    }
}
#endif // SEQAN3_HAS_SEQAN2

template <typename t>
using one_dimensional = std::conditional_t<std::same_as<std::string, t>, std::string, std::vector<t>>;
template <typename t>
using two_dimensional = std::vector<one_dimensional<t>>;

BENCHMARK_TEMPLATE(index_benchmark_seqan3, tag::fm_index,    one_dimensional<seqan3::dna4>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan3, tag::fm_index,    two_dimensional<seqan3::dna4>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan3, tag::fm_index,    one_dimensional<seqan3::aa27>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan3, tag::fm_index,    two_dimensional<seqan3::aa27>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan3, tag::fm_index,    one_dimensional<std::string> )->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan3, tag::fm_index,    two_dimensional<std::string> )->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan3, tag::bi_fm_index, one_dimensional<seqan3::dna4>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan3, tag::bi_fm_index, two_dimensional<seqan3::dna4>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan3, tag::bi_fm_index, one_dimensional<seqan3::aa27>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan3, tag::bi_fm_index, two_dimensional<seqan3::aa27>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan3, tag::bi_fm_index, one_dimensional<std::string> )->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan3, tag::bi_fm_index, two_dimensional<std::string> )->Apply(arguments);

#if SEQAN3_HAS_SEQAN2
template <typename t>
using one_dimensional2 = seqan::String<t>;
template <typename t>
using two_dimensional2 = seqan::StringSet<seqan::String<t>>;

BENCHMARK_TEMPLATE(index_benchmark_seqan2, tag::fm_index,    one_dimensional2<seqan::Dna>      )->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan2, tag::fm_index,    two_dimensional2<seqan::Dna>      )->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan2, tag::fm_index,    one_dimensional2<seqan::AminoAcid>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan2, tag::fm_index,    two_dimensional2<seqan::AminoAcid>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan2, tag::fm_index,    one_dimensional2<char>            )->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan2, tag::fm_index,    two_dimensional2<char>            )->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan2, tag::bi_fm_index, one_dimensional2<seqan::Dna>      )->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan2, tag::bi_fm_index, two_dimensional2<seqan::Dna>      )->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan2, tag::bi_fm_index, one_dimensional2<seqan::AminoAcid>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan2, tag::bi_fm_index, two_dimensional2<seqan::AminoAcid>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan2, tag::bi_fm_index, one_dimensional2<char>            )->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan2, tag::bi_fm_index, two_dimensional2<char>            )->Apply(arguments);
#endif // SEQAN3_HAS_SEQAN2

BENCHMARK_MAIN();
