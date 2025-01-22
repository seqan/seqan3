// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/views/rank_to.hpp>
#include <seqan3/search/fm_index/all.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/seqan2.hpp>
#include <seqan3/utility/range/to.hpp>
#include <seqan3/utility/type_traits/lazy_conditional.hpp>

#if SEQAN3_HAS_SEQAN2
#    include <seqan/index.h>
#endif

static constexpr int32_t max_length{50000};
static constexpr size_t seed{0x6'12'6f};

static void arguments(benchmark::internal::Benchmark * b)
{
#ifndef NDEBUG
    constexpr std::array<int32_t, 2> values{50, 5000};
#else
    constexpr std::array<int32_t, 3> values{50, 5000, 50000};
#endif // NDEBUG
    for (int32_t length : values)
    {
        if (length > max_length)
            throw std::logic_error{"Increase max_length to at least " + std::to_string(length)};

        b->Args({length, 5});
    }

    b->Args({500, 1000});
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
    std::string const char_rng{[]()
                               {
                                   std::vector<uint8_t> const ranks{
                                       seqan3::test::generate_numeric_sequence<uint8_t>(max_length, 0, 253, seed)};
                                   return ranks | seqan3::views::rank_to<char> | seqan3::ranges::to<std::string>();
                               }()};
};

sequence_store_seqan3 store{};

template <tag index_tag, typename rng_t>
void index_benchmark_seqan3(benchmark::State & state)
{
    using alphabet_t = seqan3::range_innermost_value_t<rng_t>;
    using inner_rng_t =
        std::conditional_t<seqan3::range_dimension_v<rng_t> == 1, rng_t, std::ranges::range_value_t<rng_t>>;

    rng_t sequence;
    inner_rng_t inner_sequence;
    if constexpr (std::same_as<alphabet_t, seqan3::dna4>)
        inner_sequence = store.dna4_rng | std::views::take(state.range(0)) | seqan3::ranges::to<inner_rng_t>();
    else if constexpr (std::same_as<alphabet_t, seqan3::aa27>)
        inner_sequence = store.aa27_rng | std::views::take(state.range(0)) | seqan3::ranges::to<inner_rng_t>();
    else
        inner_sequence = store.char_rng | std::views::take(state.range(0)) | seqan3::ranges::to<inner_rng_t>();

    if constexpr (seqan3::range_dimension_v<rng_t> == 1)
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
    seqan2::String<seqan2::Dna> const dna4_rng{
        seqan3::test::generate_sequence_seqan2<seqan2::Dna>(max_length, 0, seed)};
    seqan2::String<seqan2::AminoAcid> const aa27_rng{
        seqan3::test::generate_sequence_seqan2<seqan2::AminoAcid>(max_length, 0, seed)};
    seqan2::String<char> const char_rng{
        []()
        {
            std::vector<uint8_t> const ranks{
                seqan3::test::generate_numeric_sequence<uint8_t>(max_length, 0, 253, seed)};
            return ranks | seqan3::views::rank_to<char> | seqan3::ranges::to<std::string>();
        }()};
};

sequence_store_seqan2 store2{};

// Since the seqan2 alphabet has a value type, the dimension is actually one less than seqan3::range_dimension_v reports.
template <typename rng_t>
constexpr size_t seqan_dimension_v()
{
    if (seqan3::detail::is_type_specialisation_of_v<rng_t, seqan2::String>)
        return 1;
    if (seqan3::detail::is_type_specialisation_of_v<rng_t, seqan2::StringSet>)
        return 2;
    return 0;
}

template <tag index_tag, typename rng_t>
void index_benchmark_seqan2(benchmark::State & state)
{
    constexpr size_t dimension = seqan_dimension_v<rng_t>();
    static_assert(dimension != 0, "Use seqan2::String or seqan2::StringSet for SeqAn2 index benchmarks!");

    // Calling std::ranges::range_value_t twice on seqan2::String<char> is not valid.
    using alphabet_t = seqan3::detail::lazy_conditional_t<
        dimension == 1,
        std::ranges::range_value_t<rng_t>,
        seqan3::detail::lazy<std::ranges::range_value_t, std::ranges::range_value_t<rng_t>>>;

    using inner_rng_t = std::conditional_t<dimension == 1, rng_t, std::ranges::range_value_t<rng_t>>;
    using index_cfg = seqan2::FastFMIndexConfig<void, uint64_t>;
    using index_t =
        std::conditional_t<index_tag == tag::fm_index,
                           seqan2::Index<rng_t, seqan2::FMIndex<void, index_cfg>>,
                           seqan2::Index<rng_t, seqan2::BidirectionalIndex<seqan2::FMIndex<void, index_cfg>>>>;

    rng_t sequence;
    inner_rng_t inner_sequence;

    if constexpr (std::same_as<alphabet_t, seqan2::Dna>)
        inner_sequence = seqan2::prefix(store2.dna4_rng, state.range(0));
    else if constexpr (std::same_as<alphabet_t, seqan2::AminoAcid>)
        inner_sequence = seqan2::prefix(store2.aa27_rng, state.range(0));
    else
        inner_sequence = seqan2::prefix(store2.char_rng, state.range(0));

    if constexpr (dimension == 1)
    {
        sequence = std::move(inner_sequence);
    }
    else
    {
        for (int32_t i = 0; i < state.range(1); ++i)
            seqan2::appendValue(sequence, inner_sequence);
    }

    for (auto _ : state)
    {
        index_t index{sequence};
        seqan2::indexCreate(index, seqan2::FibreSALF());
    }
}
#endif // SEQAN3_HAS_SEQAN2

template <typename t>
using one_dimensional = std::conditional_t<std::same_as<std::string, t>, std::string, std::vector<t>>;
template <typename t>
using two_dimensional = std::vector<one_dimensional<t>>;

BENCHMARK_TEMPLATE(index_benchmark_seqan3, tag::fm_index, one_dimensional<seqan3::dna4>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan3, tag::fm_index, two_dimensional<seqan3::dna4>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan3, tag::fm_index, one_dimensional<seqan3::aa27>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan3, tag::fm_index, two_dimensional<seqan3::aa27>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan3, tag::fm_index, one_dimensional<std::string>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan3, tag::fm_index, two_dimensional<std::string>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan3, tag::bi_fm_index, one_dimensional<seqan3::dna4>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan3, tag::bi_fm_index, two_dimensional<seqan3::dna4>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan3, tag::bi_fm_index, one_dimensional<seqan3::aa27>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan3, tag::bi_fm_index, two_dimensional<seqan3::aa27>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan3, tag::bi_fm_index, one_dimensional<std::string>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan3, tag::bi_fm_index, two_dimensional<std::string>)->Apply(arguments);

#if SEQAN3_HAS_SEQAN2
template <typename t>
using one_dimensional2 = seqan2::String<t>;
template <typename t>
using two_dimensional2 = seqan2::StringSet<seqan2::String<t>>;

BENCHMARK_TEMPLATE(index_benchmark_seqan2, tag::fm_index, one_dimensional2<seqan2::Dna>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan2, tag::fm_index, two_dimensional2<seqan2::Dna>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan2, tag::fm_index, one_dimensional2<seqan2::AminoAcid>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan2, tag::fm_index, two_dimensional2<seqan2::AminoAcid>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan2, tag::fm_index, one_dimensional2<char>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan2, tag::fm_index, two_dimensional2<char>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan2, tag::bi_fm_index, one_dimensional2<seqan2::Dna>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan2, tag::bi_fm_index, two_dimensional2<seqan2::Dna>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan2, tag::bi_fm_index, one_dimensional2<seqan2::AminoAcid>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan2, tag::bi_fm_index, two_dimensional2<seqan2::AminoAcid>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan2, tag::bi_fm_index, one_dimensional2<char>)->Apply(arguments);
BENCHMARK_TEMPLATE(index_benchmark_seqan2, tag::bi_fm_index, two_dimensional2<char>)->Apply(arguments);
#endif // SEQAN3_HAS_SEQAN2

BENCHMARK_MAIN();
