// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <benchmark/benchmark.h>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/test/performance/units.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/seqan2.hpp>

#ifdef SEQAN3_HAS_SEQAN2
#include <seqan/align.h>
#include <seqan/align_parallel.h>
#endif

// ----------------------------------------------------------------------------
// Sequence pair generators for benchmarks
// ----------------------------------------------------------------------------

template <typename alphabet_t>
struct seqan3_sequence_pair_generator
{
    using sequence_t = std::vector<alphabet_t>;

    auto operator()(benchmark::State &) const
    {
        auto sequence1 = seqan3::test::generate_sequence<alphabet_t>(sequence_length, 0, 0);
        auto sequence2 = seqan3::test::generate_sequence<alphabet_t>(sequence_length, 0, 1);
        return std::tuple{std::move(sequence1), std::move(sequence2)};
    };

    size_t sequence_length;

    static constexpr bool is_collection = false;
};

template <typename alphabet_t>
struct seqan3_sequence_pair_collection_generator
{
    using sequence_t = std::vector<alphabet_t>;

    auto operator()(benchmark::State & state) const
    {
        size_t sequence_length_variance = state.range(0);
        return seqan3::test::generate_sequence_pairs<alphabet_t>(sequence_length,
                                                                 set_size,
                                                                 sequence_length_variance);
    };

    size_t sequence_length;
    size_t set_size;

    static constexpr bool is_collection = true;
};

#ifdef SEQAN3_HAS_SEQAN2
template <typename alphabet_t>
struct seqan2_sequence_pair_generator
{
    using sequence_t = seqan::String<alphabet_t>;

    auto operator()(benchmark::State &) const
    {
        auto sequence1 = seqan3::test::generate_sequence_seqan2<alphabet_t>(sequence_length, 0, 0);
        auto sequence2 = seqan3::test::generate_sequence_seqan2<alphabet_t>(sequence_length, 0, 1);
        return std::tuple{std::move(sequence1), std::move(sequence2)};
    };

    size_t sequence_length;

    static constexpr bool is_collection = false;
};

template <typename alphabet_t>
struct seqan2_sequence_pair_collection_generator
{
    using sequence_t = seqan::String<alphabet_t>;

    auto operator()(benchmark::State & state) const
    {
        size_t sequence_length_variance = state.range(0);
        return seqan3::test::generate_sequence_pairs_seqan2<alphabet_t>(sequence_length,
                                                                        set_size,
                                                                        sequence_length_variance);
    };

    size_t sequence_length;
    size_t set_size;

    static constexpr bool is_collection = true;
};
#endif

// ----------------------------------------------------------------------------
//  seqan3 pairwise alignment benchmarks
// ----------------------------------------------------------------------------

template <typename sequence_pair_generator_t, typename align_configs_t>
void seqan3_align_pairwise_benchmark(benchmark::State & state,
                                     sequence_pair_generator_t sequence_pair_generator,
                                     align_configs_t && align_cfg)
{
    auto sequence_pair_or_pairs = sequence_pair_generator(state);
    constexpr bool collection_benchmark = sequence_pair_generator.is_collection;

    int64_t total = 0;
    for (auto _ : state)
    {
        for (auto && result : seqan3::align_pairwise(sequence_pair_or_pairs, align_cfg))
            total += result.score();
    }

    std::conditional_t<collection_benchmark, decltype(std::views::all), decltype(std::views::single)> view_adaptor{};
    state.counters["cells"] = seqan3::test::pairwise_cell_updates(view_adaptor(sequence_pair_or_pairs), align_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
    state.counters["total"] = total;
}

#ifdef SEQAN3_HAS_SEQAN2

// ----------------------------------------------------------------------------
//  seqan2 pairwise alignment benchmarks
// ----------------------------------------------------------------------------

template <typename sequence_pair_generator_t,
          typename align_cfg_t,
          typename scoring_scheme_t,
          typename execution_policy_t,
          typename seqan3_align_cfg_t>
void seqan2_align_pairwise_benchmark(benchmark::State & state,
                                     sequence_pair_generator_t sequence_pair_generator,
                                     align_cfg_t align_cfg,
                                     scoring_scheme_t scoring_scheme,
                                     execution_policy_t execution_policy,
                                     size_t thread_count,
                                     seqan3_align_cfg_t seqan3_align_cfg)
{
    using sequence_t = typename sequence_pair_generator_t::sequence_t;
    static constexpr bool collection_benchmark = sequence_pair_generator.is_collection;
    static constexpr bool with_alignment = seqan3_align_cfg_t::template exists<seqan3::align_cfg::output_alignment>();

    auto [sequences1, sequences2] = sequence_pair_generator(state);

    auto [gap_sequences1, gap_sequences2] = [&]()
    {
        using gapped_sequence_t [[maybe_unused]] = seqan::Gaps<sequence_t>;
        if constexpr (with_alignment && collection_benchmark)
        {
            seqan::StringSet<gapped_sequence_t> gaps1{};
            seqan::StringSet<gapped_sequence_t> gaps2{};

            for (unsigned i = 0; i < seqan::length(sequences1); ++i)
            {
                appendValue(gaps1, gapped_sequence_t{sequences1[i]});
                appendValue(gaps2, gapped_sequence_t{sequences2[i]});
            }

            return std::tuple{gaps1, gaps2};
        }
        else if constexpr (with_alignment && !collection_benchmark)
        {
            return std::tuple<gapped_sequence_t, gapped_sequence_t>{sequences1, sequences2};
        }
        else // if constexpr (!with_alignment)
        {
            return std::tuple{sequences1, sequences2};
        }
    }();

    auto algorithm = [&execution_policy, &scoring_scheme, &align_cfg](auto & sequences1, auto & sequences2)
    {
        if constexpr (with_alignment && collection_benchmark)
        {
            auto scores = seqan::globalAlignment(execution_policy, sequences1, sequences2, scoring_scheme, align_cfg);
            return std::accumulate(seqan::begin(scores), seqan::end(scores), 0);
        }
        else if constexpr (with_alignment && !collection_benchmark)
        {
            // NOTE execution_policy interface can't handle single sequences
            return seqan::globalAlignment(/*execution_policy,*/ sequences1, sequences2, scoring_scheme, align_cfg);
        }
        else if constexpr (!with_alignment && collection_benchmark)
        {
            auto scores = seqan::globalAlignmentScore(execution_policy, sequences1, sequences2, scoring_scheme, align_cfg);
            return std::accumulate(seqan::begin(scores), seqan::end(scores), 0);
        }
        else // if constexpr (!with_alignment && !collection_benchmark)
        {
            // NOTE execution_policy interface can't handle single sequences
            return seqan::globalAlignmentScore(/*execution_policy,*/ sequences1, sequences2, scoring_scheme, align_cfg);
        }
    };

    seqan::setNumThreads(execution_policy, thread_count);

    int64_t total = 0;
    for (auto _ : state)
        total += algorithm(gap_sequences1, gap_sequences2);

    std::conditional_t<collection_benchmark, decltype(std::views::all), decltype(std::views::single)> view_adaptor{};
    auto sequence_pairs_view = seqan3::views::zip(view_adaptor(sequences1), view_adaptor(sequences2));
    state.counters["cells"] = seqan3::test::pairwise_cell_updates(sequence_pairs_view, seqan3_align_cfg);
    state.counters["CUPS"] = seqan3::test::cell_updates_per_second(state.counters["cells"]);
    state.counters["total"] = total;
}

#endif // SEQAN3_HAS_SEQAN2
