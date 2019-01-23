// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <iostream>
#include <memory>
#include <random>
#include <utility>
#include <vector>

#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/aminoacid/aa20.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

#if __has_include(<seqan/align.h>)
    #define SEQAN3_HAS_SEQAN2 1
#endif

#ifdef SEQAN3_HAS_SEQAN2
    #include <seqan/align.h>
    #include <seqan/basic.h>
    #include <seqan/sequence.h>
#endif

using namespace seqan3;

template <typename alphabet_t>
auto generate_sequence_seqan3(size_t const len = 500,
                              size_t const variance = 0,
                              size_t const seed = 0)
{
    std::mt19937 gen(seed);
    std::uniform_int_distribution<uint8_t> dis_alpha(0, alphabet_size_v<alphabet_t>);
    std::uniform_int_distribution<size_t> dis_length(len - variance, len + variance);

    std::vector<alphabet_t> sequence;

    size_t length = dis_length(gen);
    for (size_t l = 0; l < length; ++l)
        sequence.push_back(alphabet_t{}.assign_rank(dis_alpha(gen)));

    return sequence;
}

#ifdef SEQAN3_HAS_SEQAN2
template <typename alphabet_t>
auto generate_sequence_seqan2(size_t const len = 500,
                              size_t const variance = 0,
                              size_t const seed = 0)
{
    std::mt19937 gen(seed);
    std::uniform_int_distribution<uint8_t> dis_alpha(0, seqan::ValueSize<alphabet_t>::VALUE);
    std::uniform_int_distribution<size_t> dis_length(len - variance, len + variance);

    seqan::String<alphabet_t> sequence;
    size_t length = dis_length(gen);
    for (size_t l = 0; l < length; ++l)
        appendValue(sequence, alphabet_t{dis_alpha(gen)});

    return sequence;
}
#endif // generate seqan2 data.

// ============================================================================
//  affine; score; dna4; single
// ============================================================================

void seqan3_affine_dna4(benchmark::State & state)
{
    auto cfg = align_cfg::mode{align_cfg::global_alignment} |
               align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}} |
               align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}} |
               align_cfg::result{align_cfg::with_score};

    auto seq1 = generate_sequence_seqan3<seqan3::dna4>(500, 0, 0);
    auto seq2 = generate_sequence_seqan3<seqan3::dna4>(500, 0, 1);

    for (auto _ : state)
    {
        auto rng = align_pairwise(std::tie(seq1, seq2), cfg);
        *seqan3::begin(rng);
    }
}

BENCHMARK(seqan3_affine_dna4);

#ifdef SEQAN3_HAS_SEQAN2

void seqan2_affine_dna4(benchmark::State & state)
{
    auto seq1 = generate_sequence_seqan2<seqan::Dna>(500, 0, 0);
    auto seq2 = generate_sequence_seqan2<seqan::Dna>(500, 0, 1);

    for (auto _ : state)
    {
        // In SeqAn2 the gap open contains already the gap extension costs, that's why we use -11 here.
        seqan::globalAlignmentScore(seq1, seq2, seqan::Score<int>{4, -5, -1, -11});
    }
}

BENCHMARK(seqan2_affine_dna4);
#endif // SEQAN3_HAS_SEQAN2

// ============================================================================
//  affine; score; dna4; set
// ============================================================================

void seqan3_affine_dna4_collection(benchmark::State & state)
{
    auto cfg = align_cfg::mode{align_cfg::global_alignment} |
               align_cfg::gap{gap_scheme{gap_score{-1}, gap_open_score{-10}}} |
               align_cfg::scoring{nucleotide_scoring_scheme{match_score{4}, mismatch_score{-5}}} |
               align_cfg::result{align_cfg::with_score};

    using sequence_t = decltype(generate_sequence_seqan3<seqan3::dna4>());

    std::vector<std::pair<sequence_t, sequence_t>> vec;
    for (unsigned i = 0; i < 100; ++i)
    {
        sequence_t seq1 = generate_sequence_seqan3<seqan3::dna4>(100, 0, i);
        sequence_t seq2 = generate_sequence_seqan3<seqan3::dna4>(100, 0, i + 100);
        vec.push_back(std::pair{seq1, seq2});
    }

    for (auto _ : state)
    {
        for (auto && rng : align_pairwise(vec, cfg))
            rng.get_score();
    }
}

BENCHMARK(seqan3_affine_dna4_collection);

#ifdef SEQAN3_HAS_SEQAN2

void seqan2_affine_dna4_collection(benchmark::State & state)
{
    using sequence_t = decltype(generate_sequence_seqan2<seqan::Dna>());

    seqan::StringSet<sequence_t> vec1;
    seqan::StringSet<sequence_t> vec2;
    for (unsigned i = 0; i < 100; ++i)
    {
        sequence_t seq1 = generate_sequence_seqan2<seqan::Dna>(100, 0, i);
        sequence_t seq2 = generate_sequence_seqan2<seqan::Dna>(100, 0, i + 100);
        appendValue(vec1, seq1);
        appendValue(vec2, seq2);
    }

    for (auto _ : state)
    {
        // In SeqAn2 the gap open contains already the gap extension costs, that's why we use -11 here.
        seqan::globalAlignmentScore(vec1, vec2, seqan::Score<int>{4, -5, -1, -11});
    }
}

BENCHMARK(seqan2_affine_dna4_collection);
#endif // SEQAN3_HAS_SEQAN2

// ============================================================================
//  instantiate tests
// ============================================================================

BENCHMARK_MAIN();
