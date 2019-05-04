// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>
#include <cstring>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/search/fm_index/all.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/seqan2.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/range/view/char_to.hpp>

using namespace seqan3;

unsigned int const SEED = 1234;

// ============================================================================
// construct as many indices as possible
// ============================================================================
template<typename index_t, std::size_t const length>
void construct_index(benchmark::State & state)
{
    using alphabet_t = typename index_t::char_type;

    auto sequence = test::generate_sequence<alphabet_t>(length, 0, SEED);

    for (auto _ : state)
        index_t index{sequence};
}

// show effect of alphabet size and compare fm with bi-fm
BENCHMARK_TEMPLATE(construct_index, fm_index<std::vector<dna4>>, 25000);
BENCHMARK_TEMPLATE(construct_index, fm_index<std::vector<aa27>>, 25000);
BENCHMARK_TEMPLATE(construct_index, fm_index<std::vector<phred63>>, 25000);

BENCHMARK_TEMPLATE(construct_index, bi_fm_index<std::vector<dna4>>, 25000);
BENCHMARK_TEMPLATE(construct_index, bi_fm_index<std::vector<aa27>>, 25000);
BENCHMARK_TEMPLATE(construct_index, bi_fm_index<std::vector<phred63>>, 25000);

#if SEQAN3_HAS_SEQAN2

// ============================================================================
// export generated sequences for seqan2
// ============================================================================
std::map<std::size_t, std::vector<char>> dna_sequence_dict;

std::vector<char> get_sequence(std::size_t length)
{
    if (dna_sequence_dict.find(length) == dna_sequence_dict.end())
    {
        // generate new
        auto sequence = test::generate_sequence<dna4>(length, 0, SEED);
        dna_sequence_dict.insert({length, sequence | view::to_char});
    }

    return dna_sequence_dict[length];
}

// ============================================================================
// seqan3 index construction with dna4 only
// ============================================================================
template<typename index_t, std::size_t const length>
void construct_index_seqan3(benchmark::State & state)
{

    auto sequence_char = get_sequence(length);
    std::vector<dna4> sequence = sequence_char | view::char_to<dna4>;

    for (auto _ : state)
        index_t index{sequence};
}

// ============================================================================
// seqan2 comparison using exact same input as seqan3 above
// ============================================================================
#include <seqan/index.h>

using namespace seqan;

template<typename index_t, const std::size_t length>
void construct_index_seqan2(benchmark::State & state)
{
    std::vector<char> imported_seq = get_sequence(length);
    String<Dna> sequence{};

    for (auto c : imported_seq)
        sequence += c;

    for (auto _ : state)
    {
        index_t index(sequence);
        indexRequire(index, FibreSA());
    }
}

typedef FastFMIndexConfig<void, uint32_t> cfg;

BENCHMARK_TEMPLATE(construct_index_seqan3, fm_index<std::vector<dna4>>, 5);
BENCHMARK_TEMPLATE(construct_index_seqan3, fm_index<std::vector<dna4>>, 50);
BENCHMARK_TEMPLATE(construct_index_seqan3, fm_index<std::vector<dna4>>, 500);
BENCHMARK_TEMPLATE(construct_index_seqan3, fm_index<std::vector<dna4>>, 5000);
BENCHMARK_TEMPLATE(construct_index_seqan3, fm_index<std::vector<dna4>>, 50000);
BENCHMARK_TEMPLATE(construct_index_seqan3, fm_index<std::vector<dna4>>, 500000);
BENCHMARK_TEMPLATE(construct_index_seqan3, fm_index<std::vector<dna4>>, 5000000);

BENCHMARK_TEMPLATE(construct_index_seqan2, Index<String<Dna>, FMIndex<void, cfg>>, 5);
BENCHMARK_TEMPLATE(construct_index_seqan2, Index<String<Dna>, FMIndex<void, cfg>>, 50);
BENCHMARK_TEMPLATE(construct_index_seqan2, Index<String<Dna>, FMIndex<void, cfg>>, 500);
BENCHMARK_TEMPLATE(construct_index_seqan2, Index<String<Dna>, FMIndex<void, cfg>>, 5000);
BENCHMARK_TEMPLATE(construct_index_seqan2, Index<String<Dna>, FMIndex<void, cfg>>, 50000);
BENCHMARK_TEMPLATE(construct_index_seqan2, Index<String<Dna>, FMIndex<void, cfg>>, 500000);
BENCHMARK_TEMPLATE(construct_index_seqan2, Index<String<Dna>, FMIndex<void, cfg>>, 5000000);

#endif // seqan2 included

BENCHMARK_MAIN();
