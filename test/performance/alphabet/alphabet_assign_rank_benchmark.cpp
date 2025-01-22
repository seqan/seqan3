// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h>

#include <algorithm>
#include <cstring>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/test/performance/simd_dna4.hpp>
#include <seqan3/test/seqan2.hpp>

#if SEQAN3_HAS_SEQAN2
#    include <seqan/align.h>
#    include <seqan/basic.h>
#    include <seqan/modifier.h>
#endif

template <typename rank_t>
void fill_rank_array(std::array<rank_t, 256> & ranks, size_t const alphabet_size)
{
    uint8_t i = 0;
    for (rank_t & r : ranks)
        r = i++ % alphabet_size;
}

template <seqan3::semialphabet alphabet_t>
void assign_rank(benchmark::State & state)
{
    using rank_t = seqan3::alphabet_rank_t<alphabet_t>;

    std::array<rank_t, 256> ranks{};
    fill_rank_array<rank_t>(ranks, seqan3::alphabet_size<alphabet_t>);

    alphabet_t a{};
    for (auto _ : state)
        for (rank_t r : ranks)
            benchmark::DoNotOptimize(seqan3::assign_rank_to(r, a));
}

/* regular alphabets, sorted by size */
BENCHMARK_TEMPLATE(assign_rank, seqan3::gap);
BENCHMARK_TEMPLATE(assign_rank, seqan3::dna4);
BENCHMARK_TEMPLATE(assign_rank, seqan3::rna4);
BENCHMARK_TEMPLATE(assign_rank, seqan3::simd_dna4);
BENCHMARK_TEMPLATE(assign_rank, seqan3::dna5);
BENCHMARK_TEMPLATE(assign_rank, seqan3::rna5);
BENCHMARK_TEMPLATE(assign_rank, seqan3::dna15);
BENCHMARK_TEMPLATE(assign_rank, seqan3::rna15);
BENCHMARK_TEMPLATE(assign_rank, seqan3::aa20);
BENCHMARK_TEMPLATE(assign_rank, seqan3::aa27);
BENCHMARK_TEMPLATE(assign_rank, seqan3::phred42);
BENCHMARK_TEMPLATE(assign_rank, seqan3::phred63);
BENCHMARK_TEMPLATE(assign_rank, seqan3::phred94);
/* adaptations */
BENCHMARK_TEMPLATE(assign_rank, char);
BENCHMARK_TEMPLATE(assign_rank, char32_t);
/* alphabet variant */
BENCHMARK_TEMPLATE(assign_rank, seqan3::gapped<seqan3::dna4>);
BENCHMARK_TEMPLATE(assign_rank,
                   seqan3::alphabet_variant<seqan3::gap,
                                            seqan3::dna4,
                                            seqan3::dna5,
                                            seqan3::dna15,
                                            seqan3::rna15,
                                            seqan3::rna4,
                                            seqan3::rna5>);
BENCHMARK_TEMPLATE(assign_rank, seqan3::alphabet_variant<seqan3::dna4, char>);
/* alphabet tuple */
BENCHMARK_TEMPLATE(assign_rank, seqan3::masked<seqan3::dna4>);
BENCHMARK_TEMPLATE(assign_rank, seqan3::qualified<seqan3::dna4, seqan3::phred42>);
BENCHMARK_TEMPLATE(assign_rank, seqan3::qualified<seqan3::dna5, seqan3::phred63>);
BENCHMARK_TEMPLATE(assign_rank, seqan3::qualified<seqan3::dna5, seqan3::phred94>);

#if SEQAN3_HAS_SEQAN2
template <typename alphabet_t>
void assign_rank_seqan2(benchmark::State & state)
{
    using rank_t = uint8_t;

    std::array<rank_t, 256> ranks{};
    fill_rank_array<rank_t>(ranks, seqan2::ValueSize<alphabet_t>::VALUE);

    alphabet_t a{};
    for (auto _ : state)
    {
        for (rank_t r : ranks)
        {
            a = r;
            benchmark::DoNotOptimize(a);
        }
    }
}

BENCHMARK_TEMPLATE(assign_rank_seqan2, seqan2::Dna);
BENCHMARK_TEMPLATE(assign_rank_seqan2, seqan2::Rna);
BENCHMARK_TEMPLATE(assign_rank_seqan2, seqan2::Dna5);
BENCHMARK_TEMPLATE(assign_rank_seqan2, seqan2::Rna5);
BENCHMARK_TEMPLATE(assign_rank_seqan2, seqan2::Iupac);
BENCHMARK_TEMPLATE(assign_rank_seqan2, seqan2::AminoAcid);

BENCHMARK_TEMPLATE(assign_rank_seqan2, seqan2::Dna5Q);
BENCHMARK_TEMPLATE(assign_rank_seqan2, typename seqan2::GappedValueType<seqan2::Dna>::Type);
#endif

BENCHMARK_MAIN();
