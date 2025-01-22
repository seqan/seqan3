// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h>

#include <algorithm>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/test/performance/simd_dna4.hpp>
#include <seqan3/test/seqan2.hpp>

#if SEQAN3_HAS_SEQAN2
#    include <seqan/align.h>
#    include <seqan/basic.h>
#    include <seqan/modifier.h>
#endif

template <typename alphabet_t, bool is_seqan2>
std::array<alphabet_t, 256> create_alphabet_array(size_t const alphabet_size)
{
    std::array<alphabet_t, 256> alphabet_array;

    auto convert_to_alphabet = [](auto const rank, auto & a)
    {
        if constexpr (is_seqan2)
            a = static_cast<uint8_t>(rank);
        else
            seqan3::assign_rank_to(rank, a);
    };

    uint8_t i = 0;
    for (alphabet_t & a : alphabet_array)
        convert_to_alphabet(i++ % alphabet_size, a);

    return alphabet_array;
}

template <seqan3::semialphabet alphabet_t>
void to_rank(benchmark::State & state)
{
    std::array<alphabet_t, 256> alphs = create_alphabet_array<alphabet_t, false>(seqan3::alphabet_size<alphabet_t>);

    for (auto _ : state)
    {
        for (alphabet_t a : alphs)
        {
            auto rank = seqan3::to_rank(a);
            benchmark::DoNotOptimize(rank);
        }
    }
}

/* regular alphabets, sorted by size */
BENCHMARK_TEMPLATE(to_rank, seqan3::gap);
BENCHMARK_TEMPLATE(to_rank, seqan3::dna4);
BENCHMARK_TEMPLATE(to_rank, seqan3::rna4);
BENCHMARK_TEMPLATE(to_rank, seqan3::simd_dna4);
BENCHMARK_TEMPLATE(to_rank, seqan3::dna5);
BENCHMARK_TEMPLATE(to_rank, seqan3::rna5);
BENCHMARK_TEMPLATE(to_rank, seqan3::dna15);
BENCHMARK_TEMPLATE(to_rank, seqan3::rna15);
BENCHMARK_TEMPLATE(to_rank, seqan3::aa20);
BENCHMARK_TEMPLATE(to_rank, seqan3::aa27);
BENCHMARK_TEMPLATE(to_rank, seqan3::phred42);
BENCHMARK_TEMPLATE(to_rank, seqan3::phred63);
BENCHMARK_TEMPLATE(to_rank, seqan3::phred94);
/* adaptations */
BENCHMARK_TEMPLATE(to_rank, char);
BENCHMARK_TEMPLATE(to_rank, char32_t);
/* alphabet variant */
BENCHMARK_TEMPLATE(to_rank, seqan3::gapped<seqan3::dna4>);
BENCHMARK_TEMPLATE(to_rank,
                   seqan3::alphabet_variant<seqan3::gap,
                                            seqan3::dna4,
                                            seqan3::dna5,
                                            seqan3::dna15,
                                            seqan3::rna15,
                                            seqan3::rna4,
                                            seqan3::rna5>);
BENCHMARK_TEMPLATE(to_rank, seqan3::alphabet_variant<seqan3::dna4, char>);
/* alphabet tuple */
BENCHMARK_TEMPLATE(to_rank, seqan3::masked<seqan3::dna4>);
BENCHMARK_TEMPLATE(to_rank, seqan3::qualified<seqan3::dna4, seqan3::phred42>);
BENCHMARK_TEMPLATE(to_rank, seqan3::qualified<seqan3::dna5, seqan3::phred63>);
BENCHMARK_TEMPLATE(to_rank, seqan3::qualified<seqan3::dna5, seqan3::phred94>);

#if SEQAN3_HAS_SEQAN2
template <typename alphabet_t>
void to_rank_seqan2(benchmark::State & state)
{
    std::array<alphabet_t, 256> alphs = create_alphabet_array<alphabet_t, true>(seqan2::ValueSize<alphabet_t>::VALUE);

    for (auto _ : state)
    {
        for (alphabet_t a : alphs)
        {
            auto ord_value = seqan2::ordValue(a);
            benchmark::DoNotOptimize(ord_value);
        }
    }
}

BENCHMARK_TEMPLATE(to_rank_seqan2, seqan2::Dna);
BENCHMARK_TEMPLATE(to_rank_seqan2, seqan2::Rna);
BENCHMARK_TEMPLATE(to_rank_seqan2, seqan2::Dna5);
BENCHMARK_TEMPLATE(to_rank_seqan2, seqan2::Rna5);
BENCHMARK_TEMPLATE(to_rank_seqan2, seqan2::Iupac);
BENCHMARK_TEMPLATE(to_rank_seqan2, seqan2::AminoAcid);

BENCHMARK_TEMPLATE(to_rank_seqan2, seqan2::Dna5Q);
BENCHMARK_TEMPLATE(to_rank_seqan2, typename seqan2::GappedValueType<seqan2::Dna>::Type);
#endif

BENCHMARK_MAIN();
