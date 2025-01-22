// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h>

#include <algorithm>
#include <cstring>
#include <numeric>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/test/performance/simd_dna4.hpp>
#include <seqan3/test/seqan2.hpp>

#if SEQAN3_HAS_SEQAN2
#    include <seqan/align.h>
#    include <seqan/basic.h>
#    include <seqan/modifier.h>
#endif

template <seqan3::alphabet alphabet_t>
void assign_char(benchmark::State & state)
{
    using char_t = seqan3::alphabet_char_t<alphabet_t>;

    std::array<char_t, 256> chars{};
    std::iota(chars.begin(), chars.end(), 0);

    alphabet_t a{};
    for (auto _ : state)
        for (char_t c : chars)
            benchmark::DoNotOptimize(seqan3::assign_char_to(c, a));
}

/* regular alphabets, sorted by size */
BENCHMARK_TEMPLATE(assign_char, seqan3::gap);
BENCHMARK_TEMPLATE(assign_char, seqan3::dna4);
BENCHMARK_TEMPLATE(assign_char, seqan3::rna4);
BENCHMARK_TEMPLATE(assign_char, seqan3::simd_dna4);
BENCHMARK_TEMPLATE(assign_char, seqan3::dna5);
BENCHMARK_TEMPLATE(assign_char, seqan3::rna5);
BENCHMARK_TEMPLATE(assign_char, seqan3::dna15);
BENCHMARK_TEMPLATE(assign_char, seqan3::rna15);
BENCHMARK_TEMPLATE(assign_char, seqan3::aa20);
BENCHMARK_TEMPLATE(assign_char, seqan3::aa27);
BENCHMARK_TEMPLATE(assign_char, seqan3::phred42);
BENCHMARK_TEMPLATE(assign_char, seqan3::phred63);
BENCHMARK_TEMPLATE(assign_char, seqan3::phred94);
/* adaptations */
BENCHMARK_TEMPLATE(assign_char, char);
BENCHMARK_TEMPLATE(assign_char, char32_t);
/* alphabet variant */
BENCHMARK_TEMPLATE(assign_char, seqan3::gapped<seqan3::dna4>);
BENCHMARK_TEMPLATE(assign_char,
                   seqan3::alphabet_variant<seqan3::gap,
                                            seqan3::dna4,
                                            seqan3::dna5,
                                            seqan3::dna15,
                                            seqan3::rna15,
                                            seqan3::rna4,
                                            seqan3::rna5>);
BENCHMARK_TEMPLATE(assign_char, seqan3::alphabet_variant<seqan3::dna4, char>);
/* alphabet tuple */
BENCHMARK_TEMPLATE(assign_char, seqan3::masked<seqan3::dna4>);
BENCHMARK_TEMPLATE(assign_char, seqan3::qualified<seqan3::dna4, seqan3::phred42>);
BENCHMARK_TEMPLATE(assign_char, seqan3::qualified<seqan3::dna5, seqan3::phred63>);
BENCHMARK_TEMPLATE(assign_char, seqan3::qualified<seqan3::dna5, seqan3::phred94>);

#if SEQAN3_HAS_SEQAN2
template <typename alphabet_t>
void assign_char_seqan2(benchmark::State & state)
{
    std::array<char, 256> chars{};
    std::iota(chars.begin(), chars.end(), 0);

    alphabet_t a{};
    for (auto _ : state)
    {
        for (char c : chars)
        {
            a = c;
            benchmark::DoNotOptimize(a);
        }
    }
}

BENCHMARK_TEMPLATE(assign_char_seqan2, seqan2::Dna);
BENCHMARK_TEMPLATE(assign_char_seqan2, seqan2::Rna);
BENCHMARK_TEMPLATE(assign_char_seqan2, seqan2::Dna5);
BENCHMARK_TEMPLATE(assign_char_seqan2, seqan2::Rna5);
BENCHMARK_TEMPLATE(assign_char_seqan2, seqan2::Iupac);
BENCHMARK_TEMPLATE(assign_char_seqan2, seqan2::AminoAcid);

BENCHMARK_TEMPLATE(assign_char_seqan2, seqan2::Dna5Q);
BENCHMARK_TEMPLATE(assign_char_seqan2, typename seqan2::GappedValueType<seqan2::Dna>::Type);
#endif

BENCHMARK_MAIN();
