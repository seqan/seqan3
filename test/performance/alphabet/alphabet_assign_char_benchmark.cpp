// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <algorithm>
#include <cstring>
#include <numeric>

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/test/seqan2.hpp>

#if SEQAN3_HAS_SEQAN2
#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/modifier.h>
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
BENCHMARK_TEMPLATE(assign_char, seqan3::alphabet_variant<seqan3::gap, seqan3::dna4, seqan3::dna5, seqan3::dna15,
                                                         seqan3::rna15, seqan3::rna4, seqan3::rna5>);
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
        for (char c : chars)
            benchmark::DoNotOptimize(a = c);
}

BENCHMARK_TEMPLATE(assign_char_seqan2, seqan::Dna);
BENCHMARK_TEMPLATE(assign_char_seqan2, seqan::Rna);
BENCHMARK_TEMPLATE(assign_char_seqan2, seqan::Dna5);
BENCHMARK_TEMPLATE(assign_char_seqan2, seqan::Rna5);
BENCHMARK_TEMPLATE(assign_char_seqan2, seqan::Iupac);
BENCHMARK_TEMPLATE(assign_char_seqan2, seqan::AminoAcid);

BENCHMARK_TEMPLATE(assign_char_seqan2, seqan::Dna5Q);
BENCHMARK_TEMPLATE(assign_char_seqan2, typename seqan::GappedValueType<seqan::Dna>::Type);
#endif

BENCHMARK_MAIN();
