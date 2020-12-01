// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <algorithm>
#include <cstring>

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/test/seqan2.hpp>

#if SEQAN3_HAS_SEQAN2
#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/modifier.h>
#endif

template <typename alphabet_t, bool is_seqan2>
std::array<alphabet_t, 256> create_alphabet_array()
{
    std::array<alphabet_t, 256> alphabet_array{};

    auto convert_to_alphabet = [] (auto const c, auto & a)
    {
        if constexpr (is_seqan2)
        {
             a = (char) c;
        }
        else
        {
            using char_t = seqan3::alphabet_char_t<alphabet_t>;
            seqan3::assign_char_to(char_t(c), a);
        }
    };

    uint8_t i = 0;
    for (alphabet_t & a : alphabet_array)
        convert_to_alphabet(i++, a);

    return alphabet_array;
}

template <seqan3::alphabet alphabet_t>
void to_char(benchmark::State & state)
{
    std::array<alphabet_t, 256> alphs = create_alphabet_array<alphabet_t, false>();

    for (auto _ : state)
        for (alphabet_t a : alphs)
            benchmark::DoNotOptimize(seqan3::to_char(a));
}

/* regular alphabets, sorted by size */
BENCHMARK_TEMPLATE(to_char, seqan3::gap);
BENCHMARK_TEMPLATE(to_char, seqan3::dna4);
BENCHMARK_TEMPLATE(to_char, seqan3::rna4);
BENCHMARK_TEMPLATE(to_char, seqan3::dna5);
BENCHMARK_TEMPLATE(to_char, seqan3::rna5);
BENCHMARK_TEMPLATE(to_char, seqan3::dna15);
BENCHMARK_TEMPLATE(to_char, seqan3::rna15);
BENCHMARK_TEMPLATE(to_char, seqan3::aa20);
BENCHMARK_TEMPLATE(to_char, seqan3::aa27);
BENCHMARK_TEMPLATE(to_char, seqan3::phred42);
BENCHMARK_TEMPLATE(to_char, seqan3::phred63);
BENCHMARK_TEMPLATE(to_char, seqan3::phred94);
/* adaptations */
BENCHMARK_TEMPLATE(to_char, char);
BENCHMARK_TEMPLATE(to_char, char32_t);
/* alphabet variant */
BENCHMARK_TEMPLATE(to_char, seqan3::gapped<seqan3::dna4>);
BENCHMARK_TEMPLATE(to_char, seqan3::alphabet_variant<seqan3::gap, seqan3::dna4, seqan3::dna5, seqan3::dna15,
                                                     seqan3::rna15, seqan3::rna4, seqan3::rna5>);
BENCHMARK_TEMPLATE(to_char, seqan3::alphabet_variant<seqan3::dna4, char>);
/* alphabet tuple */
BENCHMARK_TEMPLATE(to_char, seqan3::masked<seqan3::dna4>);
BENCHMARK_TEMPLATE(to_char, seqan3::qualified<seqan3::dna4, seqan3::phred42>);
BENCHMARK_TEMPLATE(to_char, seqan3::qualified<seqan3::dna5, seqan3::phred63>);
BENCHMARK_TEMPLATE(to_char, seqan3::qualified<seqan3::dna5, seqan3::phred94>);

#if SEQAN3_HAS_SEQAN2
template <typename alphabet_t>
void to_char_seqan2(benchmark::State & state)
{
    std::array<alphabet_t, 256> alphs = create_alphabet_array<alphabet_t, true>();

    for (auto _ : state)
        for (alphabet_t a : alphs)
            benchmark::DoNotOptimize(static_cast<char>(a));
}

BENCHMARK_TEMPLATE(to_char_seqan2, seqan::Dna);
BENCHMARK_TEMPLATE(to_char_seqan2, seqan::Rna);
BENCHMARK_TEMPLATE(to_char_seqan2, seqan::Dna5);
BENCHMARK_TEMPLATE(to_char_seqan2, seqan::Rna5);
BENCHMARK_TEMPLATE(to_char_seqan2, seqan::Iupac);
BENCHMARK_TEMPLATE(to_char_seqan2, seqan::AminoAcid);

BENCHMARK_TEMPLATE(to_char_seqan2, seqan::Dna5Q);
BENCHMARK_TEMPLATE(to_char_seqan2, typename seqan::GappedValueType<seqan::Dna>::Type);
#endif

BENCHMARK_MAIN();
