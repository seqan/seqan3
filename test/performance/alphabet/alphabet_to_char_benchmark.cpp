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

template <seqan3::alphabet alphabet_t>
void to_char_(benchmark::State & state)
{
    using char_t = seqan3::alphabet_char_t<alphabet_t>;

    std::array<alphabet_t, 256> alphs{};
    size_t i = 0;
    for (alphabet_t & a : alphs)
        seqan3::assign_char_to(char_t(i++), a);

    for (auto _ : state)
        for (alphabet_t a : alphs)
            benchmark::DoNotOptimize(seqan3::to_char(a));
}

/* regular alphabets, sorted by size */
BENCHMARK_TEMPLATE(to_char_, seqan3::gap);
BENCHMARK_TEMPLATE(to_char_, seqan3::dna4);
BENCHMARK_TEMPLATE(to_char_, seqan3::rna4);
BENCHMARK_TEMPLATE(to_char_, seqan3::dna5);
BENCHMARK_TEMPLATE(to_char_, seqan3::rna5);
BENCHMARK_TEMPLATE(to_char_, seqan3::dna15);
BENCHMARK_TEMPLATE(to_char_, seqan3::rna15);
BENCHMARK_TEMPLATE(to_char_, seqan3::aa20);
BENCHMARK_TEMPLATE(to_char_, seqan3::aa27);
BENCHMARK_TEMPLATE(to_char_, seqan3::phred42);
BENCHMARK_TEMPLATE(to_char_, seqan3::phred63);
/* adaptations */
BENCHMARK_TEMPLATE(to_char_, char);
BENCHMARK_TEMPLATE(to_char_, char32_t);
/* alphabet variant */
BENCHMARK_TEMPLATE(to_char_, seqan3::gapped<seqan3::dna4>);
BENCHMARK_TEMPLATE(to_char_, seqan3::alphabet_variant<seqan3::gap, seqan3::dna4, seqan3::dna5, seqan3::dna15,
                                                      seqan3::rna15, seqan3::rna4, seqan3::rna5>);
BENCHMARK_TEMPLATE(to_char_, seqan3::alphabet_variant<seqan3::dna4, char>);
/* alphabet tuple */
BENCHMARK_TEMPLATE(to_char_, seqan3::masked<seqan3::dna4>);
BENCHMARK_TEMPLATE(to_char_, seqan3::qualified<seqan3::dna4, seqan3::phred42>);
BENCHMARK_TEMPLATE(to_char_, seqan3::qualified<seqan3::dna5, seqan3::phred63>);

BENCHMARK_MAIN();
