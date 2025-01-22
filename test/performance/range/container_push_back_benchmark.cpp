// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h>

#include <deque>
#include <list>
#include <vector>

#include <sdsl/int_vector.hpp>

#include <seqan3/alignment/decorator/gap_decorator.hpp>
#include <seqan3/alphabet/all.hpp>
#include <seqan3/alphabet/container/bitpacked_sequence.hpp>
#include <seqan3/test/seqan2.hpp>
#include <seqan3/utility/container/small_vector.hpp>

template <typename t>
using sdsl_int_vec = sdsl::int_vector<sizeof(t) * 8>;

template <typename t>
using small_vec = seqan3::small_vector<t, 10000>;

// ============================================================================
//  push_back
// ============================================================================

template <template <typename...> typename container_t, typename alphabet_t, typename... args>
void push_back(benchmark::State & state)
{
    alphabet_t letter{};

    for (auto _ : state)
    {
        container_t<alphabet_t, args...> container;
        for (size_t i = 0; i < 10000; ++i)
            container.push_back(letter);
        benchmark::DoNotOptimize(letter = container.back());
    }

    state.counters["sizeof"] = sizeof(alphabet_t);
    if constexpr (seqan3::alphabet<alphabet_t>)
        state.counters["alph_size"] = seqan3::alphabet_size<alphabet_t>;
}

BENCHMARK_TEMPLATE(push_back, std::vector, char);
BENCHMARK_TEMPLATE(push_back, std::vector, uint8_t);
BENCHMARK_TEMPLATE(push_back, std::vector, uint16_t);
BENCHMARK_TEMPLATE(push_back, std::vector, uint32_t);
BENCHMARK_TEMPLATE(push_back, std::vector, uint64_t);
BENCHMARK_TEMPLATE(push_back, std::vector, seqan3::gap);
BENCHMARK_TEMPLATE(push_back, std::vector, seqan3::dna4);
BENCHMARK_TEMPLATE(push_back, std::vector, seqan3::gapped<seqan3::dna4>);
BENCHMARK_TEMPLATE(push_back, std::vector, seqan3::dna15);
BENCHMARK_TEMPLATE(push_back, std::vector, seqan3::aa27);
BENCHMARK_TEMPLATE(push_back, std::vector, seqan3::alphabet_variant<char, seqan3::dna4>);

BENCHMARK_TEMPLATE(push_back, std::deque, char);
BENCHMARK_TEMPLATE(push_back, std::deque, uint8_t);
BENCHMARK_TEMPLATE(push_back, std::deque, uint16_t);
BENCHMARK_TEMPLATE(push_back, std::deque, uint32_t);
BENCHMARK_TEMPLATE(push_back, std::deque, uint64_t);
BENCHMARK_TEMPLATE(push_back, std::deque, seqan3::gap);
BENCHMARK_TEMPLATE(push_back, std::deque, seqan3::dna4);
BENCHMARK_TEMPLATE(push_back, std::deque, seqan3::gapped<seqan3::dna4>);
BENCHMARK_TEMPLATE(push_back, std::deque, seqan3::dna15);
BENCHMARK_TEMPLATE(push_back, std::deque, seqan3::aa27);
BENCHMARK_TEMPLATE(push_back, std::deque, seqan3::alphabet_variant<char, seqan3::dna4>);

BENCHMARK_TEMPLATE(push_back, std::list, char);
BENCHMARK_TEMPLATE(push_back, std::list, uint8_t);
BENCHMARK_TEMPLATE(push_back, std::list, uint16_t);
BENCHMARK_TEMPLATE(push_back, std::list, uint32_t);
BENCHMARK_TEMPLATE(push_back, std::list, uint64_t);
BENCHMARK_TEMPLATE(push_back, std::list, seqan3::gap);
BENCHMARK_TEMPLATE(push_back, std::list, seqan3::dna4);
BENCHMARK_TEMPLATE(push_back, std::list, seqan3::gapped<seqan3::dna4>);
BENCHMARK_TEMPLATE(push_back, std::list, seqan3::dna15);
BENCHMARK_TEMPLATE(push_back, std::list, seqan3::aa27);
BENCHMARK_TEMPLATE(push_back, std::list, seqan3::alphabet_variant<char, seqan3::dna4>);

BENCHMARK_TEMPLATE(push_back, sdsl_int_vec, uint8_t);
BENCHMARK_TEMPLATE(push_back, sdsl_int_vec, uint16_t);
BENCHMARK_TEMPLATE(push_back, sdsl_int_vec, uint32_t);
BENCHMARK_TEMPLATE(push_back, sdsl_int_vec, uint64_t);

BENCHMARK_TEMPLATE(push_back, seqan3::bitpacked_sequence, char);
BENCHMARK_TEMPLATE(push_back, seqan3::bitpacked_sequence, uint32_t);
BENCHMARK_TEMPLATE(push_back, seqan3::bitpacked_sequence, seqan3::gap);
BENCHMARK_TEMPLATE(push_back, seqan3::bitpacked_sequence, seqan3::dna4);
BENCHMARK_TEMPLATE(push_back, seqan3::bitpacked_sequence, seqan3::gapped<seqan3::dna4>);
BENCHMARK_TEMPLATE(push_back, seqan3::bitpacked_sequence, seqan3::dna15);
BENCHMARK_TEMPLATE(push_back, seqan3::bitpacked_sequence, seqan3::aa27);
BENCHMARK_TEMPLATE(push_back, seqan3::bitpacked_sequence, seqan3::alphabet_variant<char, seqan3::dna4>);

BENCHMARK_TEMPLATE(push_back, small_vec, char);
BENCHMARK_TEMPLATE(push_back, small_vec, uint32_t);
BENCHMARK_TEMPLATE(push_back, small_vec, seqan3::gap);
BENCHMARK_TEMPLATE(push_back, small_vec, seqan3::dna4);
BENCHMARK_TEMPLATE(push_back, small_vec, seqan3::gapped<seqan3::dna4>);
BENCHMARK_TEMPLATE(push_back, small_vec, seqan3::dna15);
BENCHMARK_TEMPLATE(push_back, small_vec, seqan3::aa27);
BENCHMARK_TEMPLATE(push_back, small_vec, seqan3::alphabet_variant<char, seqan3::dna4>);

// ============================================================================
//  push_back SeqAn2
// ============================================================================

#if SEQAN3_HAS_SEQAN2

#    include <seqan/sequence.h>

template <template <typename...> typename container_t, typename spec_t, typename alphabet_t>
void push_back2(benchmark::State & state)
{
    alphabet_t letter{};

    for (auto _ : state)
    {
        container_t<alphabet_t, spec_t> container;
        for (size_t i = 0; i < 10000; ++i)
            seqan2::appendValue(container, letter);
        benchmark::DoNotOptimize(letter = seqan2::back(container));
    }

    state.counters["sizeof"] = sizeof(alphabet_t);
    state.counters["alph_size"] = seqan2::ValueSize<alphabet_t>::VALUE;
}

BENCHMARK_TEMPLATE(push_back, std::vector, seqan2::Dna);
BENCHMARK_TEMPLATE(push_back, std::vector, seqan2::Dna5);
BENCHMARK_TEMPLATE(push_back, std::vector, seqan2::Iupac);
BENCHMARK_TEMPLATE(push_back, std::vector, seqan2::AminoAcid);
BENCHMARK_TEMPLATE(push_back, std::vector, seqan2::Dna5Q);

BENCHMARK_TEMPLATE(push_back2, seqan2::String, seqan2::Alloc<>, char);
BENCHMARK_TEMPLATE(push_back2, seqan2::String, seqan2::Alloc<>, uint8_t);
BENCHMARK_TEMPLATE(push_back2, seqan2::String, seqan2::Alloc<>, uint16_t);
BENCHMARK_TEMPLATE(push_back2, seqan2::String, seqan2::Alloc<>, uint32_t);
BENCHMARK_TEMPLATE(push_back2, seqan2::String, seqan2::Alloc<>, uint64_t);

BENCHMARK_TEMPLATE(push_back2, seqan2::String, seqan2::Alloc<>, seqan3::gap);
BENCHMARK_TEMPLATE(push_back2, seqan2::String, seqan2::Alloc<>, seqan3::dna4);
BENCHMARK_TEMPLATE(push_back2, seqan2::String, seqan2::Alloc<>, seqan3::gapped<seqan3::dna4>);
BENCHMARK_TEMPLATE(push_back2, seqan2::String, seqan2::Alloc<>, seqan3::dna15);
BENCHMARK_TEMPLATE(push_back2, seqan2::String, seqan2::Alloc<>, seqan3::aa27);
BENCHMARK_TEMPLATE(push_back2, seqan2::String, seqan2::Alloc<>, seqan3::alphabet_variant<char, seqan3::dna4>);

BENCHMARK_TEMPLATE(push_back2, seqan2::String, seqan2::Alloc<>, seqan2::Dna);
BENCHMARK_TEMPLATE(push_back2, seqan2::String, seqan2::Alloc<>, seqan2::Dna5);
BENCHMARK_TEMPLATE(push_back2, seqan2::String, seqan2::Alloc<>, seqan2::Iupac);
BENCHMARK_TEMPLATE(push_back2, seqan2::String, seqan2::Alloc<>, seqan2::AminoAcid);
BENCHMARK_TEMPLATE(push_back2, seqan2::String, seqan2::Alloc<>, seqan2::Dna5Q);

BENCHMARK_TEMPLATE(push_back2, seqan2::String, seqan2::Packed<>, seqan2::Dna);
BENCHMARK_TEMPLATE(push_back2, seqan2::String, seqan2::Packed<>, seqan2::Dna5);
BENCHMARK_TEMPLATE(push_back2, seqan2::String, seqan2::Packed<>, seqan2::Iupac);
BENCHMARK_TEMPLATE(push_back2, seqan2::String, seqan2::Packed<>, seqan2::AminoAcid);
// BENCHMARK_TEMPLATE(push_back2, seqan2::String, seqan2::Packed<>, seqan2::Dna5Q); // broken in SeqAn2

BENCHMARK_TEMPLATE(push_back2, seqan2::String, seqan2::Array<10000>, seqan2::Dna);
BENCHMARK_TEMPLATE(push_back2, seqan2::String, seqan2::Array<10000>, seqan2::Dna5);
BENCHMARK_TEMPLATE(push_back2, seqan2::String, seqan2::Array<10000>, seqan2::Iupac);
BENCHMARK_TEMPLATE(push_back2, seqan2::String, seqan2::Array<10000>, seqan2::AminoAcid);
BENCHMARK_TEMPLATE(push_back2, seqan2::String, seqan2::Array<10000>, seqan2::Dna5Q);

#endif // SEQAN3_HAS_SEQAN2

// ============================================================================
//  run
// ============================================================================

BENCHMARK_MAIN();
