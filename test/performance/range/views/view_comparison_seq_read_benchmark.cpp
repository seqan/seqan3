// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <list>
#include <vector>

#include <benchmark/benchmark.h>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/range/views/all.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/seqan2.hpp>

enum tag
{
    baseline,
    slice,
    comp,
    rev,
    rev_comp,
    to_char,
    to_upper,
    filter,
    pseudofilter,
    until,
    until_throw
};

// ============================================================================
//  SeqAn3
// ============================================================================

template <typename rng_t, typename counter_t = size_t>
void sequential_read_impl(benchmark::State & state, rng_t && range, counter_t && counter = 0)
{
    using alphabet_t = std::ranges::range_value_t<rng_t>;
    alphabet_t a;

    for (auto _ : state)
    {
        for (auto && c : range)
            benchmark::DoNotOptimize(a = c);

        counter = 0;
    }
}

template <template <typename> typename container_t, tag id>
void sequential_read(benchmark::State & state)
{
    auto rando = seqan3::test::generate_sequence<seqan3::dna4>(10000, 0, 0);

    container_t<seqan3::dna4> seq(rando.begin(), rando.end());

    if constexpr (id == baseline)
    {
        sequential_read_impl(state, seq);
    }
    else if constexpr (id == slice)
    {
        auto long_rando = seqan3::test::generate_sequence<seqan3::dna4>(30000, 0, 0);
        container_t<seqan3::dna4> long_seq(long_rando.begin(), long_rando.end());
        sequential_read_impl(state, long_seq | seqan3::views::slice(10000, 20000));
    }
    else if constexpr (id == comp)
    {
        sequential_read_impl(state, seq | seqan3::views::complement);
    }
    else if constexpr (id == rev)
    {
        sequential_read_impl(state, seq | std::views::reverse);
    }
    else if constexpr (id == rev_comp)
    {
        sequential_read_impl(state, seq | std::views::reverse | seqan3::views::complement);
    }
    else if constexpr (id == to_char)
    {
        sequential_read_impl(state, seq | seqan3::views::to_char);
    }
    else if constexpr (id == to_upper)
    {
        auto char_rando = seqan3::test::generate_sequence<char>(10000, 0, 0);
        container_t<char> char_seq(char_rando.begin(), char_rando.end());
        sequential_read_impl(state, char_seq | seqan3::views::to_upper);
    }
    else if constexpr (id == filter)
    {
        sequential_read_impl(state, seq | std::views::filter([] (seqan3::dna4 c) { return c.to_rank() != 1; }));
    }
    else if constexpr (id == pseudofilter)
    {
        sequential_read_impl(state, seq | std::views::filter([] (seqan3::dna4) { return true; }));
    }
    else if constexpr (id == until)
    {
        auto long_rando = seqan3::test::generate_sequence<seqan3::dna4>(20000, 0, 0);
        container_t<seqan3::dna4> long_seq(long_rando.begin(), long_rando.end());
        size_t i = 0;
        sequential_read_impl(state, long_seq | seqan3::views::take_until([&i] (seqan3::dna4 c)
        {
            // on average every fourth character an 'A' so 10K processes elements
            return ( (i+= (c.to_rank() == 0)) == 2500);
        }), i);
    }
    else if constexpr (id == until_throw)
    {
        auto long_rando = seqan3::test::generate_sequence<seqan3::dna4>(20000, 0, 0);
        container_t<seqan3::dna4> long_seq(long_rando.begin(), long_rando.end());
        size_t i = 0;
        sequential_read_impl(state, long_seq | seqan3::views::take_until_or_throw([&i] (seqan3::dna4 c)
        {
            // on average every fourth character an 'A' so 10K processes elements
            return ( (i+= (c.to_rank() == 0)) == 2500);
        }), i);
    }
}

BENCHMARK_TEMPLATE(sequential_read, std::vector, tag::baseline);
BENCHMARK_TEMPLATE(sequential_read, std::vector, tag::slice);
BENCHMARK_TEMPLATE(sequential_read, std::vector, tag::comp);
BENCHMARK_TEMPLATE(sequential_read, std::vector, tag::rev);
BENCHMARK_TEMPLATE(sequential_read, std::vector, tag::rev_comp);
BENCHMARK_TEMPLATE(sequential_read, std::vector, tag::to_char);
BENCHMARK_TEMPLATE(sequential_read, std::vector, tag::to_upper);
BENCHMARK_TEMPLATE(sequential_read, std::vector, tag::filter);
BENCHMARK_TEMPLATE(sequential_read, std::vector, tag::pseudofilter);
BENCHMARK_TEMPLATE(sequential_read, std::vector, tag::until);
BENCHMARK_TEMPLATE(sequential_read, std::vector, tag::until_throw);

BENCHMARK_TEMPLATE(sequential_read, std::list, tag::baseline);
BENCHMARK_TEMPLATE(sequential_read, std::list, tag::slice);
BENCHMARK_TEMPLATE(sequential_read, std::list, tag::comp);
BENCHMARK_TEMPLATE(sequential_read, std::list, tag::rev);
BENCHMARK_TEMPLATE(sequential_read, std::list, tag::rev_comp);
BENCHMARK_TEMPLATE(sequential_read, std::list, tag::to_char);
BENCHMARK_TEMPLATE(sequential_read, std::list, tag::to_upper);
BENCHMARK_TEMPLATE(sequential_read, std::list, tag::filter);
BENCHMARK_TEMPLATE(sequential_read, std::list, tag::pseudofilter);
BENCHMARK_TEMPLATE(sequential_read, std::list, tag::until);
BENCHMARK_TEMPLATE(sequential_read, std::list, tag::until_throw);

#if SEQAN3_HAS_SEQAN2

#include <seqan/sequence.h>
#include <seqan/modifier.h>

// ============================================================================
//  SeqAn2
// ============================================================================

template <typename rng_t, typename counter_t = size_t>
void sequential_read_impl2(benchmark::State & state, rng_t && range, counter_t && counter = 0)
{
    using alphabet_t = typename seqan::Value<rng_t>::Type;
    alphabet_t a;

    for (auto _ : state)
    {
        for (auto && c : range)
            benchmark::DoNotOptimize(a = c);

        counter = 0;
    }
}

template <template <typename> typename container_t, tag id>
void sequential_read2(benchmark::State & state)
{
    auto rando = seqan3::test::generate_sequence_seqan2<seqan::Dna>(10000, 0, 0);

    container_t<seqan::Dna> seq = rando;

    if constexpr (id == baseline)
    {
        sequential_read_impl(state, seq);
    }
    else if constexpr (id == slice)
    {
        auto long_rando = seqan3::test::generate_sequence_seqan2<seqan::Dna>(30000, 0, 0);
        container_t<seqan::Dna> long_seq = long_rando;
        sequential_read_impl(state, seqan::infix(long_seq, 10000, 20000));
    }
    else if constexpr (id == comp)
    {
        using t = seqan::ModifiedString<container_t<seqan::Dna>, seqan::ModView<seqan::FunctorComplement<seqan::Dna>>>;
        sequential_read_impl(state, t{seq});
    }
    else if constexpr (id == rev)
    {
        using t = seqan::ModifiedString<container_t<seqan::Dna>, seqan::ModReverse>;
        sequential_read_impl(state, t{seq});
    }
    else if constexpr (id == rev_comp)
    {
        using t0 = seqan::ModifiedString<container_t<seqan::Dna>, seqan::ModView<seqan::FunctorComplement<seqan::Dna>>>;
        using t = seqan::ModifiedString<t0, seqan::ModReverse>;
        sequential_read_impl(state, t{seq});
    }
    else if constexpr (id == to_char)
    {
        struct char_func
        {
            using result_type = uint8_t;
            constexpr result_type operator()(seqan::Dna d) { return static_cast<result_type>(d); }
        };

        using t = seqan::ModifiedString<container_t<seqan::Dna>, seqan::ModView<char_func>>;
        sequential_read_impl(state, t{seq});
    }
    else if constexpr (id == to_upper)
    {
        auto char_rando = seqan3::test::generate_sequence_seqan2<char>(10000, 0, 0);
        container_t<char> char_seq = char_rando;
        using t = seqan::ModifiedString<container_t<seqan::Dna>, seqan::ModView<seqan::FunctorUpcase<char>>>;
        sequential_read_impl(state, t{seq});
    }
}

BENCHMARK_TEMPLATE(sequential_read2, seqan::String, tag::baseline);
BENCHMARK_TEMPLATE(sequential_read2, seqan::String, tag::slice);
BENCHMARK_TEMPLATE(sequential_read2, seqan::String, tag::comp);
BENCHMARK_TEMPLATE(sequential_read2, seqan::String, tag::rev);
BENCHMARK_TEMPLATE(sequential_read2, seqan::String, tag::rev_comp);
BENCHMARK_TEMPLATE(sequential_read2, seqan::String, tag::to_char);
BENCHMARK_TEMPLATE(sequential_read2, seqan::String, tag::to_upper);
#endif

// ============================================================================
//  run
// ============================================================================

BENCHMARK_MAIN();
