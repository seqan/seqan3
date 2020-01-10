// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <deque>
#include <list>
#include <forward_list>
#include <string>
#include <vector>

#include <benchmark/benchmark.h>

#include <seqan3/core/char_operations/predicate.hpp>
#include <seqan3/range/views/take_until.hpp>
#include <seqan3/range/views/single_pass_input.hpp>

// ============================================================================
//  sequential_read
// ============================================================================

template <typename container_t, typename adaptor_t, bool invert, bool single_pass = false, bool one_adapt = false>
void sequential_read(benchmark::State & state)
{
    container_t c;
    c.resize(1'000'000);
    uint8_t i = 0;
    for (auto & e : c)
        e = ++i; // dummy values

    uint8_t dummy = 0;

    if constexpr (std::same_as<adaptor_t, void>)
    {
        using single_t = std::conditional_t<single_pass, decltype(c | seqan3::views::single_pass_input), container_t &>;

        for (auto _ : state)
        {
            single_t s{c};
            for (auto e : s)
                if (dummy += e; e >= 101)
                    break;
        }
    }
    else
    {
        using single_t = std::conditional_t<single_pass, decltype(c | seqan3::views::single_pass_input), container_t &>;
        auto adaptor = adaptor_t{}(seqan3::is_in_interval<invert ? 0 : 101, invert ? 100 : 255>);

        for (auto _ : state)
        {
            single_t s{c};
            if constexpr (one_adapt)
            {
                auto v = s | adaptor;;
                for (auto e : v)
                    dummy += e;
            }
            else
            {
                auto v = s | adaptor | adaptor | adaptor | adaptor | adaptor
                           | adaptor | adaptor | adaptor | adaptor | adaptor;

                for (auto e : v)
                    dummy += e;

            }
        }
    }

    [[maybe_unused]] volatile uint8_t dummy2 = dummy;

    state.counters["single-pass"] = single_pass;
    state.counters["only_one_adapt"] = one_adapt;
}

// runs with chained adaptor (cannot use or_throw here)
BENCHMARK_TEMPLATE(sequential_read, std::string, void, false);
BENCHMARK_TEMPLATE(sequential_read, std::string, decltype(std::views::take_while), true);
BENCHMARK_TEMPLATE(sequential_read, std::string, decltype(seqan3::views::take_until), false);

BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, void, false);
BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, decltype(std::views::take_while), true);
BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, decltype(seqan3::views::take_until), false);

BENCHMARK_TEMPLATE(sequential_read, std::deque<uint8_t>, void, false);
BENCHMARK_TEMPLATE(sequential_read, std::deque<uint8_t>, decltype(std::views::take_while), true);
BENCHMARK_TEMPLATE(sequential_read, std::deque<uint8_t>, decltype(seqan3::views::take_until), false);

BENCHMARK_TEMPLATE(sequential_read, std::list<uint8_t>, void, false);
BENCHMARK_TEMPLATE(sequential_read, std::list<uint8_t>, decltype(std::views::take_while), true);
BENCHMARK_TEMPLATE(sequential_read, std::list<uint8_t>, decltype(seqan3::views::take_until), false);

BENCHMARK_TEMPLATE(sequential_read, std::forward_list<uint8_t>, void, false);
BENCHMARK_TEMPLATE(sequential_read, std::forward_list<uint8_t>, decltype(std::views::take_while), true);
BENCHMARK_TEMPLATE(sequential_read, std::forward_list<uint8_t>, decltype(seqan3::views::take_until), false);

BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, void, false, true);
BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, decltype(std::views::take_while), true, true);
BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, decltype(seqan3::views::take_until), false, true);

BENCHMARK_TEMPLATE(sequential_read, std::forward_list<uint8_t>, void, false, true);
BENCHMARK_TEMPLATE(sequential_read, std::forward_list<uint8_t>, decltype(std::views::take_while), true, true);
BENCHMARK_TEMPLATE(sequential_read, std::forward_list<uint8_t>, decltype(seqan3::views::take_until), false, true);

// runs with one adaptor
BENCHMARK_TEMPLATE(sequential_read, std::string, void, false, false, true);
BENCHMARK_TEMPLATE(sequential_read, std::string, decltype(std::views::take_while), true, false, true);
BENCHMARK_TEMPLATE(sequential_read, std::string, decltype(seqan3::views::take_until), false, false, true);
BENCHMARK_TEMPLATE(sequential_read, std::string, decltype(seqan3::views::take_until_or_throw), false, false, true);

BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, void, false, false, true);
BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, decltype(std::views::take_while), true, false, true);
BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, decltype(seqan3::views::take_until), false, false, true);
BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, decltype(seqan3::views::take_until_or_throw), false, false, true);

BENCHMARK_TEMPLATE(sequential_read, std::deque<uint8_t>, void, false, false, true);
BENCHMARK_TEMPLATE(sequential_read, std::deque<uint8_t>, decltype(std::views::take_while), true, false, true);
BENCHMARK_TEMPLATE(sequential_read, std::deque<uint8_t>, decltype(seqan3::views::take_until), false, false, true);
BENCHMARK_TEMPLATE(sequential_read, std::deque<uint8_t>, decltype(seqan3::views::take_until_or_throw), false, false, true);

BENCHMARK_TEMPLATE(sequential_read, std::list<uint8_t>, void, false, false, true);
BENCHMARK_TEMPLATE(sequential_read, std::list<uint8_t>, decltype(std::views::take_while), true, false, true);
BENCHMARK_TEMPLATE(sequential_read, std::list<uint8_t>, decltype(seqan3::views::take_until), false, false, true);
BENCHMARK_TEMPLATE(sequential_read, std::list<uint8_t>, decltype(seqan3::views::take_until_or_throw), false, false, true);

BENCHMARK_TEMPLATE(sequential_read, std::forward_list<uint8_t>, void, false, false, true);
BENCHMARK_TEMPLATE(sequential_read, std::forward_list<uint8_t>, decltype(std::views::take_while), true, false, true);
BENCHMARK_TEMPLATE(sequential_read, std::forward_list<uint8_t>, decltype(seqan3::views::take_until), false, false, true);
BENCHMARK_TEMPLATE(sequential_read, std::forward_list<uint8_t>, decltype(seqan3::views::take_until_or_throw), false, false, true);

BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, void, false, true, true);
BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, decltype(std::views::take_while), true, true, true);
BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, decltype(seqan3::views::take_until), false, true, true);
BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, decltype(seqan3::views::take_until_or_throw), false, true, true);

BENCHMARK_TEMPLATE(sequential_read, std::forward_list<uint8_t>, void, false, true, true);
BENCHMARK_TEMPLATE(sequential_read, std::forward_list<uint8_t>, decltype(std::views::take_while), true, true, true);
BENCHMARK_TEMPLATE(sequential_read, std::forward_list<uint8_t>, decltype(seqan3::views::take_until), false, true, true);
BENCHMARK_TEMPLATE(sequential_read, std::forward_list<uint8_t>, decltype(seqan3::views::take_until_or_throw), false, true, true);

// ============================================================================
//  run
// ============================================================================

BENCHMARK_MAIN();
