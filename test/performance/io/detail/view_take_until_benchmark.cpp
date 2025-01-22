// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h>

#include <deque>
#include <forward_list>
#include <list>
#include <string>
#include <vector>

#include <seqan3/io/views/detail/take_until_view.hpp>
#include <seqan3/utility/char_operations/predicate.hpp>
#include <seqan3/utility/views/single_pass_input.hpp>

// ============================================================================
//  sequential_read
// ============================================================================

template <typename container_t,
          typename adaptor_t,
          bool invert_predicate,
          bool use_single_pass = false,
          bool use_multiple_adaptors = true>
void sequential_read(benchmark::State & state)
{
    container_t const container = []()
    {
        auto values = std::views::iota(0u, 1000u);
        return container_t{values.begin(), values.end()};
    }();
    size_t sum{};

    // We either use `seqan3::views::single_pass_input` on the container or access it via const lvalue reference.
    using single_pass_or_ref_t = std::
        conditional_t<use_single_pass, decltype(container | seqan3::views::single_pass_input), container_t const &>;

    if constexpr (std::same_as<adaptor_t, void>) // No adaptor
    {
        for (auto _ : state)
        {
            single_pass_or_ref_t single_pass_or_ref{container};
            for (auto elem : single_pass_or_ref)
                if (sum += elem; elem >= 101)
                    break;
        }
    }
    else // {seqan3,std}::views::take* adaptor
    {
        constexpr auto predicate = seqan3::is_in_interval < invert_predicate ? 0 : 101, invert_predicate ? 100 : 255 > ;
        auto adaptor = adaptor_t{}(predicate);

        for (auto _ : state)
        {
            single_pass_or_ref_t single_pass_or_ref{container};
            if constexpr (use_multiple_adaptors)
            {
                auto view = single_pass_or_ref | adaptor | adaptor | adaptor | adaptor;
                for (auto elem : view)
                    sum += elem;
            }
            else
            {
                auto view = single_pass_or_ref | adaptor;
                for (auto elem : view)
                    sum += elem;
            }
        }
    }

    benchmark::DoNotOptimize(sum);

    state.counters["use-single-pass"] = use_single_pass;
    state.counters["use-multiple-adaptors"] = use_multiple_adaptors;
}

// runs with chained adaptor (cannot use or_throw here)
BENCHMARK_TEMPLATE(sequential_read, std::string, void, false);
BENCHMARK_TEMPLATE(sequential_read, std::string, decltype(std::views::take_while), true);
BENCHMARK_TEMPLATE(sequential_read, std::string, decltype(seqan3::detail::take_until), false);

BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, void, false);
BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, decltype(std::views::take_while), true);
BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, decltype(seqan3::detail::take_until), false);

BENCHMARK_TEMPLATE(sequential_read, std::deque<uint8_t>, void, false);
BENCHMARK_TEMPLATE(sequential_read, std::deque<uint8_t>, decltype(std::views::take_while), true);
BENCHMARK_TEMPLATE(sequential_read, std::deque<uint8_t>, decltype(seqan3::detail::take_until), false);

BENCHMARK_TEMPLATE(sequential_read, std::list<uint8_t>, void, false);
BENCHMARK_TEMPLATE(sequential_read, std::list<uint8_t>, decltype(std::views::take_while), true);
BENCHMARK_TEMPLATE(sequential_read, std::list<uint8_t>, decltype(seqan3::detail::take_until), false);

BENCHMARK_TEMPLATE(sequential_read, std::forward_list<uint8_t>, void, false);
BENCHMARK_TEMPLATE(sequential_read, std::forward_list<uint8_t>, decltype(std::views::take_while), true);
BENCHMARK_TEMPLATE(sequential_read, std::forward_list<uint8_t>, decltype(seqan3::detail::take_until), false);

BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, void, false, true);
BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, decltype(std::views::take_while), true, true);
BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, decltype(seqan3::detail::take_until), false, true);

BENCHMARK_TEMPLATE(sequential_read, std::forward_list<uint8_t>, void, false, true);
BENCHMARK_TEMPLATE(sequential_read, std::forward_list<uint8_t>, decltype(std::views::take_while), true, true);
BENCHMARK_TEMPLATE(sequential_read, std::forward_list<uint8_t>, decltype(seqan3::detail::take_until), false, true);

// runs with one adaptor
BENCHMARK_TEMPLATE(sequential_read, std::string, void, false, false, false);
BENCHMARK_TEMPLATE(sequential_read, std::string, decltype(std::views::take_while), true, false, false);
BENCHMARK_TEMPLATE(sequential_read, std::string, decltype(seqan3::detail::take_until), false, false, false);
BENCHMARK_TEMPLATE(sequential_read, std::string, decltype(seqan3::detail::take_until_or_throw), false, false, false);

BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, void, false, false, false);
BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, decltype(std::views::take_while), true, false, false);
BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, decltype(seqan3::detail::take_until), false, false, false);
BENCHMARK_TEMPLATE(sequential_read,
                   std::vector<uint8_t>,
                   decltype(seqan3::detail::take_until_or_throw),
                   false,
                   false,
                   false);

BENCHMARK_TEMPLATE(sequential_read, std::deque<uint8_t>, void, false, false, false);
BENCHMARK_TEMPLATE(sequential_read, std::deque<uint8_t>, decltype(std::views::take_while), true, false, false);
BENCHMARK_TEMPLATE(sequential_read, std::deque<uint8_t>, decltype(seqan3::detail::take_until), false, false, false);
BENCHMARK_TEMPLATE(sequential_read,
                   std::deque<uint8_t>,
                   decltype(seqan3::detail::take_until_or_throw),
                   false,
                   false,
                   false);

BENCHMARK_TEMPLATE(sequential_read, std::list<uint8_t>, void, false, false, false);
BENCHMARK_TEMPLATE(sequential_read, std::list<uint8_t>, decltype(std::views::take_while), true, false, false);
BENCHMARK_TEMPLATE(sequential_read, std::list<uint8_t>, decltype(seqan3::detail::take_until), false, false, false);
BENCHMARK_TEMPLATE(sequential_read,
                   std::list<uint8_t>,
                   decltype(seqan3::detail::take_until_or_throw),
                   false,
                   false,
                   false);

BENCHMARK_TEMPLATE(sequential_read, std::forward_list<uint8_t>, void, false, false, false);
BENCHMARK_TEMPLATE(sequential_read, std::forward_list<uint8_t>, decltype(std::views::take_while), true, false, false);
BENCHMARK_TEMPLATE(sequential_read,
                   std::forward_list<uint8_t>,
                   decltype(seqan3::detail::take_until),
                   false,
                   false,
                   false);
BENCHMARK_TEMPLATE(sequential_read,
                   std::forward_list<uint8_t>,
                   decltype(seqan3::detail::take_until_or_throw),
                   false,
                   false,
                   false);

BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, void, false, true, false);
BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, decltype(std::views::take_while), true, true, false);
BENCHMARK_TEMPLATE(sequential_read, std::vector<uint8_t>, decltype(seqan3::detail::take_until), false, true, false);
BENCHMARK_TEMPLATE(sequential_read,
                   std::vector<uint8_t>,
                   decltype(seqan3::detail::take_until_or_throw),
                   false,
                   true,
                   false);

BENCHMARK_TEMPLATE(sequential_read, std::forward_list<uint8_t>, void, false, true, false);
BENCHMARK_TEMPLATE(sequential_read, std::forward_list<uint8_t>, decltype(std::views::take_while), true, true, false);
BENCHMARK_TEMPLATE(sequential_read,
                   std::forward_list<uint8_t>,
                   decltype(seqan3::detail::take_until),
                   false,
                   true,
                   false);
BENCHMARK_TEMPLATE(sequential_read,
                   std::forward_list<uint8_t>,
                   decltype(seqan3::detail::take_until_or_throw),
                   false,
                   true,
                   false);

// ============================================================================
//  run
// ============================================================================

BENCHMARK_MAIN();
