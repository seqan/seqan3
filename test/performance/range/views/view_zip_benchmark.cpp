// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <random>

#include <benchmark/benchmark.h>

#include <seqan3/range/views/repeat.hpp>
#include <seqan3/range/views/zip.hpp>
#include <seqan3/test/pretty_printing.hpp>

using namespace seqan3;

template <typename rng_t>
void resize_impl(rng_t &, size_t)
{}

template <typename rng_t>
    requires requires (rng_t rng, size_t size) {rng.resize(size);}
void resize_impl(rng_t & rng, size_t size)
{
    rng.resize(size);
}

template <typename ...rng_t, size_t ...N>
constexpr void resize_dispatcher(std::tuple<rng_t...> & data, std::index_sequence<N...>)
{
    (resize_impl(std::get<N>(data), 1'000'000), ...);
}

template <typename adaptor_t, typename ...rng_t, size_t ...N>
constexpr auto get_adaptor(std::tuple<rng_t...> & data, std::index_sequence<N...>)
{
    return adaptor_t{}(std::get<N>(data)...);
}

template <typename adaptor_t, typename ...rng_t>
    requires sizeof...(rng_t) > 0
void zip_factory(benchmark::State & state)
{
    std::tuple<rng_t...> data{};

    using indices = std::index_sequence_for<rng_t...>;

    resize_dispatcher(data, indices{});

    auto v = get_adaptor<adaptor_t>(data, indices{});

    for (auto _ : state)
    {
        for (auto && tup : v)
        {
            benchmark::DoNotOptimize(tup);
        }
    }
}

BENCHMARK_TEMPLATE(zip_factory, decltype(seqan3::views::zip), std::vector<size_t>);
BENCHMARK_TEMPLATE(zip_factory, decltype(seqan3::views::zip), std::vector<size_t>, std::vector<size_t>);
BENCHMARK_TEMPLATE(zip_factory, decltype(seqan3::views::zip), std::vector<size_t>, std::vector<size_t>, std::vector<size_t>);
BENCHMARK_TEMPLATE(zip_factory, decltype(seqan3::views::zip), std::vector<size_t>, std::vector<size_t>, std::vector<size_t>, std::vector<size_t>);
BENCHMARK_TEMPLATE(zip_factory, decltype(seqan3::views::zip), std::vector<size_t>, std::vector<size_t>, std::vector<char>);
BENCHMARK_TEMPLATE(zip_factory, decltype(seqan3::views::zip), std::vector<size_t>, std::vector<size_t>, std::vector<std::string>);
BENCHMARK_TEMPLATE(zip_factory, decltype(seqan3::views::zip), std::vector<size_t>, decltype(views::repeat('L')));
BENCHMARK_TEMPLATE(zip_factory, decltype(seqan3::views::zip), std::vector<size_t>, decltype(views::repeat('L')), decltype(views::repeat('L')));

BENCHMARK_TEMPLATE(zip_factory, decltype(std::views::zip), std::vector<size_t>);
BENCHMARK_TEMPLATE(zip_factory, decltype(std::views::zip), std::vector<size_t>, std::vector<size_t>);
BENCHMARK_TEMPLATE(zip_factory, decltype(std::views::zip), std::vector<size_t>, std::vector<size_t>, std::vector<size_t>);
BENCHMARK_TEMPLATE(zip_factory, decltype(std::views::zip), std::vector<size_t>, std::vector<size_t>, std::vector<size_t>, std::vector<size_t>);
BENCHMARK_TEMPLATE(zip_factory, decltype(std::views::zip), std::vector<size_t>, std::vector<size_t>, std::vector<char>);
BENCHMARK_TEMPLATE(zip_factory, decltype(std::views::zip), std::vector<size_t>, std::vector<size_t>, std::vector<std::string>);
BENCHMARK_TEMPLATE(zip_factory, decltype(std::views::zip), std::vector<size_t>, decltype(views::repeat('L')));
BENCHMARK_TEMPLATE(zip_factory, decltype(std::views::zip), std::vector<size_t>, decltype(views::repeat('L')), decltype(views::repeat('L')));

// ============================================================================
//  run
// ============================================================================

BENCHMARK_MAIN();
