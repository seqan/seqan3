// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h>

#include <fstream>
#include <iostream>
#include <iterator>

#include <seqan3/alphabet/aminoacid/all.hpp>
#include <seqan3/io/stream/detail/fast_istreambuf_iterator.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/seqan2.hpp>
#include <seqan3/test/tmp_directory.hpp>

#if SEQAN3_HAS_SEQAN2
#    include <seqan/stream.h>
#endif

enum class tag
{
    std_stream_it,
    std_streambuf_it,
    seqan3_streambuf_it,
    seqan2_stream_it
};

template <tag id>
void read_all(benchmark::State & state)
{
    /* prepare file for reading */
    seqan3::test::tmp_directory tmp{};
    auto filename = tmp.path() / "foo";

    {
        std::ofstream os{filename, std::ios::binary};

        std::vector<seqan3::aa27> cont_rando = seqan3::test::generate_sequence<seqan3::aa27>(10000, 0, 0);

        for (size_t i = 0; i < 100; ++i)
            for (auto c : cont_rando)
                os.put(seqan3::to_char(c));
    }

    /* start benchmark */
    if constexpr (id == tag::std_stream_it)
    {
        for (auto _ : state)
        {
            char c{};
            std::ifstream s{filename, std::ios::binary};
            std::istream_iterator<char> it{s};
            std::istream_iterator<char> e{};

            for (; it != e; ++it)
                c += *it;

            benchmark::DoNotOptimize(c);
        }
    }
    else if constexpr (id == tag::std_streambuf_it)
    {
        for (auto _ : state)
        {
            char c{};
            std::ifstream s{filename, std::ios::binary};
            std::istreambuf_iterator<char> it{s};
            std::istreambuf_iterator<char> e{};

            for (; it != e; ++it)
                c += *it;

            benchmark::DoNotOptimize(c);
        }
    }
    else if constexpr (id == tag::seqan3_streambuf_it)
    {
        for (auto _ : state)
        {
            char c{};
            std::ifstream s{filename, std::ios::binary};
            seqan3::detail::fast_istreambuf_iterator<char> it{*s.rdbuf()};
            std::default_sentinel_t e{};

            for (; it != e; ++it)
                c += *it;

            benchmark::DoNotOptimize(c);
        }
    }
#ifdef SEQAN3_HAS_SEQAN2
    else if constexpr (id == tag::seqan2_stream_it)
    {
        for (auto _ : state)
        {
            char c{};
            std::ifstream s{filename, std::ios::binary};
            auto it = seqan2::Iter<std::ifstream, seqan2::StreamIterator<seqan2::Input>>{s};

            for (; !seqan2::atEnd(it); ++it)
                c += *it;

            benchmark::DoNotOptimize(c);
        }
    }
#endif
}

BENCHMARK_TEMPLATE(read_all, tag::std_stream_it);
BENCHMARK_TEMPLATE(read_all, tag::std_streambuf_it);
BENCHMARK_TEMPLATE(read_all, tag::seqan3_streambuf_it);
#ifdef SEQAN3_HAS_SEQAN2
BENCHMARK_TEMPLATE(read_all, tag::seqan2_stream_it);
#endif

BENCHMARK_MAIN();
