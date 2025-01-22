// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h>

#include <fstream>
#include <iostream>
#include <iterator>

#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/io/stream/detail/fast_ostreambuf_iterator.hpp>
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
    seqan3_streambuf_it_write_range,
    seqan2_stream_it,
    seqan2_stream_it_write_range
};

template <tag id>
void write_all(benchmark::State & state)
{
    /* prepare file for writing */
    seqan3::test::tmp_directory tmp{};
    auto filename = tmp.path() / "foo";
    std::ofstream os{filename, std::ios::binary};

    // sequence to write:
    std::vector<char> sequence = seqan3::test::generate_sequence<char>(10000, 0, 0);

#ifdef SEQAN3_HAS_SEQAN2
    auto seqan2_sequence = seqan3::test::generate_sequence_seqan2<char>(10000, 0, 0);
#endif
    /* start benchmark */
    for (auto _ : state)
    {
        std::ofstream os{filename, std::ios::binary};

        auto it = [&os]()
        {
            if constexpr (id == tag::std_stream_it)
            {
                return std::ostream_iterator<char>{os};
            }
            else if constexpr (id == tag::std_streambuf_it)
            {
                return std::ostreambuf_iterator<char>{os};
            }
            else if constexpr (id == tag::seqan3_streambuf_it || id == tag::seqan3_streambuf_it_write_range)
            {
                return seqan3::detail::fast_ostreambuf_iterator<char>{*os.rdbuf()};
            }
#ifdef SEQAN3_HAS_SEQAN2
            else if constexpr (id == tag::seqan2_stream_it || id == tag::seqan2_stream_it_write_range)
            {
                return seqan2::Iter<std::ofstream, seqan2::StreamIterator<seqan2::Output>>{os};
            }
#endif
        }();

        if constexpr (id == tag::seqan3_streambuf_it_write_range)
        {
            it.write_range(sequence);
        }
#ifdef SEQAN3_HAS_SEQAN2
        else if constexpr (id == tag::seqan2_stream_it_write_range)
        {
            seqan2::write(it, seqan2_sequence);
        }
#endif
        else
        {
            for (auto chr : sequence)
                *it = chr;
        }
    }
}

BENCHMARK_TEMPLATE(write_all, tag::std_stream_it);
BENCHMARK_TEMPLATE(write_all, tag::std_streambuf_it);
BENCHMARK_TEMPLATE(write_all, tag::seqan3_streambuf_it);
BENCHMARK_TEMPLATE(write_all, tag::seqan3_streambuf_it_write_range);
#ifdef SEQAN3_HAS_SEQAN2
BENCHMARK_TEMPLATE(write_all, tag::seqan2_stream_it);
BENCHMARK_TEMPLATE(write_all, tag::seqan2_stream_it_write_range);
#endif

BENCHMARK_MAIN();
