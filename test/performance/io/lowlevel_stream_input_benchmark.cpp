// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <fstream>
#include <iostream>
#include <iterator>

#include <seqan3/alphabet/aminoacid/all.hpp>
#include <seqan3/io/stream/iterator.hpp>
#include <seqan3/test/performance/sequence_generator.hpp>
#include <seqan3/test/seqan2.hpp>
#include <seqan3/test/tmp_filename.hpp>

#if SEQAN3_HAS_SEQAN2
#include <seqan/stream.h>
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
    seqan3::test::tmp_filename filename{"foo"};

    {
        std::ofstream os{filename.get_path(), std::ios::binary};

        std::vector<seqan3::aa27> cont_rando = seqan3::test::generate_sequence<seqan3::aa27>(10'000, 0, 0);

        for (size_t i = 0; i < 100; ++i)
            for (auto c : cont_rando)
                os.put(seqan3::to_char(c));
    }

    /* start benchmark */
    char c{}; // prevents optimisation
    if constexpr (id == tag::std_stream_it)
    {
        for (auto _ : state)
        {
            std::ifstream s{filename.get_path(), std::ios::binary};
            std::istream_iterator<char> it{s};
            std::istream_iterator<char> e{};

            for (; it != e; ++it)
                c += *it;
        }
    }
    else if constexpr (id == tag::std_streambuf_it)
    {
        for (auto _ : state)
        {
            std::ifstream s{filename.get_path(), std::ios::binary};
            std::istreambuf_iterator<char> it{s};
            std::istreambuf_iterator<char> e{};

            for (; it != e; ++it)
                c += *it;
        }
    }
    else if constexpr (id == tag::seqan3_streambuf_it)
    {
        for (auto _ : state)
        {
            std::ifstream s{filename.get_path(), std::ios::binary};
            seqan3::detail::fast_istreambuf_iterator<char> it{*s.rdbuf()};
            std::ranges::default_sentinel_t e{};

            for (; it != e; ++it)
                c += *it;
        }
    }
#ifdef SEQAN3_HAS_SEQAN2
    else if constexpr (id == tag::seqan2_stream_it)
    {
        for (auto _ : state)
        {
            std::ifstream s{filename.get_path(), std::ios::binary};
            auto it = seqan::Iter<std::ifstream, seqan::StreamIterator<seqan::Input>>{s};

            for (; !seqan::atEnd(it); ++it)
                c += *it;
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
