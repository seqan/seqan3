// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#include <benchmark/benchmark.h>

#include <iostream>
#include <memory>
#include <sstream>

#ifdef SEQAN3_HAS_ZLIB
    #include <seqan3/contrib/stream/bgzf_istream.hpp>
    #include <seqan3/contrib/stream/bgzf_ostream.hpp>
    #include <seqan3/contrib/stream/gz_istream.hpp>
    #include <seqan3/contrib/stream/gz_ostream.hpp>
#endif

// only benchmark BZIP2 if explicitly requested, because slow setup
#if !defined(SEQAN3_BENCH_BZIP2) && defined(SEQAN3_HAS_BZIP2)
    #undef SEQAN3_HAS_BZIP2
#endif

#ifdef SEQAN3_HAS_BZIP2
    #include <seqan3/contrib/stream/bz2_istream.hpp>
    #include <seqan3/contrib/stream/bz2_ostream.hpp>
#endif

// SEQAN2
#if __has_include(<seqan/stream.h>)
    #define SEQAN3_HAS_SEQAN2 1

    #ifdef SEQAN3_HAS_ZLIB
        #define SEQAN_HAS_ZLIB 1
    #endif

    #ifdef SEQAN3_HAS_BZIP2
        #define SEQAN_HAS_BZIP2 1
    #endif

    #include <seqan/stream.h>
#endif

using namespace seqan3;

std::string input
{
    [] ()
    {
        std::string line{"The quick brown fox jumps over the lazy dog"};
        std::string ret;
        for (size_t i = 0; i < 10'000'000; ++i)
            ret += line;
        return ret;
    } ()
};

template <typename t>
std::string const input_comp;

#ifdef SEQAN3_HAS_SEQAN2
template <>
std::string const & input_comp<seqan::Nothing> = input;
#endif

#ifdef SEQAN3_HAS_ZLIB
template <>
std::string const input_comp<contrib::gz_istream>
{
    [] ()
    {
        std::ostringstream ret;
        contrib::gz_ostream os{ret};
        std::copy(input.begin(), input.end(), std::ostreambuf_iterator<char>(os));
        return ret.str();
    } ()
};

template <>
std::string const input_comp<contrib::bgzf_istream>
{
    [] ()
    {
        std::ostringstream ret;
        contrib::bgzf_ostream os{ret};
        std::copy(input.begin(), input.end(), std::ostreambuf_iterator<char>(os));
        return ret.str();
    } ()
};
#ifdef SEQAN3_HAS_SEQAN2
template <>
std::string const & input_comp<seqan::GZFile> = input_comp<contrib::gz_istream>;
#endif
#endif

#ifdef SEQAN3_HAS_BZIP2
template <>
std::string const input_comp<contrib::bz2_istream>
{
    [] ()
    {
        std::ostringstream ret;
        contrib::bz2_ostream os{ret};
        std::copy(input.begin(), input.end(), std::ostreambuf_iterator<char>(os));
        return ret.str();
    } ()
};
#ifdef SEQAN3_HAS_SEQAN2
template <>
std::string const & input_comp<seqan::BZ2File> = input_comp<contrib::bz2_istream>;
#endif
#endif

// ============================================================================
//  plain benchmark of ostringstream
// ============================================================================

void uncompressed(benchmark::State & state)
{
    std::istringstream s{input};
    std::istreambuf_iterator<char> it{s};
    size_t i = 0;
    for (auto _ : state)
    {
        i += *it;
        ++it;
    }

    [[maybe_unused]] volatile size_t j = i;
}

BENCHMARK(uncompressed);

// ============================================================================
//  compression applied
// ============================================================================

template <typename compressed_istream_t>
void compressed(benchmark::State & state)
{
    std::istringstream s{input_comp<compressed_istream_t>};
    compressed_istream_t comp{s};
    std::istreambuf_iterator<char> it{comp};

    size_t i = 0;
    for (auto _ : state)
    {
        i += *it;
        ++it;
    }

    [[maybe_unused]] volatile size_t j = i;
}

#ifdef SEQAN3_HAS_ZLIB
BENCHMARK_TEMPLATE(compressed, contrib::gz_istream);
BENCHMARK_TEMPLATE(compressed, contrib::bgzf_istream);
#endif

#ifdef SEQAN3_HAS_BZIP2
BENCHMARK_TEMPLATE(compressed, contrib::bz2_istream);
#endif

// ============================================================================
//  compression applied, but stuffed into plain istream
// ============================================================================

template <typename compressed_istream_t>
void compressed_type_erased(benchmark::State & state)
{
    std::istringstream s{input_comp<compressed_istream_t>};
    std::unique_ptr<std::istream> comp{new compressed_istream_t{s}};
    std::istreambuf_iterator<char> it{*comp};

    size_t i = 0;
    for (auto _ : state)
    {
        i += *it;
        ++it;
    }

    [[maybe_unused]] volatile size_t j = i;
}

#ifdef SEQAN3_HAS_ZLIB
BENCHMARK_TEMPLATE(compressed_type_erased, contrib::gz_istream);
BENCHMARK_TEMPLATE(compressed_type_erased, contrib::bgzf_istream);
#endif
#ifdef SEQAN3_HAS_BZIP2
BENCHMARK_TEMPLATE(compressed_type_erased, contrib::bz2_istream);
#endif

// ============================================================================
//  compression applied, but stuffed into plain istream, also stringstream erased
// ============================================================================

template <typename compressed_istream_t>
void compressed_type_erased2(benchmark::State & state)
{
    std::unique_ptr<std::istream> s{new std::istringstream{input_comp<compressed_istream_t>}};
    std::unique_ptr<std::istream> comp{new compressed_istream_t{*s}};
    std::istreambuf_iterator<char> it{*comp};

    size_t i = 0;
    for (auto _ : state)
    {
        i += *it;
        ++it;
    }

    [[maybe_unused]] volatile size_t j = i;
}

#ifdef SEQAN3_HAS_ZLIB
BENCHMARK_TEMPLATE(compressed_type_erased2, contrib::gz_istream);
BENCHMARK_TEMPLATE(compressed_type_erased2, contrib::bgzf_istream);
#endif
#ifdef SEQAN3_HAS_BZIP2
BENCHMARK_TEMPLATE(compressed_type_erased2, contrib::bz2_istream);
#endif

// ============================================================================
//  seqan2 virtual stream
// ============================================================================

#ifdef SEQAN3_HAS_SEQAN2
template <typename compression_type>
void seqan2_compressed(benchmark::State & state)
{
    std::istringstream s{input_comp<compression_type>};
    seqan::VirtualStream<char, seqan::Input> comp;
    compression_type tag;
    open(comp, s, tag);

    auto it = seqan::directionIterator(comp, seqan::Input());

    size_t i = 0;
    for (auto _ : state)
    {
        i += *it;
        ++it;
    }

    [[maybe_unused]] volatile size_t j = i;
}
BENCHMARK_TEMPLATE(seqan2_compressed, seqan::Nothing);

#ifdef SEQAN_HAS_ZLIB
BENCHMARK_TEMPLATE(seqan2_compressed, seqan::GZFile);
BENCHMARK_TEMPLATE(seqan2_compressed, seqan::BgzfFile);
#endif
#ifdef SEQAN_HAS_BZIP2
BENCHMARK_TEMPLATE(seqan2_compressed, seqan::BZ2File);
#endif

#endif // SEQAN3_HAS_SEQAN2

// ============================================================================
//  instantiate tests
// ============================================================================

BENCHMARK_MAIN();
