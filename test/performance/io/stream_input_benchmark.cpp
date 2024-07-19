// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h>

#include <iostream>
#include <memory>
#include <sstream>

#include <seqan3/io/stream/detail/fast_istreambuf_iterator.hpp>

#if defined(SEQAN3_HAS_ZLIB)
#    include <seqan3/contrib/stream/bgzf_istream.hpp>
#    include <seqan3/contrib/stream/bgzf_ostream.hpp>
#    include <seqan3/contrib/stream/gz_istream.hpp>
#    include <seqan3/contrib/stream/gz_ostream.hpp>
#endif

// only benchmark BZIP2 if explicitly requested, because slow setup
#if !defined(SEQAN3_BENCH_BZIP2) && defined(SEQAN3_HAS_BZIP2)
#    undef SEQAN3_HAS_BZIP2
#endif

#if defined(SEQAN3_HAS_BZIP2)
#    include <seqan3/contrib/stream/bz2_istream.hpp>
#    include <seqan3/contrib/stream/bz2_ostream.hpp>
#endif

// SEQAN2
#if __has_include(<seqan/stream.h>)
#    define SEQAN3_HAS_SEQAN2 1

#    if defined(SEQAN3_HAS_ZLIB)
#        define SEQAN_HAS_ZLIB 1
#    endif

#    if defined(SEQAN3_HAS_BZIP2)
#        define SEQAN_HAS_BZIP2 1
#    endif

#    include <seqan/stream.h>
#endif

#ifndef NDEBUG
inline constexpr size_t input_size = 10000;
#else
inline constexpr size_t input_size = 10'000'000;
#endif // NDEBUG

std::string input{[]()
                  {
                      std::string line{"The quick brown fox jumps over the lazy dog"};
                      std::string ret;
                      for (size_t i = 0; i < input_size; ++i)
                          ret += line;
                      return ret;
                  }()};

template <typename t>
std::string const input_comp;

#ifdef SEQAN3_HAS_SEQAN2
template <>
std::string const & input_comp<seqan2::Nothing> = input;
#endif

#if defined(SEQAN3_HAS_ZLIB)
template <>
std::string const input_comp<seqan3::contrib::gz_istream>{
    []()
    {
        std::ostringstream ret;
        { // In scope to force flush of ostream on destruction.
            seqan3::contrib::gz_ostream os{ret};
            std::copy(input.begin(), input.end(), std::ostreambuf_iterator<char>(os));
        }
        return ret.str();
    }()};

template <>
std::string const input_comp<seqan3::contrib::bgzf_istream>{
    []()
    {
        std::ostringstream ret;
        { // In scope to force flush of ostream on destruction.
            seqan3::contrib::bgzf_ostream os{ret};
            std::copy(input.begin(), input.end(), std::ostreambuf_iterator<char>(os));
        }
        return ret.str();
    }()};
#    ifdef SEQAN3_HAS_SEQAN2
template <>
std::string const & input_comp<seqan2::GZFile> = input_comp<seqan3::contrib::gz_istream>;

template <>
std::string const & input_comp<seqan2::BgzfFile> = input_comp<seqan3::contrib::bgzf_istream>;
#    endif
#endif

#if defined(SEQAN3_HAS_BZIP2)
template <>
std::string const input_comp<seqan3::contrib::bz2_istream>{
    []()
    {
        std::ostringstream ret;
        { // In scope to force flush of ostream on destruction.
            seqan3::contrib::bz2_ostream os{ret};
            std::copy(input.begin(), input.end(), std::ostreambuf_iterator<char>(os));
        }
        return ret.str();
    }()};
#    ifdef SEQAN3_HAS_SEQAN2
template <>
std::string const & input_comp<seqan2::BZ2File> = input_comp<seqan3::contrib::bz2_istream>;
#    endif
#endif

// ============================================================================
//  plain benchmark of ostringstream
// ============================================================================

void uncompressed(benchmark::State & state)
{
    std::istringstream s{input};
    size_t i = 0;
    for (auto _ : state)
    {
        s.clear();
        s.seekg(0, std::ios::beg);
        seqan3::detail::fast_istreambuf_iterator<char> it{*s.rdbuf()};

        for (; it != std::default_sentinel; ++it)
            i += *it;
    }

    state.counters["iterations_per_run"] = i;
}

BENCHMARK(uncompressed);

// ============================================================================
//  compression applied
// ============================================================================

template <typename compressed_istream_t>
void compressed(benchmark::State & state)
{
    std::istringstream s{input_comp<compressed_istream_t>};

    size_t i = 0;
    for (auto _ : state)
    {
        s.clear();
        s.seekg(0, std::ios::beg);
        compressed_istream_t comp{s};
        seqan3::detail::fast_istreambuf_iterator<char> it{*comp.rdbuf()};

        for (size_t v = 0; v < input_comp<compressed_istream_t>.size(); ++v)
        {
            i += *it;
            ++it;
        }
    }

    state.counters["iterations_per_run"] = i;
}

#if defined(SEQAN3_HAS_ZLIB)
BENCHMARK_TEMPLATE(compressed, seqan3::contrib::gz_istream);
BENCHMARK_TEMPLATE(compressed, seqan3::contrib::bgzf_istream);
#endif

#if defined(SEQAN3_HAS_BZIP2)
BENCHMARK_TEMPLATE(compressed, seqan3::contrib::bz2_istream);
#endif

// ============================================================================
//  compression applied, but stuffed into plain istream
// ============================================================================

template <typename compressed_istream_t>
void compressed_type_erased(benchmark::State & state)
{
    std::istringstream s{input_comp<compressed_istream_t>};
    size_t i = 0;
    for (auto _ : state)
    {
        s.clear();
        s.seekg(0, std::ios::beg);
        std::unique_ptr<std::istream> comp{new compressed_istream_t{s}};
        seqan3::detail::fast_istreambuf_iterator<char> it{*comp->rdbuf()};

        for (size_t v = 0; v < input_comp<compressed_istream_t>.size(); ++v)
        {
            i += *it;
            ++it;
        }
    }

    state.counters["iterations_per_run"] = i;
}

#if defined(SEQAN3_HAS_ZLIB)
BENCHMARK_TEMPLATE(compressed_type_erased, seqan3::contrib::gz_istream);
BENCHMARK_TEMPLATE(compressed_type_erased, seqan3::contrib::bgzf_istream);
#endif
#if defined(SEQAN3_HAS_BZIP2)
BENCHMARK_TEMPLATE(compressed_type_erased, seqan3::contrib::bz2_istream);
#endif

// ============================================================================
//  compression applied, but stuffed into plain istream, also stringstream erased
// ============================================================================

template <typename compressed_istream_t>
void compressed_type_erased2(benchmark::State & state)
{
    std::unique_ptr<std::istream> s{new std::istringstream{input_comp<compressed_istream_t>}};
    size_t i = 0;
    for (auto _ : state)
    {
        s->clear();
        s->seekg(0, std::ios::beg);
        std::unique_ptr<std::istream> comp{new compressed_istream_t{*s}};
        seqan3::detail::fast_istreambuf_iterator<char> it{*comp->rdbuf()};

        for (size_t v = 0; v < input_comp<compressed_istream_t>.size(); ++v)
        {
            i += *it;
            ++it;
        }
    }

    state.counters["iterations_per_run"] = i;
}

#if defined(SEQAN3_HAS_ZLIB)
BENCHMARK_TEMPLATE(compressed_type_erased2, seqan3::contrib::gz_istream);
BENCHMARK_TEMPLATE(compressed_type_erased2, seqan3::contrib::bgzf_istream);
#endif
#if defined(SEQAN3_HAS_BZIP2)
BENCHMARK_TEMPLATE(compressed_type_erased2, seqan3::contrib::bz2_istream);
#endif

// ============================================================================
//  seqan2 virtual stream
// ============================================================================

#ifdef SEQAN3_HAS_SEQAN2

void seqan2_uncompressed(benchmark::State & state)
{
    std::istringstream s{input};
    size_t i = 0;
    for (auto _ : state)
    {
        s.clear();
        s.seekg(0, std::ios::beg);
        auto it = seqan2::Iter<std::istringstream, seqan2::StreamIterator<seqan2::Input>>(s);

        for (size_t v = 0; v < input.size(); ++v)
        {
            i += *it;
            ++it;
        }
    }

    state.counters["iterations_per_run"] = i;
}

template <typename compression_t, typename stream_t>
void seqan2_compressed_impl(benchmark::State & state)
{
    std::istringstream s{input_comp<compression_t>};
    size_t i = 0;
    for (auto _ : state)
    {
        s.clear();
        s.seekg(0, std::ios::beg);
        stream_t comp{s};
        auto it = seqan2::Iter<std::istringstream, seqan2::StreamIterator<seqan2::Input>>(comp);

        for (size_t v = 0; v < input_comp<compression_t>.size(); ++v)
        {
            i += *it;
            ++it;
        }
    }

    state.counters["iterations_per_run"] = i;
}

template <typename compression_type>
void seqan2_compressed(benchmark::State & state)
{
#    ifdef SEQAN_HAS_ZLIB
    if constexpr (std::is_same_v<compression_type, seqan2::GZFile>)
        seqan2_compressed_impl<seqan2::GZFile, zlib_stream::zip_istream>(state);
    else if constexpr (std::is_same_v<compression_type, seqan2::BgzfFile>)
        seqan2_compressed_impl<seqan2::BgzfFile, seqan2::bgzf_istream>(state);
#    endif // SEQAN_HAS_ZLIB

#    ifdef SEQAN_HAS_BZIP2
    if constexpr (std::is_same_v<compression_type, seqan2::BZ2File>)
        seqan2_compressed_impl<seqan2::BZ2File, bzip2_stream::bzip2_istream>(state);
#    endif // SEQAN_HAS_BZIP2

    if constexpr (std::is_same_v<compression_type, seqan2::Nothing>)
        seqan2_uncompressed(state);
}

BENCHMARK_TEMPLATE(seqan2_compressed, seqan2::Nothing);

#    ifdef SEQAN_HAS_ZLIB
BENCHMARK_TEMPLATE(seqan2_compressed, seqan2::GZFile);
BENCHMARK_TEMPLATE(seqan2_compressed, seqan2::BgzfFile);
#    endif
#    ifdef SEQAN_HAS_BZIP2
BENCHMARK_TEMPLATE(seqan2_compressed, seqan2::BZ2File);
#    endif

#endif // SEQAN3_HAS_SEQAN2

// ============================================================================
//  instantiate tests
// ============================================================================

BENCHMARK_MAIN();
