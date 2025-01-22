// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <benchmark/benchmark.h>

#include <iostream>
#include <memory>
#include <sstream>

#if defined(SEQAN3_HAS_ZLIB)
#    include <seqan3/contrib/stream/bgzf_ostream.hpp>
#    include <seqan3/contrib/stream/gz_ostream.hpp>
#endif

#if defined(SEQAN3_HAS_BZIP2)
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

// ============================================================================
//  plain benchmark of ostringstream
// ============================================================================

void uncompressed(benchmark::State & state)
{
    std::ostringstream os;

    std::ostreambuf_iterator<char> oit{os};

    size_t i = 0;
    for (auto _ : state)
        oit = static_cast<char>(i++ % 128);
}
BENCHMARK(uncompressed);

// ============================================================================
//  compression applied
// ============================================================================

template <typename compressed_ostream_t>
void compressed(benchmark::State & state)
{
    std::ostringstream os;

    compressed_ostream_t ogzf{os};

    std::ostreambuf_iterator<char> oit{ogzf};

    size_t i = 0;
    for (auto _ : state)
        oit = static_cast<char>(i++ % 128);
}

#if defined(SEQAN3_HAS_ZLIB)
BENCHMARK_TEMPLATE(compressed, seqan3::contrib::gz_ostream);
BENCHMARK_TEMPLATE(compressed, seqan3::contrib::bgzf_ostream);
#endif
#if defined(SEQAN3_HAS_BZIP2)
BENCHMARK_TEMPLATE(compressed, seqan3::contrib::bz2_ostream);
#endif

// ============================================================================
//  compression applied, but stuffed into plain ostream
// ============================================================================

template <typename compressed_ostream_t>
void compressed_type_erased(benchmark::State & state)
{
    std::ostringstream os;

    std::unique_ptr<std::ostream> ogzf{new compressed_ostream_t{os}};

    std::ostreambuf_iterator<char> oit{*ogzf};

    size_t i = 0;
    for (auto _ : state)
        oit = static_cast<char>(i++ % 128);
}

#if defined(SEQAN3_HAS_ZLIB)
BENCHMARK_TEMPLATE(compressed_type_erased, seqan3::contrib::gz_ostream);
BENCHMARK_TEMPLATE(compressed_type_erased, seqan3::contrib::bgzf_ostream);
#endif
#if defined(SEQAN3_HAS_BZIP2)
BENCHMARK_TEMPLATE(compressed_type_erased, seqan3::contrib::bz2_ostream);
#endif

// ============================================================================
//  compression applied, but stuffed into plain ostream, also stringstream erased
// ============================================================================

template <typename compressed_ostream_t>
void compressed_type_erased2(benchmark::State & state)
{
    std::unique_ptr<std::ostream> os{new std::ostringstream{}};

    std::unique_ptr<std::ostream> ogzf{new compressed_ostream_t{*os}};

    std::ostreambuf_iterator<char> oit{*ogzf};

    size_t i = 0;
    for (auto _ : state)
        oit = static_cast<char>(i++ % 128);
}

#if defined(SEQAN3_HAS_ZLIB)
BENCHMARK_TEMPLATE(compressed_type_erased2, seqan3::contrib::gz_ostream);
BENCHMARK_TEMPLATE(compressed_type_erased2, seqan3::contrib::bgzf_ostream);
#endif
#if defined(SEQAN3_HAS_BZIP2)
BENCHMARK_TEMPLATE(compressed_type_erased2, seqan3::contrib::bz2_ostream);
#endif

// ============================================================================
//  seqan2 virtual stream
// ============================================================================

#ifdef SEQAN3_HAS_SEQAN2
template <typename compression_type>
void seqan2_compressed(benchmark::State & state)
{
    std::ostringstream os;

    seqan2::VirtualStream<char, seqan2::Output> ogzf;
    compression_type tag;
    open(ogzf, os, tag);

    size_t i = 0;
    for (auto _ : state)
        write(ogzf, static_cast<char>(i++ % 128));
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
