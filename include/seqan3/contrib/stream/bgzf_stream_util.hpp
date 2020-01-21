// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides stream compression utilities.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <array>
#include <cstring>
#include <memory>
#include <thread>

#if SEQAN3_HAS_ZLIB
// Zlib headers
#include <zlib.h>
#else
#error "This file cannot be used when building without GZip-support."
#endif  // SEQAN3_HAS_ZLIB

#include <seqan3/core/bit_manipulation.hpp>
#include <seqan3/core/type_traits/range.hpp>
#include <seqan3/io/detail/magic_header.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/span>

namespace seqan3::contrib
{

/*!\brief A static variable indicating the number of threads to use for the bgzf-streams.
 *       Defaults to std::thread::hardware_concurrency.
 */
inline static uint64_t bgzf_thread_count = std::thread::hardware_concurrency();

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Classes
// ============================================================================

// Special end-of-file marker defined by the BGZF compression format.
// See: https://samtools.github.io/hts-specs/SAMv1.pdf
static constexpr std::array<char, 28> BGZF_END_OF_FILE_MARKER {{'\x1f', '\x8b', '\x08', '\x04',
                                                                '\x00', '\x00', '\x00', '\x00',
                                                                '\x00', '\xff', '\x06', '\x00',
                                                                '\x42', '\x43', '\x02', '\x00',
                                                                '\x1b', '\x00', '\x03', '\x00',
                                                                '\x00', '\x00', '\x00', '\x00',
                                                                '\x00', '\x00', '\x00', '\x00'}};

template <typename TAlgTag>
struct CompressionContext {};

template <typename TAlgTag>
struct DefaultPageSize;

template <>
struct CompressionContext<detail::gz_compression>
{
    z_stream strm;

    CompressionContext()
    {
        std::memset(&strm, 0, sizeof(z_stream));
    }
};

template <>
struct CompressionContext<detail::bgzf_compression>:
    CompressionContext<detail::gz_compression>
{
    static constexpr size_t BLOCK_HEADER_LENGTH = detail::bgzf_compression::magic_header.size();
    unsigned char headerPos;
};

template <>
struct DefaultPageSize<detail::bgzf_compression>
{
    static const unsigned MAX_BLOCK_SIZE = 64 * 1024;
    static const unsigned BLOCK_FOOTER_LENGTH = 8;
    // 5 bytes block overhead (see 3.2.4. at https://tools.ietf.org/html/rfc1951)
    static const unsigned ZLIB_BLOCK_OVERHEAD = 5;

    // Reduce the maximal input size, such that the compressed data
    // always fits in one block even for level Z_NO_COMPRESSION.
    enum { BLOCK_HEADER_LENGTH = CompressionContext<detail::bgzf_compression>::BLOCK_HEADER_LENGTH };
    static const unsigned VALUE = MAX_BLOCK_SIZE - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH - ZLIB_BLOCK_OVERHEAD;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function compressInit()
// ----------------------------------------------------------------------------

inline void
compressInit(CompressionContext<detail::gz_compression> & ctx)
{
    const int GZIP_WINDOW_BITS = -15;   // no zlib header
    const int Z_DEFAULT_MEM_LEVEL = 8;

    ctx.strm.zalloc = NULL;
    ctx.strm.zfree = NULL;

    // (weese:) We use Z_BEST_SPEED instead of Z_DEFAULT_COMPRESSION as it turned out
    //          to be 2x faster and produces only 7% bigger output
//    int status = deflateInit2(&ctx.strm, Z_DEFAULT_COMPRESSION, Z_DEFLATED,
//                              GZIP_WINDOW_BITS, Z_DEFAULT_MEM_LEVEL, Z_DEFAULT_STRATEGY);
    int status = deflateInit2(&ctx.strm, Z_BEST_SPEED, Z_DEFLATED,
                              GZIP_WINDOW_BITS, Z_DEFAULT_MEM_LEVEL, Z_DEFAULT_STRATEGY);
    if (status != Z_OK)
        throw io_error("Calling deflateInit2() failed for gz file.");
}

// ----------------------------------------------------------------------------
// Function compressInit()
// ----------------------------------------------------------------------------

inline void
compressInit(CompressionContext<detail::bgzf_compression> & ctx)
{
    compressInit(static_cast<CompressionContext<detail::gz_compression> &>(ctx));
    ctx.headerPos = 0;
}

// ----------------------------------------------------------------------------
// Helper Function _bgzfUnpackXX()
// ----------------------------------------------------------------------------

inline uint16_t
_bgzfUnpack16(char const * buffer)
{
    uint16_t tmp;
    std::uninitialized_copy(buffer, buffer + sizeof(uint16_t), reinterpret_cast<char *>(&tmp));
    return detail::to_little_endian(tmp);
}

inline uint32_t
_bgzfUnpack32(char const * buffer)
{
    uint32_t tmp;
    std::uninitialized_copy(buffer, buffer + sizeof(uint32_t), reinterpret_cast<char *>(&tmp));
    return detail::to_little_endian(tmp);
}

// ----------------------------------------------------------------------------
// Helper Function _bgzfPackXX()
// ----------------------------------------------------------------------------

inline void
_bgzfPack16(char * buffer, uint16_t value)
{
    value = detail::to_little_endian(value);
    std::uninitialized_copy(reinterpret_cast<char *>(&value),
                            reinterpret_cast<char *>(&value) + sizeof(uint16_t),
                            buffer);
}

inline void
_bgzfPack32(char * buffer, uint32_t value)
{
    value = detail::to_little_endian(value);
    std::uninitialized_copy(reinterpret_cast<char *>(&value),
                            reinterpret_cast<char *>(&value) + sizeof(uint32_t),
                            buffer);
}

// ----------------------------------------------------------------------------
// Function _compressBlock()
// ----------------------------------------------------------------------------

template <typename TDestValue, typename TDestCapacity, typename TSourceValue, typename TSourceLength>
inline TDestCapacity
_compressBlock(TDestValue *dstBegin,   TDestCapacity dstCapacity,
               TSourceValue *srcBegin, TSourceLength srcLength, CompressionContext<detail::bgzf_compression> & ctx)
{
    const size_t BLOCK_HEADER_LENGTH = DefaultPageSize<detail::bgzf_compression>::BLOCK_HEADER_LENGTH;
    const size_t BLOCK_FOOTER_LENGTH = DefaultPageSize<detail::bgzf_compression>::BLOCK_FOOTER_LENGTH;

    assert(dstCapacity > BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH);
    assert(sizeof(TDestValue) == 1u);
    assert(sizeof(unsigned) == 4u);

    // 1. COPY HEADER
    std::ranges::copy(detail::bgzf_compression::magic_header, dstBegin);

    // 2. COMPRESS
    compressInit(ctx);
    ctx.strm.next_in = (Bytef *)(srcBegin);
    ctx.strm.next_out = (Bytef *)(dstBegin + BLOCK_HEADER_LENGTH);
    ctx.strm.avail_in = srcLength * sizeof(TSourceValue);
    ctx.strm.avail_out = dstCapacity - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH;

    int status = deflate(&ctx.strm, Z_FINISH);
    if (status != Z_STREAM_END)
    {
        deflateEnd(&ctx.strm);
        throw io_error("Deflation failed. Compressed BGZF data is too big.");
    }

    status = deflateEnd(&ctx.strm);
    if (status != Z_OK)
        throw io_error("BGZF deflateEnd() failed.");


    // 3. APPEND FOOTER

    // Set compressed length into buffer, compute CRC and write CRC into buffer.

    size_t len = dstCapacity - ctx.strm.avail_out;
    _bgzfPack16(dstBegin + 16, len - 1);

    dstBegin += len - BLOCK_FOOTER_LENGTH;
    _bgzfPack32(dstBegin, crc32(crc32(0u, NULL, 0u), (Bytef *)(srcBegin), srcLength * sizeof(TSourceValue)));
    _bgzfPack32(dstBegin + 4, srcLength * sizeof(TSourceValue));

    return dstCapacity - ctx.strm.avail_out;
}

// ----------------------------------------------------------------------------
// Function decompressInit() - GZIP
// ----------------------------------------------------------------------------

inline void
decompressInit(CompressionContext<detail::gz_compression> & ctx)
{
    const int GZIP_WINDOW_BITS = -15;   // no zlib header

    ctx.strm.zalloc = NULL;
    ctx.strm.zfree = NULL;
    int status = inflateInit2(&ctx.strm, GZIP_WINDOW_BITS);
    if (status != Z_OK)
        throw io_error("GZip inflateInit2() failed.");
}

// ----------------------------------------------------------------------------
// Function decompressInit() - BGZF
// ----------------------------------------------------------------------------

inline void
decompressInit(CompressionContext<detail::bgzf_compression> & ctx)
{
    decompressInit(static_cast<CompressionContext<detail::gz_compression> &>(ctx));
    ctx.headerPos = 0;
}

// ----------------------------------------------------------------------------
// Function _decompressBlock()
// ----------------------------------------------------------------------------

template <typename TDestValue, typename TDestCapacity, typename TSourceValue, typename TSourceLength>
inline TDestCapacity
_decompressBlock(TDestValue *dstBegin,   TDestCapacity dstCapacity,
                 TSourceValue *srcBegin, TSourceLength srcLength, CompressionContext<detail::bgzf_compression> & ctx)
{
    const size_t BLOCK_HEADER_LENGTH = DefaultPageSize<detail::bgzf_compression>::BLOCK_HEADER_LENGTH;
    const size_t BLOCK_FOOTER_LENGTH = DefaultPageSize<detail::bgzf_compression>::BLOCK_FOOTER_LENGTH;

    assert(sizeof(TSourceValue) == 1u);
    assert(sizeof(unsigned) == 4u);

    // 1. CHECK HEADER

    if (srcLength <= BLOCK_HEADER_LENGTH + BLOCK_FOOTER_LENGTH)
        throw io_error("BGZF block too short.");

    if (!detail::bgzf_compression::validate_header(std::span{srcBegin, srcLength}))
        throw io_error("Invalid BGZF block header.");

    size_t compressedLen = _bgzfUnpack16(srcBegin + 16) + 1u;
    if (compressedLen != srcLength)
        throw io_error("BGZF compressed size mismatch.");


    // 2. DECOMPRESS

    decompressInit(ctx);
    ctx.strm.next_in = (Bytef *)(srcBegin + BLOCK_HEADER_LENGTH);
    ctx.strm.next_out = (Bytef *)(dstBegin);
    ctx.strm.avail_in = srcLength - BLOCK_HEADER_LENGTH - BLOCK_FOOTER_LENGTH;
    ctx.strm.avail_out = dstCapacity * sizeof(TDestValue);

    int status = inflate(&ctx.strm, Z_FINISH);
    if (status != Z_STREAM_END)
    {
        inflateEnd(&ctx.strm);
        throw io_error("Inflation failed. Decompressed BGZF data is too big.");
    }

    status = inflateEnd(&ctx.strm);
    if (status != Z_OK)
        throw io_error("BGZF inflateEnd() failed.");


    // 3. CHECK FOOTER

    // Check compressed length in buffer, compute CRC and compare with CRC in buffer.

    unsigned crc = crc32(crc32(0u, NULL, 0u), (Bytef *)(dstBegin), dstCapacity - ctx.strm.avail_out);

    srcBegin += compressedLen - BLOCK_FOOTER_LENGTH;
    if (_bgzfUnpack32(srcBegin) != crc)
        throw io_error("BGZF wrong checksum.");

    if (_bgzfUnpack32(srcBegin + 4) != dstCapacity - ctx.strm.avail_out)
        throw io_error("BGZF size mismatch.");

    return (dstCapacity - ctx.strm.avail_out) / sizeof(TDestValue);
}

}  // namespace seqan3::contrib
