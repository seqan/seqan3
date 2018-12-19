// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

/*!\file
 * \brief Meta-header for all IO related functionality.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \author René Rahn <rene.rahn AT fu-berlin.de>
 */

/*!\defgroup io IO
 * \brief The IO module contains concepts, data structures and functions related to reading and writing formatted files,
 * streams, and serialisation.
 *
 *
 * # Operations on characters and strings {#char_ops}
 *
 * \todo write me!
 *
 *
 * # Formatted files {#io_files}
 *
 * \todo write me!
 *
 *
 * # Compression and decompression {#io_compression}
 *
 * SeqAn supports wrapping iostreams in streams that transparently compress /
 * decompress the data sent to / read from the stream.
 * The following compression formats are supported:
 *
 * | **Format** | **Extension**   | **Dependency** | **Description**                                   |
 * |:-----------|:----------------|:----------------|:--------------------------------------------------|
 * | GZip       | `.gz`¹          | [zlib](https://zlib.net/)        | GNU-Zip, most common format on UNIX               |
 * | BGZF       | `.gz`, `.bgzf`² | [zlib](https://zlib.net/)        | [Blocked GZip](http://samtools.github.io/hts-specs/SAMv1.pdf), compatible extension to GZip, features parallelisation|
 * | BZip2      | `.bz2`          | [libbz2](https://www.bzip.org)   | Stronger compression than GZip, slower to compress |
 *
 * <small>¹ SeqAn always assumes GZip and does not handle pure `.Z`.<br>
 * ²Some file formats like `.bam` or `.bcf` are implicitly BGZF-compressed without showing this in the
 * extension.</small>
 *
 * Whether or not a format is supported depends on whether the respective dependency was available and
 * linked against your program (if you use CMake, this should happen automatically).
 *
 * Formatted files employ compression/decompression streams transparently, i.e. if the given file-extension or
 * "magic-header" of a file suggest this, the respective stream is automatically (de-)compressed.
 *
 * # Serialisation {#serialisation}
 *
 * \todo write me!
 */

#pragma once

#include <seqan3/io/alignment_file/all.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/io/filesystem.hpp>
#include <seqan3/io/record.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/io/stream/all.hpp>
#include <seqan3/io/structure_file/all.hpp>
