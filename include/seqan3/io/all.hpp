// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

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
#include <seqan3/io/record.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/io/stream/all.hpp>
#include <seqan3/io/structure_file/all.hpp>
