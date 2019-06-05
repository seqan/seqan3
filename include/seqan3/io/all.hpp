// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
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
 * # Streams and (de-)compression {#io_compression}
 *
 * SeqAn works with regular iostreams as provided by the standard library, but it also handles compressed streams:
 *
 * | **Format** | **Extension**   | **Dependency**                   | **Description**                                   |
 * |:-----------|:----------------|:---------------------------------|:--------------------------------------------------|
 * | GZip       | `.gz`¹          | [zlib](https://zlib.net/)        | GNU-Zip, most common format on UNIX               |
 * | BGZF       | `.gz`, `.bgzf`² | [zlib](https://zlib.net/)        | [Blocked GZip](http://samtools.github.io/hts-specs/SAMv1.pdf), compatible extension to GZip, features parallelisation|
 * | BZip2      | `.bz2`          | [libbz2](https://www.bzip.org)   | Stronger compression than GZip, slower to compress |
 *
 * <small>¹ SeqAn always assumes GZip and does not handle pure `.Z`.<br>
 * ² Some file formats like `.bam` or `.bcf` are implicitly BGZF-compressed without showing this in the
 * extension.</small>
 *
 * Support for these compression format is **optional** and depends on whether the respective dependency is available
 * when you build your program (if you use CMake, this should happen automatically).
 *
 * SeqAn file types apply compression/decompression streams transparently, i.e. if the given file-extension or
 * "magic-header" of a file suggest this, the respective stream is automatically (de-)compressed.
 *
 * The (de)compression stream wrappers are currently only used internally and not part of the API.
 *
 * # Files and formats {#io_files}
 *
 * SeqAn has the notion of *files* and *formats*. *File* is an abstraction level higher than *format*.
 * A file describes a common use-case and it typically supports multiple *formats*.
 * The developer needs to know which kind of file they want to read/write, this choice is made at compile-time.
 * The format, on the other hand, is automatically detected based on the file provided by the user to the program.
 *
 * For example, seqan3::sequence_file_input handles reading sequence files.
 * It can be created directly from an input stream, or from a file name. After opening the file it will detect whether
 * the format is seqan3::format_fasta or seqan3::format_fastq (or another supported format) automatically – normally
 * by comparing the extension.
 *
 * | **File**                      | **Formats**                                                                                                  |
 * |:------------------------------|:-------------------------------------------------------------------------------------------------------------|
 * | seqan3::alignment_file_input  | seqan3::format_sam, seqan3::format_bam                                                                       |
 * | seqan3::alignment_file_output | seqan3::format_sam, seqan3::format_bam                                                                       |
 * | seqan3::sequence_file_input   | seqan3::format_embl, seqan3::format_fasta, seqan3::format_fastq, seqan3::format_genbank, seqan3::format_sam  |
 * | seqan3::sequence_file_output  | seqan3::format_embl, seqan3::format_fasta, seqan3::format_fastq, seqan3::format_genbank, seqan3::format_sam  |
 * | seqan3::structure_file_input  | seqan3::format_vienna                                                                                        |
 * | seqan3::structure_file_output | seqan3::format_vienna                                                                                        |
 *
 *
 * Some formats are available in multiple files, e.g. seqan3::format_sam can be read by seqan3::sequence_file_input
 * and by seqan3::alignment_file_input. This represents different use-cases of the same file.
 *
 * Typically formats are supported for reading and writing, but this does not always have to be the case. See the above
 * links for more information.
 *
 * # Records and fields
 *
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
