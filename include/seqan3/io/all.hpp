// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Meta-header for all IO related functionality.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \author René Rahn <rene.rahn AT fu-berlin.de>
 */

/*!\defgroup io IO
 * \brief The IO module provides stream handling formatted I/O.
 *
 * # stream_REMOVEMEs and (de-)compression {#io_compression}
 *
 * SeqAn works with regular iostreams as provided by the standard library, but it also handles compressed streams:
 *
 * | **Format** | **Extension**   | **Dependency**                   | **Description**                                   |
 * |:-----------|:----------------|:---------------------------------|:--------------------------------------------------|
 * | GZip       | `.gz`¹          | [zlib](https://zlib.net/)        | GNU-Zip, most common format on UNIX               |
 * | BGZF       | `.gz`, `.bgzf`² | [zlib](https://zlib.net/)        | [Blocked GZip](https://samtools.github.io/hts-specs/SAMv1.pdf), compatible extension to GZip, features parallelisation|
 * | BZip2      | `.bz2`          | [libbz2](https://www.bzip.org)   | Stronger compression than GZip, slower to compress |
 *
 * <small>¹ SeqAn always assumes GZip and does not handle pure `.Z`.<br>
 * ² Some file formats like `.bam` or `.bcf` are implicitly BGZF-compressed without showing this in the
 * extension.</small>
 *
 * Support for these compression formats is **optional** and depends on whether the respective dependency is available
 * when you build your program (if you use CMake, this should happen automatically).
 *
 * SeqAn file types apply compression/decompression streams transparently, i.e. if the given file-extension or
 * "magic-header" of a file suggest this, the respective stream is automatically (de-)compressed.
 *
 * The (de)compression stream wrappers are currently only used internally and not part of the API.
 *
 * # Formatted I/O
 *
 * ## Files and formats {#io_files}
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
 * Some formats are available in multiple files, e.g. seqan3::format_sam can be read by seqan3::sequence_file_input
 * and by seqan3::alignment_file_input. This represents different use-cases of the same file format.
 *
 * Typically formats are supported for reading and writing, but this does not always have to be the case. See the above
 * links for more information.
 *
 * ## Records and fields
 *
 * The main file interface that SeqAn offers is *record-based*, i.e. every file conceptionally is a range of
 * records. And each record in turn behaves as a tuple of fields.
 *
 * The record type of all files is based on seqan3::record, but the composition of fields is different **per file**.
 *
 * In particular this means:
 *
 *   * You can iterate over a seqan3::sequence_file_input just like you iterate over a std::vector.
 *   * The element type is seqan3::record and for seqan3::sequence_file_input the records typically consist of three
 *     fields: ID, sequence and qualities.
 *   * Not every *format* may provide every field (e.g. all three fields can be extracted from a FastQ file, but only
 *     ID and sequence can be extracted from a FastA file).
 *   * This is not a problem when reading files; the format is detected at run-time and either it provides all desired
 *     fields, or some of them will be empty.
 *   * When writing files it depends on the format whether it allows certain fields to be unset; e.g. you can convert
 *     a FastQ file to a FastA file, but not the other way around (at least not without providing placeholders for the
 *     qualitites).
 *
 * Please have a look the tutorial for \ref tutorial_sequence_file and the API docs for seqan3::sequence_file_input
 * to learn about this design in practice.
 *
 * # Serialisation {#serialisation}
 *
 * Besides formatted I/O which is realised via files and formats, SeqAn also supports object-level serialisation.
 * This enables you to store data structures like indexes or sequences directly to disk.
 *
 * We use the [cereal library](https://github.com/USCiLab/cereal) to accomplish this.
 * For more information see [cereal's documentation](https://uscilab.github.io/cereal/) or our tutorial on
 * \ref tutorial_index_search which contains an example.
 */

#pragma once

#include <seqan3/io/alignment_file/all.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/io/record.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/io/stream/all.hpp>
#include <seqan3/io/structure_file/all.hpp>
