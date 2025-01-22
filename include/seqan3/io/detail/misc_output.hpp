// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides various utility functions required only for output.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <filesystem>
#include <functional>
#include <iostream>
#include <string>
#include <tuple>

#if defined(SEQAN3_HAS_BZIP2)
#    include <seqan3/contrib/stream/bz2_ostream.hpp>
#endif
#if defined(SEQAN3_HAS_ZLIB)
#    include <seqan3/contrib/stream/bgzf_ostream.hpp>
#    include <seqan3/contrib/stream/gz_ostream.hpp>
#endif
#include <seqan3/contrib/stream/bgzf.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/utility/concept.hpp>

namespace seqan3::detail
{

/*!\brief Depending on the given filename/extension, create a compression stream or just forward the primary stream.
 * \ingroup io
 * \param[in] primary_stream The primary (uncompressed) stream for writing.
 * \param[in,out] filename  The associated filename; compression extensions will be stripped.
 * \returns A pointer to the secondary stream with a default deleter or a nop-deleter.
 * \throws seqan3::file_open_error If a compression-extension is used, but is not supported/available.
 */
template <builtin_character char_t>
inline auto make_secondary_ostream(std::basic_ostream<char_t> & primary_stream, std::filesystem::path & filename)
    -> std::unique_ptr<std::basic_ostream<char_t>, std::function<void(std::basic_ostream<char_t> *)>>
{
    // don't assume ownership
    constexpr auto stream_deleter_noop = [](std::basic_ostream<char_t> *) {};
    // assume ownership
    [[maybe_unused]] constexpr auto stream_deleter_default = [](std::basic_ostream<char_t> * ptr)
    {
        delete ptr;
    };

    std::string extension = filename.extension().string();

    if (extension == ".gz")
    {
#if defined(SEQAN3_HAS_ZLIB)
        filename.replace_extension("");
        return {new contrib::basic_gz_ostream<char_t>{primary_stream}, stream_deleter_default};
#else
        throw file_open_error{"Trying to write a gzipped file, but no ZLIB available."};
#endif
    }
    else if ((extension == ".bgzf") || (extension == ".bam"))
    {
#if defined(SEQAN3_HAS_ZLIB)
        if (extension != ".bam") // remove extension except for bam
            filename.replace_extension("");

        return {new contrib::basic_bgzf_ostream<char_t>{primary_stream}, stream_deleter_default};
#else
        throw file_open_error{"Trying to write a bgzf'ed file, but no ZLIB available."};
#endif
    }
    else if (extension == ".bz2")
    {
#if defined(SEQAN3_HAS_BZIP2)
        filename.replace_extension("");
        return {new contrib::basic_bz2_ostream<char_t>{primary_stream}, stream_deleter_default};
#else
        throw file_open_error{"Trying to write a bzipped file, but no libbz2 available."};
#endif
    }
    else if (extension == ".zst")
    {
        throw file_open_error{"Trying to write a zst'ed file, but SeqAn does not yet support this."};
    }

    return {&primary_stream, stream_deleter_noop};
}

} // namespace seqan3::detail
