// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides various utility functions required only for input.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <iostream>
#include <string>
#include <tuple>

#include <seqan3/core/concept/core_language.hpp>
#ifdef SEQAN3_HAS_BZIP2
    #include <seqan3/contrib/stream/bz2_istream.hpp>
#endif
#ifdef SEQAN3_HAS_ZLIB
    #include <seqan3/contrib/stream/bgzf_istream.hpp>
    #include <seqan3/contrib/stream/bgzf_stream_util.hpp>
    #include <seqan3/contrib/stream/gz_istream.hpp>
#endif
#include <seqan3/io/detail/magic_header.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/concepts>
#include <seqan3/std/filesystem>
#include <seqan3/std/ranges>
#include <seqan3/std/span>

namespace seqan3::detail
{

/*!\brief Check whether the query range is a prefix of the reference range.
 * \param[in] reference The range that is expected to be the longer one.
 * \param[in] query     The range that is expected to be the shorter one.
 */
template <std::ranges::forward_range ref_t, std::ranges::forward_range query_t>
inline bool starts_with(ref_t && reference, query_t && query)
//!\cond
    requires std::equality_comparable_with<std::ranges::range_reference_t<ref_t>,
                                           std::ranges::range_reference_t<query_t>>
//!\endcond
{
    auto rit  = std::ranges::begin(reference);
    auto rend = std::ranges::end(reference);

    auto qit  = std::ranges::begin(query);
    auto qend = std::ranges::end(query);

    while (true)
    {
        if (qit == qend)
            return true;

        if (rit == rend)
            return false;

        if (*qit != *rit)
            return false;

        ++qit;
        ++rit;
    }
}

/*!\brief Depending on the magic bytes of the given stream, return a decompression stream or forward the primary stream.
 * \param[in] primary_stream The primary (device) stream for reading.
 * \param[in,out] filename  The associated filename; compression extensions will be stripped. [optional]
 * \returns A pointer to the secondary stream with defaulted or NOP'ed deleter.
 * \throws seqan3::file_open_error If the magic bytes suggest compression, but is not supported/available.
 */
template <builtin_character char_t>
inline auto make_secondary_istream(std::basic_istream<char_t> & primary_stream, std::filesystem::path & filename)
    -> std::unique_ptr<std::basic_istream<char_t>, std::function<void(std::basic_istream<char_t>*)>>
{
    assert(primary_stream.good());

    // don't assume ownership
    constexpr auto stream_deleter_noop     = [] (std::basic_istream<char_t> *) {};
    // assume ownership
    [[maybe_unused]] constexpr auto stream_deleter_default  = [] (std::basic_istream<char_t> * ptr) { delete ptr; };

    // extract "magic header"
    std::istreambuf_iterator<char_t> it{primary_stream};
    std::array<char, bgzf_compression::magic_header.size()> magic_number{}; // Largest magic header from bgzf
    size_t read_chars = 0;
    for (; read_chars < magic_number.size(); ++read_chars)
    {
        if (it == std::istreambuf_iterator<char_t>{})
            break;

        magic_number[read_chars] = *it;
        ++it;
    }

    // unget all read chars.
    for (size_t i = 0 ; i < read_chars; ++i)
        primary_stream.unget();

    std::string extension{};
    if (filename.has_extension())
        extension = filename.extension().string().substr(1);

    // tests whether the given extension matches with one of the given compression tags.
    [[maybe_unused]] auto contains_extension = [] (auto compression_tag, auto const & extension) constexpr
    {
        return std::ranges::find(decltype(compression_tag)::file_extensions, extension) !=
               std::ranges::end(decltype(compression_tag)::file_extensions);
    };

    // set return value appropriately
    if (read_chars == magic_number.size() && bgzf_compression::validate_header(std::span{magic_number})) // BGZF
    {
    #ifdef SEQAN3_HAS_ZLIB
        if (contains_extension(gz_compression{}, extension) || contains_extension(bgzf_compression{}, extension))
            filename.replace_extension();

        return {new contrib::basic_bgzf_istream<char_t>{primary_stream},
                stream_deleter_default};
    #else
        throw file_open_error{"Trying to read from a bgzf file, but no ZLIB available."};
    #endif
    }
    else if (starts_with(magic_number, gz_compression::magic_header)) // GZIP
    {
    #ifdef SEQAN3_HAS_ZLIB
        if (contains_extension(gz_compression{}, extension) || contains_extension(bgzf_compression{}, extension))
            filename.replace_extension();

        return {new contrib::basic_gz_istream<char_t>{primary_stream}, stream_deleter_default};
    #else
        throw file_open_error{"Trying to read from a gzipped file, but no ZLIB available."};
    #endif
    }
    else if (starts_with(magic_number, bz2_compression::magic_header)) // BZip2
    {
    #ifdef SEQAN3_HAS_BZIP2
        if (contains_extension(bz2_compression{}, extension))
            filename.replace_extension();

        return {new contrib::basic_bz2_istream<char_t>{primary_stream}, stream_deleter_default};
    #else
        throw file_open_error{"Trying to read from a bzipped file, but no libbz2 available."};
    #endif
    }
    else if (starts_with(magic_number, zstd_compression::magic_header)) // ZStd
    {
        throw file_open_error{"Trying to read from a zst'ed file, but SeqAn does not yet support this."};
    }

    return {&primary_stream, stream_deleter_noop};
}

//!\overload
template <builtin_character char_t>
inline auto make_secondary_istream(std::basic_istream<char_t> & primary_stream)
{
    std::filesystem::path p;
    return make_secondary_istream(primary_stream, p);
}

} // namespace seqan3::detail
