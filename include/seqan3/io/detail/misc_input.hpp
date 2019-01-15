// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
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
    #include <seqan3/contrib/stream/gz_istream.hpp>
#endif
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

/*!\brief Check whether the query range is a prefix of the reference range.
 * \param[in] reference The range that is expected to be the longer one.
 * \param[in] query     The range that is expected to be the shorter one.
 */
template <std::ranges::ForwardRange ref_t, std::ranges::ForwardRange query_t>
inline bool starts_with(ref_t && reference, query_t && query)
//!\cond
    requires std::EqualityComparableWith<reference_t<ref_t>, reference_t<query_t>>
//!\endcond
{
    auto rit  = begin(reference);
    auto rend = end(reference);

    auto qit  = begin(query);
    auto qend = end(query);

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
template <char_concept char_t>
inline auto make_secondary_istream(std::basic_istream<char_t> & primary_stream, filesystem::path & filename)
    -> std::unique_ptr<std::basic_istream<char_t>, std::function<void(std::basic_istream<char_t>*)>>
{
    // don't assume ownership
    constexpr auto stream_deleter_noop     = [] (std::basic_istream<char_t> *) {};
    // assume ownership
    [[maybe_unused]] constexpr auto stream_deleter_default  = [] (std::basic_istream<char_t> * ptr) { delete ptr; };

    // extract "magic header"
    std::istreambuf_iterator<char_t> it{primary_stream};
    constexpr size_t magic_header_size = 4;
    std::array<char, magic_header_size> magic_number;
    for (size_t i = 0; i < magic_header_size; ++i)
    {
        if (it == std::istreambuf_iterator<char_t>{}) // file is smaller than 4 bytes -> not compressed, could be legal
            return {&primary_stream, stream_deleter_noop};

        magic_number[i] = *it;
        ++it;
    }

    // restore original stream content
    for (size_t i = 0; i < magic_header_size; ++i)
        primary_stream.putback(magic_number[magic_header_size - 1 - i]);

    std::string const extension = filename.extension().string();

    // set return value appropriately
    if (starts_with(magic_number, std::array{'\x1f', '\x8b', '\x08'})) // GZIP
    {
    #ifdef SEQAN3_HAS_ZLIB
        if ((extension == ".gz") || (extension == ".bgzf"))
            filename.replace_extension("");

        return {new contrib::basic_gz_istream<char_t>{primary_stream}, stream_deleter_default};
    #else
        throw file_open_error{"Trying to read from a gzipped file, but no ZLIB available."};
    #endif
    }
    else if (starts_with(magic_number, std::array{'\x42', '\x5a', '\x68'})) // BZip2
    {
    #ifdef SEQAN3_HAS_BZIP2
        if (extension == ".bz2")
            filename.replace_extension("");

        return {new contrib::basic_bz2_istream<char_t>{primary_stream}, stream_deleter_default};
    #else
        throw file_open_error{"Trying to read from a bzipped file, but no libbz2 available."};
    #endif
    }
    else if (starts_with(magic_number, std::array{'\x28', '\xb5', '\x2f', '\xfd'})) // ZStd
    {
        throw file_open_error{"Trying to read from a zst'ed file, but SeqAn does not yet support this."};
    }

    return {&primary_stream, stream_deleter_noop};
}

//!\overload
template <char_concept char_t>
inline auto make_secondary_istream(std::basic_istream<char_t> & primary_stream)
{
    filesystem::path p;
    return make_secondary_istream(primary_stream, p);
}

} // namespace seqan3::detail
