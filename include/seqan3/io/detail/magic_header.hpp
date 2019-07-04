// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::magic_header.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <array>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

//!\brief Defines a magic byte sequence to disambiguate different compression formats. Default is empty.
//!\ingroup io
template <typename header_tag_t>
inline constexpr std::array<char, 0> magic_header{};

//!\brief A tag signifying a gz compressed file.
//!\ingroup io
struct gz_compression
{};

/*!\brief The magic byte sequence to disambiguate gz compressed files.
 * \ingroup io
 *
 * \details
 *
 * Specialises seqan3::detail::magic_header for seqan3::detail::gz_compression.
 */
template <>
inline constexpr std::array<char, 3> magic_header<gz_compression>{'\x1f', '\x8b', '\x08'};

//!\brief A tag signifying a bz2 compressed file.
//!\ingroup io
struct bz2_compression
{};

/*!\brief The magic byte sequence to disambiguate bz2 compressed files.
 * \ingroup io
 *
 * \details
 *
 * Specialises seqan3::detail::magic_header for seqan3::detail::bz2_compression.
 */
template <>
inline constexpr std::array<char, 3> magic_header<bz2_compression>{'\x42', '\x5a', '\x68'};

//!\brief A tag signifying a zstd compressed file.
//!\ingroup io
struct zstd_compression
{};

/*!\brief The magic byte sequence to disambiguate zstd compressed files.
 * \ingroup io
 *
 * \details
 *
 * Specialises seqan3::detail::magic_header for seqan3::detail::zstd_compression.
 */
template <>
inline constexpr std::array<char, 4> magic_header<zstd_compression>{'\x28', '\xb5', '\x2f', '\xfd'};

//!\brief A tag signifying a bgzf compressed file.
//!\ingroup io
struct bgzf_compression
{};

/*!\brief The magic byte sequence to disambiguate bgzf compressed files.
 * \ingroup io
 *
 * \details
 *
 * Specialises seqan3::detail::magic_header for seqan3::detail::bgzf_compression.
 */
template <>
inline constexpr std::array<char, 18> magic_header<bgzf_compression>
{
//  ID1                              ID2                              CM
    magic_header<gz_compression>[0], magic_header<gz_compression>[1], magic_header<gz_compression>[2],
//  FLG     [MTIME                       ]  XFL     OS      [XLEN        ]
    '\x04', '\x00', '\x00', '\x00', '\x00', '\x00', '\xff', '\x06', '\x00',
//  B       C       [SLEN        ]  [BSIZE       ]
    '\x42', '\x43', '\x02', '\x00', '\x00', '\x00'
};

} // namespace seqan3::detail
