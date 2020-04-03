// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::magic_header.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <array>
#include <string>
#include <vector>

#include <seqan3/core/bit_manipulation.hpp>
#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/core/type_list/type_list.hpp>
#include <seqan3/core/type_list/traits.hpp>
#include <seqan3/std/type_traits>
#include <seqan3/std/span>

namespace seqan3::detail
{

//!\brief A tag signifying a gz compressed file.
//!\ingroup io
struct gz_compression
{
    //!\brief The valid file extension for gz compression.
    static inline std::vector<std::string> file_extensions
    {
        {"gz"}
    };

    //!\brief The magic byte sequence to disambiguate gz compressed files.
    static constexpr std::array<char, 3> magic_header{'\x1f', '\x8b', '\x08'};
};

//!\brief A tag signifying a bz2 compressed file.
//!\ingroup io
struct bz2_compression
{
    //!\brief The valid file extension for bz2 compression.
    static inline std::vector<std::string> file_extensions
    {
        {"bz2"}
    };

    //!\brief The magic byte sequence to disambiguate bz2 compressed files.
    static constexpr std::array<char, 3> magic_header{'\x42', '\x5a', '\x68'};
};

//!\brief A tag signifying a zstd compressed file.
//!\ingroup io
struct zstd_compression
{
    //!\brief The valid file extension for zstd compression.
    static inline std::vector<std::string> file_extensions
    {
        {"zst"}
    };

    //!\brief The magic byte sequence to disambiguate zstd compressed files.
    static constexpr std::array<char, 4> magic_header{'\x28', '\xb5', '\x2f', '\xfd'};
};

//!\brief A tag signifying a bgzf compressed file.
//!\ingroup io
struct bgzf_compression
{
    //!\brief The valid file extension for bgzf compression.
    static inline std::vector<std::string> file_extensions
    {
        {"bgzf"}
    };

    //!\brief The magic byte sequence to disambiguate bgzf compressed files.
    static constexpr std::array<char, 18> magic_header
    {
    //  ID1                              ID2                              CM
        gz_compression::magic_header[0], gz_compression::magic_header[1], gz_compression::magic_header[2],
    //  FLG     [MTIME                       ]  XFL     OS      [XLEN        ]
        '\x04', '\x00', '\x00', '\x00', '\x00', '\x00', '\xff', '\x06', '\x00',
    //  B       C       [SLEN        ]  [BSIZE       ]
        '\x42', '\x43', '\x02', '\x00', '\x00', '\x00'
    };

    /*!\brief Checks if the given header is a bgzf header.
     * \param[in] header The header to validate.
     * \returns `true` if it is a bgzf header, otherwise `false`.
     */
    template <typename char_t, size_t extend>
    static bool validate_header(std::span<char_t, extend> header)
    {
        static_assert(seqan3::detail::weakly_equality_comparable_with<char_t, char>,
                      "The given char type of the span must be comparable with char.");

        return (header[0] == magic_header[0] &&                                                     // GZ_ID1
                header[1] == magic_header[1] &&                                                     // GZ_ID2
                header[2] == magic_header[2] &&                                                     // GZ_CM
                (header[3] & magic_header[3]) != 0 &&                                                // FLG_FEXTRA
                to_little_endian(*reinterpret_cast<uint16_t *>(&header[10])) == magic_header[10] &&  // BGZF_ID1
                header[12] == magic_header[12] &&                                                    // BGZF_ID2
                header[13] == magic_header[13] &&                                                    // BGZF_SLEN
                to_little_endian(*reinterpret_cast<uint16_t *>(&header[14])) == magic_header[14]);   // BGZF_XLEN
    }
};

/*!\brief A seqan3::type_list containing the available compression formats.
 * \ingroup io
 */
using compression_formats = pack_traits::drop_front<void
                                                    #if SEQAN3_HAS_ZLIB
                                                    , gz_compression
                                                    , bgzf_compression
                                                    #endif // SEQAN3_HAS_ZLIB
                                                    #if SEQAN3_HAS_BZIP2
                                                    , bz2_compression
                                                    #endif // SEQAN3_HAS_BZIP2
                                                    #if SEQAN3_HAS_ZSTD
                                                    , zstd_compression
                                                    #endif // SEQAN3_HAS_ZSTD
                                                    >;

} // namespace seqan3::detail
