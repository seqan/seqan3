// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <cstddef>
#include <cstdint>

/*!\file
 * \brief Provides SeqAn version macros and global variables.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

//!\brief The major version as MACRO.
#define SEQAN3_VERSION_MAJOR 3
//!\brief The minor version as MACRO.
#define SEQAN3_VERSION_MINOR 0
//!\brief The patch version as MACRO.
#define SEQAN3_VERSION_PATCH 3

//!\brief The full version as MACRO (number).
#define SEQAN3_VERSION (SEQAN3_VERSION_MAJOR * 10000 \
                     + SEQAN3_VERSION_MINOR * 100 \
                     + SEQAN3_VERSION_PATCH)

/*!\brief Converts a number to a string. Preprocessor needs this indirection to
 * properly expand the values to strings.
 */
#define SEQAN3_VERSION_CSTRING_HELPER_STR(str) #str

//!\brief Converts version numbers to string.
#define SEQAN3_VERSION_CSTRING_HELPER_FUNC(MAJOR, MINOR, PATCH) \
    SEQAN3_VERSION_CSTRING_HELPER_STR(MAJOR) "."\
    SEQAN3_VERSION_CSTRING_HELPER_STR(MINOR) "."\
    SEQAN3_VERSION_CSTRING_HELPER_STR(PATCH)

//!\brief The full version as null terminated string.
#define SEQAN3_VERSION_CSTRING \
    SEQAN3_VERSION_CSTRING_HELPER_FUNC(SEQAN3_VERSION_MAJOR, \
                                       SEQAN3_VERSION_MINOR, \
                                       SEQAN3_VERSION_PATCH)

namespace seqan3
{

//!\brief The major version.
constexpr uint8_t seqan3_version_major = SEQAN3_VERSION_MAJOR;
//!\brief The minor version.
constexpr uint8_t seqan3_version_minor = SEQAN3_VERSION_MINOR;
//!\brief The patch version.
constexpr uint8_t seqan3_version_patch = SEQAN3_VERSION_PATCH;

//!\brief The full version as `std::size_t`.
constexpr std::size_t seqan3_version = SEQAN3_VERSION;

//!\brief The full version as null terminated string.
constexpr char const* seqan3_version_cstring = SEQAN3_VERSION_CSTRING;

} // namespace seqan3

#undef SEQAN3_VERSION_CSTRING_HELPER_STR
#undef SEQAN3_VERSION_CSTRING_HELPER_FUNC
