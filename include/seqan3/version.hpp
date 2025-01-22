// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

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
#define SEQAN3_VERSION_MINOR 4
//!\brief The patch version as MACRO.
#define SEQAN3_VERSION_PATCH 0
//!\brief The release candidate number. 0 means stable release, >= 1 means release candidate.
#define SEQAN3_RELEASE_CANDIDATE 4

//!\brief The full version as MACRO (number).
#define SEQAN3_VERSION (SEQAN3_VERSION_MAJOR * 10000 + SEQAN3_VERSION_MINOR * 100 + SEQAN3_VERSION_PATCH)

/*!\brief Converts a number to a string. Preprocessor needs this indirection to
 * properly expand the values to strings.
 */
#define SEQAN3_VERSION_CSTRING_HELPER_STR(str) #str

//!\brief Converts version numbers to string.
#define SEQAN3_VERSION_CSTRING_HELPER_FUNC(MAJOR, MINOR, PATCH)                                                        \
    SEQAN3_VERSION_CSTRING_HELPER_STR(MAJOR)                                                                           \
    "." SEQAN3_VERSION_CSTRING_HELPER_STR(MINOR) "." SEQAN3_VERSION_CSTRING_HELPER_STR(PATCH)

#if (SEQAN3_RELEASE_CANDIDATE > 0)
//!\brief A helper function that expands to a suitable release candidate suffix.
#    define SEQAN3_RELEASE_CANDIDATE_HELPER(RC) "-rc." SEQAN3_VERSION_CSTRING_HELPER_STR(RC)
#else
//!\brief A helper function that expands to a suitable release candidate suffix.
#    define SEQAN3_RELEASE_CANDIDATE_HELPER(RC) ""
#endif

//!\brief The full version as null terminated string.
#define SEQAN3_VERSION_CSTRING                                                                                         \
    SEQAN3_VERSION_CSTRING_HELPER_FUNC(SEQAN3_VERSION_MAJOR, SEQAN3_VERSION_MINOR, SEQAN3_VERSION_PATCH)               \
    SEQAN3_RELEASE_CANDIDATE_HELPER(SEQAN3_RELEASE_CANDIDATE)

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
constexpr char const * seqan3_version_cstring = SEQAN3_VERSION_CSTRING;

} // namespace seqan3

#undef SEQAN3_VERSION_CSTRING_HELPER_STR
#undef SEQAN3_VERSION_CSTRING_HELPER_FUNC
#undef SEQAN3_RELEASE_CANDIDATE_HELPER
