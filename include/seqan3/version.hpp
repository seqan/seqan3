// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

#include <string>

#include <seqan3/core/platform.hpp>

/*!\file
 * \brief Provides SeqAn version macros and global variables.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

//!\brief The major version as MACRO.
#define SEQAN3_VERSION_MAJOR 3
//!\brief The minor version as MACRO.
#define SEQAN3_VERSION_MINOR 0
//!\brief The patch version as MACRO.
#define SEQAN3_VERSION_PATCH 2

//!\brief The full version as MACRO (number).
#define SEQAN3_VERSION (SEQAN3_VERSION_MAJOR * 10000 \
                     + SEQAN3_VERSION_MINOR * 100 \
                     + SEQAN3_VERSION_PATCH)

namespace seqan3
{

//!\brief The major version.
constexpr uint8_t seqan3_version_major = SEQAN3_VERSION_MAJOR;
//!\brief The minor version.
constexpr uint8_t seqan3_version_minor = SEQAN3_VERSION_MINOR;
//!\brief The patch version.
constexpr uint8_t seqan3_version_patch = SEQAN3_VERSION_PATCH;

//!\brief The full version as `std::string`.
std::string const seqan3_version = std::to_string(seqan3_version_major) + "." +
                                   std::to_string(seqan3_version_minor) + "." +
                                   std::to_string(seqan3_version_patch);

} // namespace seqan3
