// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

#pragma once

#include <string>

/*!\file version.hpp
 * \brief Contains SeqAn version macros and global variables.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

//!\brief The major version as MACRO.
#define SEQAN_VERSION_MAJOR 3
//!\brief The minor version as MACRO.
#define SEQAN_VERSION_MINOR 0
//!\brief The patch version as MACRO.
#define SEQAN_VERSION_PATCH 0

//!\brief The full version as MACRO (number).
#define SEQAN_VERSION (SEQAN_VERSION_MAJOR * 10000 \
                     + SEQAN_VERSION_MINOR * 100 \
                     + SEQAN_VERSION_PATCH)

namespace seqan3
{

//!\brief The major version.
constexpr uint8_t seqan_version_major = SEQAN_VERSION_MAJOR;
//!\brief The minor version.
constexpr uint8_t seqan_version_minor = SEQAN_VERSION_MINOR;
//!\brief The patch version.
constexpr uint8_t seqan_version_patch = SEQAN_VERSION_PATCH;

//!\brief The full version as `std::string`.
std::string const seqan_version = std::to_string(seqan_version_major) + "." +
                                  std::to_string(seqan_version_minor) + "." +
                                  std::to_string(seqan_version_patch);

} // namespace seqan3
