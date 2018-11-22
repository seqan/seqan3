// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
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

/*!\file
 * \brief Meta-header for the \link core core module \endlink.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

// ============================================================================
// External concept implementations
// ============================================================================

#include <seqan3/core/algorithm/all.hpp>
#include <seqan3/core/add_enum_bitwise_operators.hpp>
#include <seqan3/core/concept/all.hpp>
#include <seqan3/core/detail/all.hpp>
#include <seqan3/core/metafunction/all.hpp>
#include <seqan3/core/platform.hpp>
#include <seqan3/core/pod_tuple.hpp>
#include <seqan3/core/tuple_utility.hpp>
#include <seqan3/core/type_list.hpp>

/*!\defgroup core Core
 * \brief Contains core functionality used by multiple modules.
 *
 * The core module contains concepts, functions and some classes that
 * are used by multiple other modules, but that usually are not relevant
 * to most users of the library.
 */

/*!\namespace seqan3
 * \brief The main SeqAn3 namespace.
 */

/*!\cond DEV
 * \namespace seqan3::detail
 * \brief The internal SeqAn3 namespace.
 * \details
 * The contents of this namespace are not visible to consumers of the library and the documentation is
 * only generated for developers.
 * \sa https://github.com/seqan/seqan3/wiki/Documentation
 * \endcond
 */

/*!\namespace std
 * \brief SeqAn specific customisations in the standard namespace.
 */

/*!\mainpage SeqAn3 API documentation
 *
 * This is the API documentation for the SeqAn3 library.
 *
 * We recommend browsing the documentation via the "Modules" menu on the left.
 *
 * Installation instructions are available on [GitHub](https://github.com/seqan/seqan3).
 *
 * A manual will be added in the future.
 */
