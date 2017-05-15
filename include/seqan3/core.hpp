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

/*!\file core.hpp
 * \ingroup core
 * \brief Meta-header for the \link core core module \endlink.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

// ============================================================================
// External concept implementations
// ============================================================================

#include <seqan3/core/concept/core.hpp>
#include <seqan3/core/concept/iterator.hpp>
#include <seqan3/core/detail/int_types.hpp>
#include <seqan3/core/pod_tuple.hpp>
#include <seqan3/core/platform.hpp>

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

/*!\namespace seqan3::literal
 * \brief The SeqAn3 namespace for literals.
 *
 * SeqAn implements "user defined" literals in multiple places, e.g. `auto foo = "ACGTG"_dna4`. These
 * make working with small examples and tests a lot easier, but the risk of having a name collision with
 * another library is higher so follow the example of the standard library and define all our literals
 * in the namespace `seqan3::literal`.
 *
 * \attention
 * This means you cannot use them, unless you explicitly add `using namespace seqan3::literal;` (in addition
 * to `using namespace seqan3;`).
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
