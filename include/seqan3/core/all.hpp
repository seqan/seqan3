// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Meta-header for the \link core core module \endlink.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

// ============================================================================
// External concept implementations
// ============================================================================

#include <seqan3/core/add_enum_bitwise_operators.hpp>
#include <seqan3/core/algorithm/all.hpp>
#include <seqan3/core/concept/all.hpp>
#include <seqan3/core/configuration/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/core/detail/all.hpp>
#include <seqan3/core/platform.hpp>
#include <seqan3/core/range/type_traits.hpp>
#include <seqan3/core/tuple_utility.hpp>
#include <seqan3/utility/tuple/pod_tuple.hpp>

/*!\defgroup core Core
 * \brief Provides core functionality used by multiple modules.
 *
 * The core module contains concepts, functions and some classes that
 * are used by multiple other modules, but that usually are not relevant
 * to most users of the library.
 */

/*!\namespace seqan3
 * \brief The main SeqAn3 namespace.
 */

/*!\namespace seqan3::custom
 * \brief A namespace for third party and standard library specialisations of SeqAn customisation points.
 * \see \ref about_customisation
 */

 /*!\namespace seqan3::exposition_only
  * \brief A namespace for SeqAn entities that are intended for documentation purposes only.
  *
  * \details
  *
  * We try to guarantee a strong \ref api_stability "API-Stability". This has the downside that we sometimes opt to not
  * document entities, because
  *
  * 1. we do not want to guarantee that the entity will always keep the same name or
  * 2. we need a helper entity to define the actual entity or
  * 3. we are not completely sure if the entity is general enough to mark it stable,
  *
  * but we still want to express a general look-and-feel of the API in our documentation.
  *
  * We therefore use the same trick as the C++ standard that defines entities, s.a. concepts, in a partial defined
  * state to express a general intend of the API without being explicit about it.
  *
  * For example, see https://eel.is/c++draft/iterator.concept.readable where `indirectly-readable-impl` describes the
  * general intention of the concept, but does not name it since it is a helper-entity for the std::indirectly_­readable
  * concept.
  */

/*!\if DEV
 * \namespace seqan3::detail
 * \brief The internal SeqAn3 namespace.
 * \details
 * The contents of this namespace are not visible to consumers of the library and the documentation is
 * only generated for developers.
 * \sa https://github.com/seqan/seqan3/wiki/Documentation
 * \endif
 */

/*!\namespace std
 * \brief SeqAn specific customisations in the standard namespace.
 */
