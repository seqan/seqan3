// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::type_list.
 */

#pragma once

#include <meta/meta.hpp>

#include <seqan3/core/type_traits/template_inspection.hpp>

namespace seqan3
{

// ----------------------------------------------------------------------------
// type_list class
// ----------------------------------------------------------------------------

/*!\brief Type that contains multiple types, an alias for
 * [meta::list](https://ericniebler.github.io/range-v3/structmeta_1_1list.html).
 * \ingroup type_list
 */
template <typename ...types>
using type_list = meta::list<types...>;

} // namespace seqan3

namespace seqan3::detail
{

// ----------------------------------------------------------------------------
// type_list_specialisation concept
// ----------------------------------------------------------------------------

/*!\brief Auxiliary concept that checks whether a type is a specialisation of seqan3::type_list.
 * \ingroup type_list
 */
template <typename t>
SEQAN3_CONCEPT type_list_specialisation = is_type_specialisation_of_v<t, type_list>;

} // namespace seqan3::detail
