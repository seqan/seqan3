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

#include <seqan3/core/detail/template_inspection.hpp>

namespace seqan3
{

// ----------------------------------------------------------------------------
// type_list class
// ----------------------------------------------------------------------------

/*!\brief Type that contains multiple types.
 * \ingroup type_list
 */
template <typename ...types>
struct type_list
{
    //!\brief The type list itself
    using type = type_list;

    //!\brief The number of types contained in the type list
    static constexpr size_t size() noexcept
    {
        return sizeof...(types);
    }
};

} // namespace seqan3
