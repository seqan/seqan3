// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

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
 * \ingroup utility_type_list
 */
template <typename... types>
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
