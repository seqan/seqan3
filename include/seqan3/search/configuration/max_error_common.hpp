// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides the error types for maximum number of errors.
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \author Lydia Buntrock <lydia.buntrock AT fu-berlin.de>
 */

#pragma once

#include <concepts>

#include <seqan3/core/detail/strong_type.hpp>
#include <seqan3/utility/concept.hpp>

namespace seqan3::search_cfg
{

/*!\brief A strong type of underlying type `uint8_t` that represents the number of errors.
 * \ingroup search_configuration
 * \see search_configuration
 * \tparam value_t The underlying type.
 */
struct error_count :
    detail::
        strong_type<uint8_t, error_count, detail::strong_type_skill::convert | detail::strong_type_skill::comparable>
{
    //!!\brief Import base class constructor.
    using detail::strong_type<uint8_t,
                              error_count,
                              detail::strong_type_skill::convert | detail::strong_type_skill::comparable>::strong_type;
};

/*!\brief A strong type of underlying type `double` that represents the rate of errors.
 * \ingroup search_configuration
 * \see search_configuration
 * \tparam value_t The underlying type.
 */
struct error_rate :
    detail::strong_type<double, error_rate, detail::strong_type_skill::convert | detail::strong_type_skill::comparable>
{
    //!!\brief Import base class constructor.
    using detail::strong_type<double,
                              error_rate,
                              detail::strong_type_skill::convert | detail::strong_type_skill::comparable>::strong_type;
};

} // namespace seqan3::search_cfg
