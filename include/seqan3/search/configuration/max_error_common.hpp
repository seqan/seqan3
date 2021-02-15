// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the error types for maximum number of errors.
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \author Lydia Buntrock <lydia.buntrock AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/concepts>

#include <seqan3/core/detail/strong_type.hpp>
#include <seqan3/utility/detail/exposition_only_concept.hpp>

namespace seqan3::search_cfg
{

/*!\brief A strong type of underlying type `uint8_t` that represents the number of errors.
 * \ingroup search_configuration
 * \tparam value_t The underlying type.
 */
struct error_count : detail::strong_type<uint8_t,
                                         error_count,
                                         detail::strong_type_skill::convert | detail::strong_type_skill::comparable>
{
    //!!\brief Import base class constructor.
    using detail::strong_type<uint8_t,
                              error_count,
                              detail::strong_type_skill::convert | detail::strong_type_skill::comparable>::strong_type;
};

/*!\brief A strong type of underlying type `double` that represents the rate of errors.
 * \ingroup search_configuration
 * \tparam value_t The underlying type.
 */
struct error_rate : detail::strong_type<double,
                                        error_rate,
                                        detail::strong_type_skill::convert | detail::strong_type_skill::comparable>
{
    //!!\brief Import base class constructor.
    using detail::strong_type<double,
                              error_rate,
                              detail::strong_type_skill::convert | detail::strong_type_skill::comparable>::strong_type;
};

} // namespace seqan3::search_cfg
