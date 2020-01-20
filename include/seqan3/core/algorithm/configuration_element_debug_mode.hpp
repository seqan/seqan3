// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::debug_mode.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/algorithm/pipeable_config_element.hpp>

namespace seqan3::detail
{
/*!\brief A global configuration type used to enabled debugging of algorithms.
 * \ingroup algorithm
 * \tparam wrapped_config_id_t The algorithm specific configuration id wrapped in a std::integral_constant.
 *
 * \details
 *
 * This type is used to enable specific debugging behaviour of the algorithms, e.g. to output the score and the
 * trace matrix of the alignment algorithm.
 */
template <typename wrapped_config_id_t>
struct debug_mode : public pipeable_config_element<debug_mode<wrapped_config_id_t>, bool>
{
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr typename wrapped_config_id_t::value_type id{wrapped_config_id_t::value};
};
} // namespace seqan3::detail
