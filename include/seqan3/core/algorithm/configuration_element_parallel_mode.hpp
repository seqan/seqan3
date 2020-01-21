// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::parallel_mode.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/algorithm/pipeable_config_element.hpp>

namespace seqan3::detail
{
/*!\brief A global configuration type used to enable parallel execution of algorithms.
 * \ingroup algorithm
 * \tparam wrapped_config_id_t The algorithm specific configuration id wrapped in a std::integral_constant.
 *
 * \details
 *
 * This type is used to enable the parallel mode of the algorithms.
 */
template <typename wrapped_config_id_t>
struct parallel_mode : public pipeable_config_element<parallel_mode<wrapped_config_id_t>, uint32_t>
{
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr typename wrapped_config_id_t::value_type id{wrapped_config_id_t::value};
};
} // namespace seqan3::detail
