// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::search_cfg::itv_threshold.
 * \author Sven Bönigk <sven.boenigk AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/type_traits/basic.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>
#include <seqan3/core/algorithm/parameter_pack.hpp>
#include <seqan3/search/algorithm/configuration/detail.hpp>

namespace seqan3::search_cfg
{
    //TODO
/*!\brief A configuration element for the in text verification in the index search.
 * \ingroup search_configuration
 * \tparam gap_scheme_t The type of the underlying gap scheme; must be of type seqan3::gap_scheme.
 *
 * \details
 *
 * Configures the in text verification the search algorithm. The threshold and minimum determines at what point the in text verification should be used during the search with indeces.
 * The in text verification is used as soon as the search range on the index is smaller than given threshold and more than minimum step backtracking steps were done.
 * The in text verification is only used after more than minimum steps
 * If the itv threshold is not configured, a threshold of 10 is taken. The minimum step parameter is calculated by the following formula: log(text_length)/log(4) + 3.
 *
 * ### Example
 *
 * \include test/snippet/alignment/configuration/align_cfg_gap_example.cpp
 *
 * \see seqan3::gap_scheme
 */

struct itv_threshold : public pipeable_config_element<itv_threshold, std::pair<uint16_t, uint16_t>>
{
    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::search_config_id id{detail::search_config_id::itv_threshold};
};


} // namespace seqan3::search_cfg
