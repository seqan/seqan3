// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides the default configuration for the seqan3::search() interface.
 */

#pragma once

#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/search/configuration/detail.hpp>
#include <seqan3/search/configuration/max_error.hpp>
#include <seqan3/search/configuration/max_error_rate.hpp>
#include <seqan3/search/configuration/mode.hpp>
#include <seqan3/search/configuration/output.hpp>

namespace seqan3::search_cfg
{

/*!\brief The default configuration.
 * \ingroup search_configuration
 */
inline constexpr configuration default_configuration = max_error{total{0}, substitution{0}, insertion{0}, deletion{0}} |
                                                       output{text_position} |
                                                       mode{all};

} // namespace seqan3::search_cfg
