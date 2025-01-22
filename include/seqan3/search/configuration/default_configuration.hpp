// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides the default configuration for the seqan3::search() interface.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \author Lydia Buntrock <lydia.buntrock AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/search/configuration/detail.hpp>
#include <seqan3/search/configuration/hit.hpp>
#include <seqan3/search/configuration/max_error.hpp>
#include <seqan3/search/configuration/output.hpp>

namespace seqan3::search_cfg
{

/*!\brief The default configuration: Compute all exact matches.
 * \ingroup search_configuration
 * \see search_configuration
 */
constexpr configuration default_configuration = max_error_total{error_count{0}} | max_error_substitution{error_count{0}}
                                              | max_error_insertion{error_count{0}} | max_error_deletion{error_count{0}}
                                              | output_query_id{} | output_reference_id{}
                                              | output_reference_begin_position{} | hit_all{};
} // namespace seqan3::search_cfg
