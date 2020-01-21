// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::align_cfg::max_error configuration.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/detail.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>

namespace seqan3::align_cfg
{
/*!\brief Sets the maximal errors allowed during an edit distance computation.
 * \ingroup alignment_configuration
 *
 * \details
 *
 * This configuration can only be used for computing the \ref seqan3::align_cfg::edit "edit distance".
 * It restricts the number of substitutions, insertions, and deletions within the alignment to the given value and
 * can thereby speed up the edit distance computation.
 * A typical use case is to verify a candidate region during read mapping where the number of maximal errors is given
 * beforehand. If this configuration is used for an alignment algorithm that does not compute the edit distance, a
 * seqan3::invalid_alignment_configuration exception will be thrown.
 *
 * ### Example
 *
 * \include test/snippet/alignment/configuration/align_cfg_max_error_example.cpp
 */
struct max_error : public pipeable_config_element<max_error, uint32_t>
{
    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::align_config_id id{detail::align_config_id::max_error};
};

} // namespace seqan3::align_cfg
