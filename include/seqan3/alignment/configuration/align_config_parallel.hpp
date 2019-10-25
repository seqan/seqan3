// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::align_cfg::parallel configuration.
 * \author Lydia Buntrock <lydia.buntrock AT fu-berlin.de>
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/detail.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>

namespace seqan3::align_cfg
{
/*!\brief Sets a flag whether the alignment computation should be performed in parallel.
 * \ingroup alignment_configuration
 *
 * \details
 *
 * Parallel computing allows many calculations or the execution of processes to be performed simultaneously.
 * Depending on your processor architecture you can gain a significant speed-up.
 *
 * The value represents the number of threads to be used. The default value (0) uses all available cores of your CPU.
 *
 * ### Example
 *
 * \include test/snippet/alignment/configuration/align_cfg_parallel_example.cpp
 */
struct parallel : public pipeable_config_element<parallel, uint32_t>
{
    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::align_config_id id{detail::align_config_id::parallel};
};

} // namespace seqan3::align_cfg
