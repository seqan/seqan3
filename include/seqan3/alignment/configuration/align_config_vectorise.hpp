// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::align_cfg::vectorise configuration.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 * \author Lydia Buntrock <lydia.buntrock AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/detail.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>
#include <seqan3/core/detail/empty_type.hpp>

namespace seqan3::detail
{

/*!\brief A global configuration type used to enable vectorisation of the alignment algorithm.
 * \ingroup alignment_configuration
 */
struct vectorise_tag : public pipeable_config_element<vectorise_tag, empty_type>
{
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::align_config_id id{detail::align_config_id::vectorise};
};

} // namespace seqan3::detail

namespace seqan3::align_cfg
{

/*!\brief Enables that the alignment computation is performed with SIMD vectors.
 * \ingroup alignment_configuration
 *
 * \details
 *
 * SIMD vectors allow the computation of several pairwise alignments simultaneously in one process.
 * Depending on your processor architecture you can gain a significant speed-up.
 *
 * \sa For further information on SIMD see https://en.wikipedia.org/wiki/SIMD.
 *
 * ### Example
 *
 * \include test/snippet/alignment/configuration/align_cfg_vectorise_example.cpp
 */
inline constexpr detail::vectorise_tag vectorise{};

} // namespace seqan3::align_cfg
