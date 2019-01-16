// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides global alignment configurations.
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/detail.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>

namespace seqan3::detail
{
//!\brief Selects the global alignment mode.
//!\ingroup configuration
struct global_alignment_type
{
    //!\privatesection
    //!\brief An internal id used to check for a valid alignment configuration.
    static constexpr detail::align_config_id id{detail::align_config_id::global};
};

} // namespace seqan3::detail

namespace seqan3::align_cfg
{

//!\brief Selects global alignment mode.
//!\ingroup configuration
inline constexpr detail::global_alignment_type global_alignment;

/*!\brief A configuration element for global alignment.
 * \ingroup configuration
 */
template <typename alignment_t>
//!\cond
    requires std::Same<remove_cvref_t<alignment_t>, detail::global_alignment_type>
//!\endcond
struct mode : public pipeable_config_element<mode<alignment_t>, alignment_t>
{
    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::align_config_id id{alignment_t::id};
};

/*!\name Type deduction guides
 * \relates seqan3::align_cfg::mode
 * \{
 */
//!\brief Deduces the alignment mode from the given constructor argument.
template <typename alignment_t>
mode(alignment_t) -> mode<alignment_t>;
//!}
} // namespace seqan3::align_cfg
