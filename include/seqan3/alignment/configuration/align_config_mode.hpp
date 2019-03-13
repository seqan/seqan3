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
//!\brief A strong type to select the global alignment mode.
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

/*!\brief Helper variable to select the global alignment.
 * \ingroup configuration
 *
 * \details
 *
 * ### Example
 *
 * \snippet snippet/alignment/configuration/align_cfg_mode_example.cpp global
 */
inline constexpr detail::global_alignment_type global_alignment;

/*!\brief Sets the alignment mode.
 * \ingroup configuration
 * \tparam mode_type The type of the alignment mode.
 *
 * \details
 *
 * The alignment algorithm can be categorised in different modes. For example, the local and the
 * \ref align_cfg::global_alignment "global" alignment are two different modes, while the semi-global alignment
 * is a variation of the global alignment. This differentiation makes it possible to define a subset of configurations
 * that can work with a particular mode. Since it is not possible to guess what the desired mode for a user is, this
 * configuration must be provided for the alignment algorithm and cannot be defaulted.
 *
 * ### Example
 *
 * \include test/snippet/alignment/configuration/minimal_alignment_config.cpp
 */
template <typename mode_type>
//!\cond
    requires std::Same<remove_cvref_t<mode_type>, detail::global_alignment_type>
//!\endcond
struct mode : public pipeable_config_element<mode<mode_type>, mode_type>
{
    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::align_config_id id{mode_type::id};
};

/*!\name Type deduction guides
 * \relates seqan3::align_cfg::mode
 * \{
 */
//!\brief Deduces the alignment mode from the given constructor argument.
template <typename mode_type>
mode(mode_type) -> mode<mode_type>;
//!}
} // namespace seqan3::align_cfg
