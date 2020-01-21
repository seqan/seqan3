// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides global alignment configurations.
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/detail.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>

namespace seqan3::detail
{
//!\brief A strong type to select the global alignment mode.
//!\ingroup alignment_configuration
struct global_alignment_type
{
    //!\privatesection
    //!\brief An internal id used to check for a valid alignment configuration.
    static constexpr detail::align_config_id id{detail::align_config_id::global};
};

//!\brief A strong type to select the local alignment mode.
//!\ingroup alignment_configuration
struct local_alignment_type
{
    //!\privatesection
    //!\brief An internal id used to check for a valid alignment configuration.
    static constexpr detail::align_config_id id{detail::align_config_id::local};
};

} // namespace seqan3::detail

namespace seqan3
{

/*!\brief Helper variable to select the global alignment.
 * \ingroup alignment_configuration
 *
 * \details
 *
 * ### Example
 *
 * \include snippet/alignment/configuration/align_cfg_global_mode_example.cpp
 */
inline constexpr detail::global_alignment_type global_alignment;

/*!\brief Helper variable to select the local alignment.
 * \ingroup alignment_configuration
 *
 * \details
 *
 * ### Example
 *
 * \include snippet/alignment/configuration/align_cfg_local_mode_example.cpp
 */
inline constexpr detail::local_alignment_type local_alignment;

} // namespace seqan3

namespace seqan3::align_cfg
{

/*!\brief Sets the alignment mode.
 * \ingroup alignment_configuration
 * \tparam mode_type The type of the alignment mode.
 *
 * \details
 *
 * The alignment algorithm can be categorised in different modes. For example, the
 * \ref seqan3::local_alignment "local" and the
 * \ref seqan3::global_alignment "global" alignment are two different modes, while the semi-global alignment
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
    requires std::same_as<remove_cvref_t<mode_type>, detail::global_alignment_type> ||
             std::same_as<remove_cvref_t<mode_type>, detail::local_alignment_type>
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
