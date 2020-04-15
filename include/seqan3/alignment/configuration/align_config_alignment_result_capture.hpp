// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::align_cfg::alignment_result_capture.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/detail.hpp>
#include <seqan3/alignment/pairwise/alignment_result.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>
#include <seqan3/std/type_traits>

namespace seqan3::detail
{

/*!\brief The configuration element storing the captured alignment result type.
 * \ingroup alignment_configuration
 *
 * \tparam alignment_result_t The type of the alignment result to capture; must be a type specialisation of
 *                            seqan3::alignment_result.
 *
 * \details
 *
 * Implementation of the alignment result capture config.
 *
 * \see seqan3::align_cfg::alignment_result_capture
 */
template <typename alignment_result_t>
//!\cond
    requires is_type_specialisation_of_v<alignment_result_t, alignment_result>
//!\endcond
struct alignment_result_capture_element :
    public pipeable_config_element<alignment_result_capture_element<alignment_result_t>,
                                   std::type_identity<alignment_result_t>>
{
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::align_config_id id{detail::align_config_id::alignment_result_capture};
};

} // namespace seqan3::detail

namespace seqan3::align_cfg
{
/*!\if DEV
 * \brief Configuration element capturing the configured seqan3::alignment_result for the alignment algorithm.
 * \ingroup alignment_configuration
 * \tparam alignment_result_t The alignment result type to capture; must be a type specialisation of
 *                            seqan3::alignment_result.
 *
 * \details
 *
 * This configuration element allows to capture the concrete seqan3::alignment_result type after configuring the
 * alignment algorithm with the seqan3::detail::alignment_configurator. The actual result type is wrapped in
 * std::type_identity to preserve the trivial type properties of the configuration element. Thus, on access the
 * actual type needs to be unwrapped using the member typedef `type` before it can be used.
 * The result type can be accessed via the seqan3::detail::alignment_configuration_traits over the corresponding
 * alignment configuration type.
 * If the captured alignment result wasn't added yet to the alignment configuration the corresponding
 * result type member will deduce to seqan3::detail::empty_type.
 *
 * \note This configuration element is only added internally during the alignment configuration and is not intended for
 *       public use.
 * \endif
 */
template <typename alignment_result_t>
//!\cond
    requires detail::is_type_specialisation_of_v<alignment_result_t, alignment_result>
//!\endcond
inline constexpr detail::alignment_result_capture_element<alignment_result_t> alignment_result_capture{};
}  // namespace seqan3::align_cfg

