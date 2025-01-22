// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::align_cfg::detail::result_type.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <seqan3/alignment/configuration/detail.hpp>
#include <seqan3/alignment/pairwise/alignment_result.hpp>
#include <seqan3/core/configuration/pipeable_config_element.hpp>

namespace seqan3::align_cfg::detail
{
/*!\brief Configuration element capturing the configured seqan3::alignment_result for the alignment algorithm.
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
 */
template <typename alignment_result_t>
    requires seqan3::detail::is_type_specialisation_of_v<alignment_result_t, seqan3::alignment_result>
class result_type : private pipeable_config_element
{
public:
    //!\brief The result type.
    using type = alignment_result_t;

    /*!\name Constructor, destructor and assignment
     * \{
     */
    constexpr result_type() = default;                                //!< Defaulted.
    constexpr result_type(result_type const &) = default;             //!< Defaulted.
    constexpr result_type(result_type &&) = default;                  //!< Defaulted.
    constexpr result_type & operator=(result_type const &) = default; //!< Defaulted.
    constexpr result_type & operator=(result_type &&) = default;      //!< Defaulted.
    ~result_type() = default;                                         //!< Defaulted.

    //!\}

    //!\brief Internal id to check for consistent configuration settings.
    static constexpr seqan3::detail::align_config_id id{seqan3::detail::align_config_id::result_type};
};
} // namespace seqan3::align_cfg::detail
