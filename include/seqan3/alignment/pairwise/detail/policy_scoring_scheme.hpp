// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::policy_scoring_scheme.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/align_config_scoring.hpp>

namespace seqan3::detail
{

/*!\brief Stores the configured scoring scheme used for this algorithm.
 * \ingroup pairwise_alignment
 *
 * \tparam alignment_configuration_t The type of the alignment configuration; must be a type specialisation of
 *                                   seqan3::configuration.
 * \tparam scoring_scheme_t The type of the scoring scheme.
 *
 * \details
 *
 * Stores and initialises the configured scoring scheme from the given alignment configuration.
 */
template <typename alignment_configuration_t, typename scoring_scheme_t>
class policy_scoring_scheme
{
protected:
    //!\brief The scoring scheme used for this alignment algorithm.
    scoring_scheme_t m_scoring_scheme{};

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr policy_scoring_scheme() noexcept = default; //!< Defaulted.
    constexpr policy_scoring_scheme(policy_scoring_scheme const &) noexcept = default; //!< Defaulted.
    constexpr policy_scoring_scheme(policy_scoring_scheme &&) noexcept = default; //!< Defaulted.
    constexpr policy_scoring_scheme & operator=(policy_scoring_scheme const &) noexcept = default; //!< Defaulted.
    constexpr policy_scoring_scheme & operator=(policy_scoring_scheme &&) noexcept = default; //!< Defaulted.
    ~policy_scoring_scheme() noexcept = default; //!< Defaulted.

    /*!\brief Construction and initialisation using the alignment configuration.
     * \param[in] config The alignment configuration with the stored scoring scheme.
     */
    policy_scoring_scheme(alignment_configuration_t const & config) :
        m_scoring_scheme{seqan3::get<align_cfg::scoring>(config).value}
    {}
    //!\}
};

} // namespace seqan3::detail
