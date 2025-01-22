// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::detail::policy_scoring_scheme.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/align_config_scoring_scheme.hpp>
#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/utility/simd/concept.hpp>

namespace seqan3::detail
{

/*!\brief Stores the configured scoring scheme used for this algorithm.
 * \ingroup alignment_pairwise
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
    scoring_scheme_t scoring_scheme{};

    /*!\name Constructors, destructor and assignment
     * \{
     */
    policy_scoring_scheme() = default;                                          //!< Defaulted.
    policy_scoring_scheme(policy_scoring_scheme const &) = default;             //!< Defaulted.
    policy_scoring_scheme(policy_scoring_scheme &&) = default;                  //!< Defaulted.
    policy_scoring_scheme & operator=(policy_scoring_scheme const &) = default; //!< Defaulted.
    policy_scoring_scheme & operator=(policy_scoring_scheme &&) = default;      //!< Defaulted.
    ~policy_scoring_scheme() = default;                                         //!< Defaulted.

    /*!\brief Construction and initialisation using the alignment configuration.
     * \param[in] config The alignment configuration with the stored scoring scheme.
     */
    explicit policy_scoring_scheme(alignment_configuration_t const & config) :
        scoring_scheme{seqan3::get<align_cfg::scoring_scheme>(config).scheme}
    {}
    //!\}

    /*!\brief Maybe converts the given sequence value to a specific profile used by the underlying scoring scheme.
     *
     * \tparam alphabet_t The type of the actual alphabet; must model either seqan3::simd::simd_concept or
     *                    seqan3::semialphabet.
     *
     * \param[in] alphabet The alphabet value to get a score profile for.
     *
     * \details
     *
     * In the vectorised alignment the scoring scheme might transform the sequence values of the first sequence into
     * a profile for a more efficient comparison of the sequence characters in simd mode.
     *
     * If the given sequence type models seqan3::semialphabet the function becomes a no-op function and returns the
     * unmodified value.
     */
    template <typename alphabet_t>
        requires simd_concept<std::remove_cvref_t<alphabet_t>>
    auto scoring_scheme_profile_column(alphabet_t && alphabet) const noexcept
    {
        return scoring_scheme.make_score_profile(std::forward<alphabet_t>(alphabet));
    }

    //!\overload
    template <semialphabet alphabet_t>
    alphabet_t scoring_scheme_profile_column(alphabet_t && alphabet) const noexcept
    {
        return std::forward<alphabet_t>(alphabet);
    }
};

} // namespace seqan3::detail
