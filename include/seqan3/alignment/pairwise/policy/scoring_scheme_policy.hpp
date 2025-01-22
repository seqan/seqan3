// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::detail::scoring_scheme_policy.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/utility/simd/concept.hpp>
#include <seqan3/utility/type_traits/basic.hpp>

namespace seqan3::detail
{

/*!\brief The CRTP-policy that stores the scoring scheme used for this alignment algorithm.
 * \ingroup alignment_pairwise_policy
 * \tparam alignment_algorithm_t The derived type (seqan3::detail::alignment_algorithm) to be augmented with this
 *                               CRTP-policy.
 * \tparam scoring_scheme_t The type of the scoring scheme.
 *
 * \remarks The template parameters of this CRTP-policy are selected in the
 *          seqan3::detail::alignment_configurator::configure_scoring_scheme when selecting the alignment for the given
 *          configuration.
 */
template <typename alignment_algorithm_t, typename scoring_scheme_t>
class scoring_scheme_policy
{
private:
    //!\brief Befriends the derived class to grant it access to the private members.
    friend alignment_algorithm_t;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Defaulted.
    constexpr scoring_scheme_policy() = default;
    //!\brief Defaulted.
    constexpr scoring_scheme_policy(scoring_scheme_policy const &) = default;
    //!\brief Defaulted.
    constexpr scoring_scheme_policy(scoring_scheme_policy &&) = default;
    //!\brief Defaulted.
    constexpr scoring_scheme_policy & operator=(scoring_scheme_policy const &) = default;
    //!\brief Defaulted.
    constexpr scoring_scheme_policy & operator=(scoring_scheme_policy &&) = default;
    //!\brief Defaulted.
    ~scoring_scheme_policy() = default;

    //!\brief Initialise the policy.
    template <typename configuration_t>
    scoring_scheme_policy(configuration_t const & /*config*/)
    {}
    //!\}

    //!\brief The scoring scheme used for this alignment algorithm.
    scoring_scheme_t scoring_scheme{};

    /*!\brief Maybe converts the given sequence value to a specific profile used by the underlying scoring scheme.
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
