// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::policy_affine_gap_recursion_simd.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <tuple>

#include <seqan3/alignment/pairwise/detail/policy_affine_gap_recursion.hpp>
#include <seqan3/core/simd/simd_algorithm.hpp>

namespace seqan3::detail
{

/*!\brief Implements the alignment recursion function for the vectorised alignment algorithm using affine gap costs.
 * \ingroup pairwise_alignment
 *
 * \copydetails seqan3::detail::policy_affine_gap_recursion
 */
template <typename alignment_configuration_t>
class policy_affine_gap_recursion_simd : protected policy_affine_gap_recursion<alignment_configuration_t>
{
protected:
    //!\brief The type of the base class.
    using base_policy_t = policy_affine_gap_recursion<alignment_configuration_t>;
    // Import the score type from the base.
    using typename base_policy_t::score_type;

    // Import parameters.
    using base_policy_t::gap_extension_score;
    using base_policy_t::gap_open_score;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    policy_affine_gap_recursion_simd() = default; //!< Defaulted.
    policy_affine_gap_recursion_simd(policy_affine_gap_recursion_simd const &) = default; //!< Defaulted.
    policy_affine_gap_recursion_simd(policy_affine_gap_recursion_simd &&) = default; //!< Defaulted.
    policy_affine_gap_recursion_simd & operator=(policy_affine_gap_recursion_simd const &) = default; //!< Defaulted.
    policy_affine_gap_recursion_simd & operator=(policy_affine_gap_recursion_simd &&) = default; //!< Defaulted.
    ~policy_affine_gap_recursion_simd() = default; //!< Defaulted.

    //!\copydoc seqan3::detail::policy_affine_gap_recursion::policy_affine_gap_recursion
    explicit policy_affine_gap_recursion_simd(alignment_configuration_t const & config)
    {
        // Get the gap scheme from the config or choose -1 and -10 as default.
        auto && selected_gap_scheme = config.template value_or<align_cfg::gap>(gap_scheme{gap_score{-1},
                                                                                          seqan3::gap_open_score{-10}});

        gap_extension_score = simd::fill<score_type>(selected_gap_scheme.get_gap_score());
        gap_open_score = simd::fill<score_type>(selected_gap_scheme.get_gap_open_score()) + gap_extension_score;
    }
    //!\}
};
} // namespace seqan3::detail
