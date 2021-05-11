// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides alignment configuration seqan3::align_cfg::score_type.
 * \author Lydia Buntrock <lydia.buntrock AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/detail.hpp>
#include <seqan3/core/configuration/pipeable_config_element.hpp>
#include <seqan3/utility/detail/exposition_only_concept.hpp>

namespace seqan3::align_cfg
{
/*!\brief A configuration element to set the score type used in the alignment algorithm.
 * \ingroup alignment_configuration
 * \tparam score_t The type to use for the computed alignment score; must model seqan3::arithmetic.
 *
 * \details
 *
 * This option configures the score type of the alignment algorithm.
 * By default, the alignment algorithm will only compute the score with score type `int32_t`.
 *
 * ### Example
 *
 * \include test/snippet/alignment/configuration/align_cfg_score_type.cpp
 */
template <arithmetic score_t>
class score_type : private pipeable_config_element
{
public:

    static_assert(std::floating_point<score_t> || std::signed_integral<score_t>,
                  "The selected score type must be a signed integral type or floating point type.");

    //!\brief The selected score type.
    using type = score_t;

    /*!\name Constructor, destructor and assignment
     * \{
     */
    constexpr score_type() = default; //!< Defaulted.
    constexpr score_type(score_type const &) = default; //!< Defaulted.
    constexpr score_type(score_type &&) = default; //!< Defaulted.
    constexpr score_type & operator=(score_type const &) = default; //!< Defaulted.
    constexpr score_type & operator=(score_type &&) = default; //!< Defaulted.
    ~score_type() = default; //!< Defaulted.

    //!\}

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr seqan3::detail::align_config_id id{seqan3::detail::align_config_id::score_type};
};

} // namespace seqan3::align_cfg
