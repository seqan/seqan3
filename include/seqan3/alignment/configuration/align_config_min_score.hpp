// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::align_cfg::min_score configuration.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <limits>

#include <seqan3/alignment/configuration/detail.hpp>
#include <seqan3/core/configuration/pipeable_config_element.hpp>

namespace seqan3::align_cfg
{
/*!\brief Sets the minimal score (maximal errors) allowed during an distance computation e.g. edit distance.
 * \ingroup alignment_configuration
 *
 * \details
 *
 * This configuration can only be used for computing the \ref seqan3::align_cfg::edit_scheme "edit distance".
 * It restricts the number of substitutions, insertions, and deletions within the alignment to the given value and
 * can thereby speed up the edit distance computation.
 * A typical use case is to verify a candidate region during read mapping where the number of maximal errors is given
 * beforehand. If this configuration is used for an alignment algorithm that does not compute the edit distance, a
 * seqan3::invalid_alignment_configuration exception will be thrown.
 *
 * ### Example
 *
 * \include test/snippet/alignment/configuration/align_cfg_min_score_example.cpp
 */
class min_score : private pipeable_config_element
{
public:
    //!\brief Minimal score for the distance computation [default: -infinity].
    int32_t score{std::numeric_limits<int32_t>::lowest()};

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr min_score() noexcept = default;                              //!< Defaulted
    constexpr min_score(min_score const &) noexcept = default;             //!< Defaulted
    constexpr min_score(min_score &&) noexcept = default;                  //!< Defaulted
    constexpr min_score & operator=(min_score const &) noexcept = default; //!< Defaulted
    constexpr min_score & operator=(min_score &&) noexcept = default;      //!< Defaulted
    ~min_score() noexcept = default;                                       //!< Defaulted

    /*!\brief Initialises the minimal score.
     *
     * \param score \copybrief score
     */
    constexpr min_score(int32_t const score) : score{score}
    {}
    //!\}

    //!\brief Internal id to check for consistent configuration settings.
    static constexpr seqan3::detail::align_config_id id{seqan3::detail::align_config_id::min_score};
};

} // namespace seqan3::align_cfg
