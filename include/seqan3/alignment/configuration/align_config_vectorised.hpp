// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::align_cfg::vectorised configuration.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 * \author Lydia Buntrock <lydia.buntrock AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/detail.hpp>
#include <seqan3/core/configuration/pipeable_config_element.hpp>
#include <seqan3/core/detail/empty_type.hpp>

namespace seqan3::align_cfg
{

/*!\brief Enables the vectorised alignment computation if possible for the current configuration.
 * \ingroup alignment_configuration
 *
 * \details
 *
 * In the vectorised alignment computation several pairwise sequence alignments are processed simultaneously in one
 * invocation. To do so, we pack the alignments in so called extended SIMD registers which allow to compute a single
 * instruction on multiple data at the same time. Depending on your processor architecture you can gain a significant
 * speed-up, e.g. by running up to 64 alignments in parallel on the latest intel CPUs. In our mode we vectorise
 * multiple alignments and not a single alignment. This means that you should provide many sequences to compute as
 * one batch rather than computing them separately as there won't be performance gains.
 *
 * \sa For further information on SIMD see https://en.wikipedia.org/wiki/SIMD.
 *
 * ### Example
 *
 * \include test/snippet/alignment/configuration/align_cfg_vectorised_example.cpp
 */
class vectorised : private pipeable_config_element
{
public:
    /*!\name Constructor, destructor and assignment
     * \{
     */
    constexpr vectorised() = default;                               //!< Defaulted.
    constexpr vectorised(vectorised const &) = default;             //!< Defaulted.
    constexpr vectorised(vectorised &&) = default;                  //!< Defaulted.
    constexpr vectorised & operator=(vectorised const &) = default; //!< Defaulted.
    constexpr vectorised & operator=(vectorised &&) = default;      //!< Defaulted.
    ~vectorised() = default;                                        //!< Defaulted.

    //!\}

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr seqan3::detail::align_config_id id{seqan3::detail::align_config_id::vectorised};
};

} // namespace seqan3::align_cfg
