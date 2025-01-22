// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::align_config::gap_cost_affine.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \author Wiep van der Toorn <w.vandertoorn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/detail.hpp>
#include <seqan3/core/configuration/pipeable_config_element.hpp>
#include <seqan3/core/detail/strong_type.hpp>

namespace seqan3::align_cfg
{
// ------------------------------------------------------------------
// seqan3::align_cfg::open_score
// ------------------------------------------------------------------

/*!\brief A strong type of underlying type `int32_t` that represents a score (usually negative) that
 *        is incurred once per stretch of consecutive gaps.
 * \ingroup alignment_configuration
 * \see seqan3::align_cfg::gap_cost_affine
 */
struct open_score : seqan3::detail::strong_type<int32_t, open_score, seqan3::detail::strong_type_skill::convert>
{
    //!\brief The type of the strong type base class.
    using base_t = seqan3::detail::strong_type<int32_t, open_score, seqan3::detail::strong_type_skill::convert>;
    using base_t::base_t; // Import the base class constructors
};

// ------------------------------------------------------------------
// seqan3::align_cfg::extension_score
// ------------------------------------------------------------------

/*!\brief A strong type of underlying type `int32_t` that represents the score (usually negative) of any character
  *        against a gap character.
  * \ingroup alignment_configuration
  * \see seqan3::align_cfg::gap_cost_affine
  */
struct extension_score :
    seqan3::detail::strong_type<int32_t, extension_score, seqan3::detail::strong_type_skill::convert>
{
    //!\brief The type of the strong type base class.
    using base_t = seqan3::detail::strong_type<int32_t, extension_score, seqan3::detail::strong_type_skill::convert>;
    using base_t::base_t; // Import the base class constructors
};

/*!\brief A configuration element for the affine gap cost scheme.
 * \ingroup alignment_configuration
 *
 * \details
 *
 * Configures the gap scheme for the alignment algorithm. The gap scheme determines how gaps are penalised inside
 * of the alignment algorithm. If the gap scheme is not configured, it will default to a linear gap scheme initialised
 * with edit distance. Note that the \ref seqan3::align_cfg::open_score "gap open score" is used as an additional score.
 * This means that the score for opening a gap during the affine alignment execution is the sum of the gap score and the
 * gap open score.
 *
 * \see seqan3::align_cfg::method_global
 * ### Example
 *
 * \include test/snippet/alignment/configuration/align_cfg_gap_cost_affine_example.cpp
 *
 */
class gap_cost_affine : private pipeable_config_element
{
public:
    //!\brief The score per consecutive sequence of gaps. Defaults to 0.
    int32_t open_score{0};
    //!\brief The cost per gap character. Defaults to -1.
    int32_t extension_score{-1};

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr gap_cost_affine() = default;                                    //!< Defaulted
    constexpr gap_cost_affine(gap_cost_affine const &) = default;             //!< Defaulted
    constexpr gap_cost_affine(gap_cost_affine &&) = default;                  //!< Defaulted
    constexpr gap_cost_affine & operator=(gap_cost_affine const &) = default; //!< Defaulted
    constexpr gap_cost_affine & operator=(gap_cost_affine &&) = default;      //!< Defaulted
    ~gap_cost_affine() = default;                                             //!< Defaulted

    /*!\brief Construction from strongly typed open score and extension score.
     * \param open_score The cost per consecutive sequence of gaps (of type seqan3::open_score).
     * \param extension_score The cost of each gap character (of type seqan3::extension_score).
     *
     * \details
     *
     * The score for a sequence of `n` gap characters is computed as `open_score + n * extension_score`.
     *
     * \attention This is the formula used most commonly in the literature, but it is different from SeqAn2 where the
     * formula was `(n-1) * extension_score + open_score`.
     */
    constexpr gap_cost_affine(seqan3::align_cfg::open_score open_score,
                              seqan3::align_cfg::extension_score extension_score) :
        open_score(std::move(open_score)),
        extension_score(std::move(extension_score))
    {}
    //!\}

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr seqan3::detail::align_config_id id{seqan3::detail::align_config_id::gap};
};

} // namespace seqan3::align_cfg
