// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::align_cfg::scoring_scheme.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <concepts>

#include <seqan3/alignment/configuration/detail.hpp>
#include <seqan3/alignment/scoring/scoring_scheme_concept.hpp>
#include <seqan3/core/configuration/pipeable_config_element.hpp>
#include <seqan3/utility/type_traits/basic.hpp>

namespace seqan3::align_cfg
{

/*!\brief Sets the scoring scheme for the alignment algorithm.
 * \ingroup alignment_configuration
 * \tparam scoring_scheme_t The type of the scoring scheme; must model std::semiregular.
 *
 * \details
 *
 * The scoring scheme allows to specify how two symbols of an alphabet are scored inside of the alignment algorithm.
 * The scheme depends on the alphabet type of the passed sequences and must be chosen accordingly.
 * During the configuration of the pairwise alignment algorithm a static assert is triggered if the scoring scheme
 * is not compatible with the given alphabet types (see seqan3::scoring_scheme_for). Accordingly,
 * there is no default for this configuration since it depends on the sequences and it must be given as a minimal
 * configuration.
 *
 * ### Example
 *
 * \include test/snippet/alignment/configuration/minimal_alignment_config.cpp
 */
template <std::semiregular scoring_scheme_t>
class scoring_scheme : private pipeable_config_element
{
public:
    //!\brief The scoring scheme to be used in the alignment algorithm.
    scoring_scheme_t scheme{};

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr scoring_scheme() = default;                                   //!< Defaulted
    constexpr scoring_scheme(scoring_scheme const &) = default;             //!< Defaulted
    constexpr scoring_scheme(scoring_scheme &&) = default;                  //!< Defaulted
    constexpr scoring_scheme & operator=(scoring_scheme const &) = default; //!< Defaulted
    constexpr scoring_scheme & operator=(scoring_scheme &&) = default;      //!< Defaulted
    ~scoring_scheme() = default;                                            //!< Defaulted

    /*!\brief Initialises the scoring scheme config with the given scheme.
     * \param[in] scheme The scoring scheme to be used in the alignment algorithm.
     *
     * \details
     *
     * This config stores a copy of the provided scheme.
     */
    explicit constexpr scoring_scheme(scoring_scheme_t scheme) : scheme{std::move(scheme)}
    {}
    //!\}

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr seqan3::detail::align_config_id id{seqan3::detail::align_config_id::scoring};
};

/*!\name Type deduction guides
 * \relates seqan3::align_cfg::scoring_scheme
 * \{
 */

//!\brief Deduces the scoring scheme type from the constructor argument.
template <typename scheme_t>
scoring_scheme(scheme_t) -> scoring_scheme<std::remove_cvref_t<scheme_t>>;
//!\}

} // namespace seqan3::align_cfg
