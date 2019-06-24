// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides configuration for alignment output.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/detail.hpp>
#include <seqan3/alignment/pairwise/alignment_result.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>

namespace seqan3::detail
{

/*!\brief Triggers score-only computation of the sequence alignment.
 * \ingroup alignment_configuration
 */
struct with_score_type
{};

/*!\brief Triggers score computation and determines the end position of the sequence alignment.
 * \ingroup alignment_configuration
 */
struct with_back_coordinate_type
{};

/*!\brief Triggers score computation and determines begin and end position of the sequence alignment.
 * \ingroup alignment_configuration
 */
struct with_front_coordinate_type
{};

/*!\brief Triggers score computation and determines the end position of the sequence alignment as well as the
 *        full traceback.
 * \ingroup alignment_configuration
 */
struct with_alignment_type
{};

} // namespace seqan3::detail

namespace seqan3
{

//!\brief Helper variable used to select score-only computation.
//!\relates seqan3::align_cfg::result
inline constexpr detail::with_score_type with_score{};
//!\brief Helper variable used to select end-position computation.
//!\relates seqan3::align_cfg::result
inline constexpr detail::with_back_coordinate_type with_back_coordinate{};
//!\brief Helper variable used to select begin position computation.
//!\relates seqan3::align_cfg::result
inline constexpr detail::with_front_coordinate_type with_front_coordinate{};
//!\brief Helper Variable used to select trace computation.
//!\relates seqan3::align_cfg::result
inline constexpr detail::with_alignment_type with_alignment{};

} // namespace seqan3

namespace seqan3::align_cfg
{

/*!\brief Sets the result of the alignment computation.
 * \ingroup alignment_configuration
 * \tparam with_type  The type used to specify which feature should be computed during the pairwise alignment.
 *                    Defaults to seqan3::detail::with_score_type.
 *
 * \details
 *
 * The output of the pairwise alignment can be configured using the result configuration element. Depending on the
 * settings, the most efficient implementation is chosen to compute the result. Currently four different modes
 * can be configured: computing only the \ref seqan3::align_cfg::result::with_score "score",
 * computing in addition the \ref seqan3::align_cfg::result::with_back_coordinate "end position", computing in
 * addition the \ref seqan3::align_cfg::result::with_front_coordinate "begin position", and finally also
 * computing the \ref seqan3::align_cfg::result::with_alignment "alignment".
 * These settings will directly affect the contents of the seqan3::alignment_result object which is returned by the
 * alignment algorithm.
 *
 * ### Example
 *
 * \include test/snippet/alignment/configuration/align_cfg_result_example.cpp
 */
template <typename with_type = detail::with_score_type>
//!\cond
    requires std::Same<with_type, detail::with_score_type> ||
             std::Same<with_type, detail::with_back_coordinate_type> ||
             std::Same<with_type, detail::with_front_coordinate_type> ||
             std::Same<with_type, detail::with_alignment_type>
//!\endcond
class result : public pipeable_config_element<result<with_type>, with_type>
{
public:
    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::align_config_id id{detail::align_config_id::result};
};

/*!\name Type deduction guides
 * \relates seqan3::align_cfg::result
 * \{
 */
//!\brief The default constructor defaults to score-only computation.
result() -> result<detail::with_score_type>;

//!\brief Deduces the alignment result from the given constructor argument.
template <typename with_type>
result(with_type) -> result<remove_cvref_t<with_type>>;
//!\}
} // namespace seqan3::align_cfg
