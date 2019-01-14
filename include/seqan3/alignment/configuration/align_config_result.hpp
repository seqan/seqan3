// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides configuration for alignment output.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/detail.hpp>
#include <seqan3/alignment/pairwise/align_result.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>

namespace seqan3::detail
{

/*!\brief Triggers score-only computation of the sequence alignment.
 * \ingroup configuration
 */
struct with_score_type
{};

/*!\brief Triggers score computation and determines the end position of the sequence alignment.
 * \ingroup configuration
 */
struct with_end_position_type
{};

/*!\brief Triggers score computation and determines begin and end position of the sequence alignment.
 * \ingroup configuration
 */
struct with_begin_position_type
{};

/*!\brief Triggers score computation and determines the end position of the sequence alignment as well as the
 *        full traceback.
 * \ingroup configuration
 */
struct with_trace_type
{};

} // namespace seqan3

namespace seqan3::align_cfg
{

//!\brief A tag used to select score-only computation.
//!\relates seqan3::detail::with_score_type;
inline constexpr detail::with_score_type with_score{};
//!\brief A tag used to select end-position computation.
//!\relates seqan3::detail::with_end_position_type;
inline constexpr detail::with_end_position_type with_end_position{};
//!\brief A tag used to select begin position computation.
//!\relates seqan3::detail::with_begin_position_type;
inline constexpr detail::with_begin_position_type with_begin_position{};
//!\brief A tag used to select trace computation.
//!\relates seqan3::detail::with_trace_type;
inline constexpr detail::with_trace_type with_trace{};


/*!\brief A configuration element for specifying the result of the alignment.
 * \ingroup configuration
 * \tparam with_tag_t The tag used to specify which feature should be computed during the pairwise alignment.
 */
template <typename with_tag_t = detail::with_score_type>
//!\cond
    requires std::Same<with_tag_t, detail::with_score_type> || std::Same<with_tag_t, detail::with_end_position_type> ||
             std::Same<with_tag_t, detail::with_begin_position_type> || std::Same<with_tag_t, detail::with_trace_type>
//!\endcond
class result : public pipeable_config_element<result<with_tag_t>, with_tag_t>
{
public:
    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::align_config_id id{detail::align_config_id::result};
};

/*!\name Type deduction guides
 * \relates seqan3::align_cfg::mode
 * \{
 */
//!\brief The default constructor defaults to score-only computation.
result() -> result<detail::with_score_type>;

//!\brief Deduces the alignment result from the given constructor argument.
template <typename with_tag_t>
result(with_tag_t) -> result<remove_cvref_t<with_tag_t>>;
//!\}
} //namespace seqan3::align_cfg
