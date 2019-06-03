// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Forwards for seqan3::edit_distance_unbanded related types.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
*/

#pragma once

#include <seqan3/core/platform.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

namespace seqan3::detail
{
/*!\todo Document me
 * \ingroup pairwise_alignment
 */
template <typename traits_type>
SEQAN3_CONCEPT EditDistanceTrait = requires
{
    typename std::remove_reference_t<traits_type>::word_type;
    typename std::remove_reference_t<traits_type>::is_semi_global_type;

    // Must be a boolean integral constant.
    requires std::Same<typename std::remove_reference_t<traits_type>::is_semi_global_type::value_type, bool>;
};

/*!\brief The default traits type for the edit distance algorithm.
 * \ingroup pairwise_alignment
 */
struct default_edit_distance_trait_type
{
    //!\brief The default word type.
    using word_type = uint64_t;
    //!\brief Semi global alignment is disabled by default.
    using is_semi_global_type = std::false_type;
};

//!\brief Store no state for state_t.
template <typename state_t>
struct empty_state
{};

//!\brief If enabled is true state_t will be added to state_type.
template <bool enabled, typename state_t>
using enable_state_t = std::conditional_t<enabled, state_t, empty_state<state_t>>;

//!\cond
template <std::ranges::ViewableRange database_t,
          std::ranges::ViewableRange query_t,
          typename align_config_t,
          EditDistanceTrait traits_t = default_edit_distance_trait_type>
class pairwise_alignment_edit_distance_unbanded; //forward declaration

template <typename word_t, typename score_t, bool is_semi_global, bool use_max_errors>
class edit_distance_score_matrix_full; //forward declaration

template <typename word_t, bool is_semi_global, bool use_max_errors>
class edit_distance_trace_matrix_full; //forward declaration
//!\endcond

} // namespace seqan3::detail
