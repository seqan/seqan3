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

#include <seqan3/alignment/configuration/align_config_max_error.hpp>
#include <seqan3/alignment/configuration/align_config_result.hpp>
#include <seqan3/alignment/pairwise/align_result_selector.hpp>
#include <seqan3/core/bit_manipulation.hpp>
#include <seqan3/core/platform.hpp>
#include <seqan3/core/type_traits/deferred_crtp_base.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

namespace seqan3::detail
{
template <typename word_t, typename score_t, bool is_semi_global, bool use_max_errors>
class edit_distance_score_matrix_full; //forward declaration

template <typename word_t, bool is_semi_global, bool use_max_errors>
class edit_distance_trace_matrix_full; //forward declaration

//!\brief Store no state for state_t.
template <typename state_t, typename...>
struct empty_state
{};

//!\brief If enabled is true state_t will be added to state_type.
template <bool enabled, typename state_t>
using enable_state_t = std::conditional_t<enabled, state_t, empty_state<state_t>>;

/*!\brief The same as std::conditional but for template template parameters.
 * \details
 * If `B` is true, <tt>selector<B, T, F>::template select</tt> inherits `T`, otherwise `F`.
 */
template<bool B, template<typename...> typename T, template<typename...> typename F>
struct selector
{
    //!\brief Depending on `B`, `select` is the template template parameter `T` or `F`.
    template <typename ...args_t>
    struct select : public std::conditional_t<B, T<args_t...>, F<args_t...>>
    {};
};

/*!\brief The default traits type for the edit distance algorithm.
 * \ingroup pairwise_alignment
 */
template <std::ranges::ViewableRange database_t,
          std::ranges::ViewableRange query_t,
          typename align_config_t,
          typename is_semi_global_t,
          typename word_t = uint_fast64_t>
struct default_edit_distance_trait_type
{
    //!\brief The type of one machine word.
    using word_type = word_t;
    static_assert(std::is_unsigned_v<word_type>, "the word type of edit_distance_unbanded must be unsigned.");
    //!\brief The type of the score.
    using score_type = int;
    //!\brief The type of the database sequence.
    using database_type = std::remove_reference_t<database_t>;
    //!\brief The type of the query sequence.
    using query_type = std::remove_reference_t<query_t>;
    //!\brief The type of the alignment config.
    using align_config_type = std::remove_reference_t<align_config_t>;

    //!\brief The size of one machine word.
    static constexpr uint8_t word_size = sizeof_bits<word_type>;
    static_assert(sizeof_bits<word_type> <= 64u, "we assume at most uint64_t as word_type");

    //!\brief The type of an iterator of the database sequence.
    using database_iterator = std::ranges::iterator_t<database_type>;
    //!\brief The alphabet type of the query sequence.
    using query_alphabet_type = std::remove_reference_t<reference_t<query_type>>;
    //!\brief The intermediate result type of the execution of this function object.
    using result_value_type = typename align_result_selector<database_type, query_type, align_config_type>::type;

    //!\brief When true the computation will use the ukkonen trick with the last active cell and bounds the error to
    //!       config.max_errors.
    static constexpr bool use_max_errors = align_config_type::template exists<align_cfg::max_error>();
    //!\brief Whether the alignment is a semi-global alignment or not.
    static constexpr bool is_semi_global = is_semi_global_t::value;
    //!\brief Whether the alignment is a global alignment or not.
    static constexpr bool is_global = !is_semi_global;

    //!\brief Whether the alignment configuration indicates to compute and/or store the score.
    static constexpr bool compute_score = align_config_type::template exists<align_cfg::result<with_score_type>>() ||
                                          !std::Same<decltype(result_value_type{}.back_coordinate), std::nullopt_t *>;
    //!\brief Whether the alignment configuration indicates to compute and/or store the back coordinate.
    static constexpr bool compute_back_coordinate = !std::Same<decltype(result_value_type{}.back_coordinate),
                                                               std::nullopt_t *>;
    //!\brief Whether the alignment configuration indicates to compute and/or store the front coordinate.
    static constexpr bool compute_front_coordinate = !std::Same<decltype(result_value_type{}.front_coordinate),
                                                                std::nullopt_t *>;
    //!\brief Whether the alignment configuration indicates to compute and/or store the alignment of the sequences.
    static constexpr bool compute_sequence_alignment = !std::Same<decltype(result_value_type{}.alignment),
                                                                  std::nullopt_t *>;
    //!\brief Whether the alignment configuration indicates to compute and/or store the score matrix.
    static constexpr bool compute_score_matrix = false;
    //!\brief Whether the alignment configuration indicates to compute and/or store the trace matrix.
    static constexpr bool compute_trace_matrix = compute_front_coordinate || compute_sequence_alignment;
    //!\brief Whether the alignment configuration indicates to compute and/or store the score or trace matrix.
    static constexpr bool compute_matrix = compute_score_matrix || compute_trace_matrix;

    //!\brief The type of the trace matrix.
    using trace_matrix_type = edit_distance_trace_matrix_full<word_type, is_semi_global, use_max_errors>;
    //!\brief The type of the score matrix.
    using score_matrix_type = edit_distance_score_matrix_full<word_type, score_type, is_semi_global, use_max_errors>;
};

//!\brief A base class for edit_distance_unbanded.
template <bool enable_policy,
          template<typename...> typename policy_t,
          typename edit_traits,
          typename derived_t>
using edit_distance_base =
    invoke_deferred_crtp_base<deferred_crtp_base<selector<enable_policy, policy_t, empty_state>::template select,
                                                 edit_traits>,
                              derived_t>;

//!\cond
template <std::ranges::ViewableRange database_t,
          std::ranges::ViewableRange query_t,
          typename align_config_t,
          typename traits_t = default_edit_distance_trait_type<database_t, query_t, align_config_t, std::false_type>>
class edit_distance_unbanded; //forward declaration
//!\endcond

} // namespace seqan3::detail
