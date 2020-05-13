// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::alignment_selector.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <functional>
#include <tuple>
#include <utility>
#include <vector>

#include <seqan3/alignment/configuration/align_config_alignment_result_capture.hpp>
#include <seqan3/alignment/matrix/detail/alignment_score_matrix_one_column.hpp>
#include <seqan3/alignment/matrix/detail/alignment_score_matrix_one_column_banded.hpp>
#include <seqan3/alignment/matrix/detail/alignment_trace_matrix_full.hpp>
#include <seqan3/alignment/matrix/detail/alignment_trace_matrix_full_banded.hpp>
#include <seqan3/alignment/pairwise/detail/policy_affine_gap_recursion.hpp>
#include <seqan3/alignment/pairwise/detail/policy_optimum_tracker.hpp>
#include <seqan3/alignment/pairwise/policy/affine_gap_policy.hpp>
#include <seqan3/alignment/pairwise/policy/affine_gap_init_policy.hpp>
#include <seqan3/alignment/pairwise/policy/alignment_matrix_policy.hpp>
#include <seqan3/alignment/pairwise/policy/find_optimum_policy.hpp>
#include <seqan3/alignment/pairwise/policy/scoring_scheme_policy.hpp>
#include <seqan3/alignment/pairwise/policy/simd_affine_gap_policy.hpp>
#include <seqan3/alignment/pairwise/policy/simd_find_optimum_policy.hpp>
#include <seqan3/alignment/pairwise/alignment_algorithm.hpp>
#include <seqan3/alignment/pairwise/align_result_selector.hpp>
#include <seqan3/alignment/pairwise/alignment_result.hpp>
#include <seqan3/alignment/pairwise/detail/pairwise_alignment_algorithm.hpp>
#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/alignment/pairwise/detail/concept.hpp>
#include <seqan3/alignment/pairwise/edit_distance_algorithm.hpp>
#include <seqan3/alignment/pairwise/execution/alignment_executor_two_way.hpp>
#include <seqan3/alignment/scoring/detail/simd_match_mismatch_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/core/concept/tuple.hpp>
#include <seqan3/core/simd/simd.hpp>
#include <seqan3/core/type_traits/deferred_crtp_base.hpp>
#include <seqan3/core/type_traits/lazy.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/core/type_list/type_list.hpp>
#include <seqan3/range/views/type_reduce.hpp>
#include <seqan3/range/views/zip.hpp>

namespace seqan3::detail
{

/*!\brief A transformation trait to extract the score type used within the seqan3::align_cfg::result object.
 * \implements seqan3::transformation_trait
 * \ingroup pairwise_alignment
 *
 * \tparam config_t The alignment configuration type; must be of type seqan3::configuration and must contain
 *                  a seqan3::align_cfg::result object.
 */
template <typename config_t>
//!\cond
    requires is_type_specialisation_of_v<config_t, configuration> && config_t::template exists<align_cfg::result>()
//!\endcond
struct align_config_result_score
{
private:
    //!\brief Helper type definition to store the type of the result config.
    using result_config_t = std::remove_reference_t<decltype(seqan3::get<align_cfg::result>(std::declval<config_t>()))>;

public:
    //!\brief The score type used for the alignment result configuration.
    using type = typename result_config_t::score_type;
};

/*!\brief Provides several contracts to test when configuring the alignment algorithm.
 * \ingroup pairwise_alignment
 * \tparam range_type            The type of the range containing sequences to be aligned.
 * \tparam alignment_config_type The type of the alignment configuration.
 *
 * \details
 *
 * This stateless helper class provides several contract testing functions for the alignment configuration.
 */
template <typename range_type,
          typename alignment_config_type>
struct alignment_contract
{
private:
    /*!\brief Auxiliary member types
     * \{
     */
    //!\brief Range type with removed references.
    using unref_range_type = std::remove_reference_t<range_type>;
    //!\brief The type of the first sequence.
    using first_seq_t  = std::tuple_element_t<0, std::ranges::range_value_t<unref_range_type>>;
    //!\brief The type of the second sequence.
    using second_seq_t = std::tuple_element_t<1, std::ranges::range_value_t<unref_range_type>>;
    //!\}

public:
    //!\brief Tests whether the value type of `range_type` is a tuple with exactly 2 members.
    constexpr static bool expects_tuple_like_value_type()
    {
        return tuple_like<alignment_config_type> &&
               std::tuple_size_v<std::ranges::range_value_t<unref_range_type>> == 2;
    }

    //!\brief Tests whether the scoring scheme is set and can be invoked with the sequences passed.
    constexpr static bool expects_valid_scoring_scheme()
    {
        if constexpr (alignment_config_type::template exists<align_cfg::scoring>())
        {
            using scoring_type = std::remove_reference_t<
                                    decltype(get<align_cfg::scoring>(std::declval<alignment_config_type>()).value)
                                 >;
            return static_cast<bool>(scoring_scheme<scoring_type,
                                                   std::ranges::range_value_t<first_seq_t>,
                                                   std::ranges::range_value_t<second_seq_t>>);
        }
        else
        {
            return false;
        }
    }

    //!\brief Expects alignment configurations.
    constexpr static bool expects_alignment_configuration()
    {
        return alignment_config_type::template exists<align_cfg::mode>();
    }
};

/*!\brief Configures the alignment algorithm given the sequences and the configuration object.
 * \implements seqan3::transformation_trait
 * \ingroup pairwise_alignment
 */
struct alignment_configurator
{
private:

    /*!\brief Transformation trait that chooses the correct matrix policy.
     * \tparam traits_t The alignment configuration traits type.
     */
    template <typename traits_t>
    struct select_matrix_policy
    {
    private:
        //!\brief Indicates whether only the coordinate is required to compute the alignment.
        static constexpr bool only_coordinates = traits_t::result_type_rank < with_front_coordinate_type::rank;

        //!\brief The selected score matrix for either banded or unbanded alignments.
        using score_matrix_t = std::conditional_t<traits_t::is_banded,
                                                  alignment_score_matrix_one_column_banded<typename traits_t::score_type>,
                                                  alignment_score_matrix_one_column<typename traits_t::score_type>>;
        //!\brief The selected trace matrix for either banded or unbanded alignments.
        using trace_matrix_t = std::conditional_t<traits_t::is_banded,
                                                  alignment_trace_matrix_full_banded<typename traits_t::trace_type,
                                                                                     only_coordinates>,
                                                  alignment_trace_matrix_full<typename traits_t::trace_type,
                                                                              only_coordinates>>;

    public:
        //!\brief The matrix policy based on the configurations given by `config_type`.
        using type = deferred_crtp_base<alignment_matrix_policy, score_matrix_t, trace_matrix_t>;
    };

    /*!\brief Transformation trait that chooses the correct gap policy.
     * \tparam traits_t The alignment configuration traits type.
     */
    template <typename traits_t>
    struct select_gap_policy
    {
    private:
        //!\brief The score type for the alignment computation.
        using score_t = typename traits_t::score_type;
        //!\brief The is_local constant converted to a type.
        using is_local_t = std::bool_constant<traits_t::is_local>;

    public:
        //!\brief The matrix policy based on the configurations given by `config_type`.
        using type = std::conditional_t<traits_t::is_vectorised,
                                        deferred_crtp_base<simd_affine_gap_policy, score_t, is_local_t>,
                                        deferred_crtp_base<affine_gap_policy, score_t, is_local_t>>;
    };

    /*!\brief Transformation trait that chooses the correct find optimum policy.
     * \implements seqan3::transformation_trait
     *
     * \tparam tarits_t The alignment algorithm traits.
     * \tparam policy_traits_t The configured traits for the policy.
     */
    template <typename traits_t, typename policy_traits_t>
    struct select_find_optimum_policy
    {
    private:
        //!\brief The score type for the alignment computation.
        using score_t = typename traits_t::score_type;
        //!\brief A bool constant to disambiguate true global alignments.
        static constexpr bool is_global_alignment = traits_t::is_global && !traits_t::is_aligned_ends;

    public:
        //!\brief The find optimum policy for either scalar or vectorised alignment.
        using type = std::conditional_t<traits_t::is_vectorised,
                                        deferred_crtp_base<simd_find_optimum_policy,
                                                           score_t,
                                                           std::bool_constant<is_global_alignment>,
                                                           policy_traits_t>,
                                        deferred_crtp_base<find_optimum_policy, policy_traits_t>>;
    };

public:
    /*!\brief Configures the algorithm.
     * \tparam sequences_t The range type containing the sequence pairs; must model std::ranges::forward_range.
     * \tparam config_t    The alignment configuration type; must be a specialisation of seqan3::configuration.
     * \param[in] cfg      The configuration object.
     *
     * \returns a std::pair over std::function wrapper of the configured alignment algorithm and the adapted
     *          alignment configuration.
     *
     * \details
     *
     * This function reads the seqan3::configuration object and generates the corresponding alignment algorithm type.
     * During this process some runtime configurations are converted to static configurations if required.
     * In case of a missing configuration that has a default, e.g. the seqan3::align_cfg::result option, the
     * default version of this configuration element is added to the passed configuration object.
     * The return type is a std::pair over a std::function object and the adapted configuration object. Thus, the
     * calling function has access to the possibly modified configuration object.
     * The function object type is determined using the following type trait:
     *
     * \include snippet/alignment/pairwise/alignment_configurator.cpp
     *
     *
     * The arguments to the function object are two ranges, which always need to be passed as lvalue references.
     * Note that even if they are not passed as const lvalue reference (which is not possible, since not all views are
     * const-iterable), they are not modified within the alignment algorithm.
     */
    template <align_pairwise_range_input sequences_t, typename config_t>
    //!\cond
        requires is_type_specialisation_of_v<config_t, configuration>
    //!\endcond
    static constexpr auto configure(config_t const & cfg)
    {

        if constexpr (!config_t::template exists<align_cfg::result>())
        {
            // Set the default result value to be computed.
            return configure<sequences_t>(cfg | align_cfg::result{with_score});
        }
        else
        {
            // ----------------------------------------------------------------------------
            // Configure the type-erased alignment function.
            // ----------------------------------------------------------------------------

            using first_seq_t = std::tuple_element_t<0, std::ranges::range_value_t<sequences_t>>;
            using second_seq_t = std::tuple_element_t<1, std::ranges::range_value_t<sequences_t>>;

            using wrapped_first_t  = type_reduce_view<first_seq_t &>;
            using wrapped_second_t = type_reduce_view<second_seq_t &>;

            // The alignment executor passes a chunk over an indexed sequence pair range to the alignment algorithm.
            using indexed_sequence_pair_range_t = typename chunked_indexed_sequence_pairs<sequences_t>::type;
            using indexed_sequence_pair_chunk_t = std::ranges::range_value_t<indexed_sequence_pair_range_t>;

            // Select the result type based on the sequences and the configuration.
            using alignment_result_value_t = typename align_result_selector<std::remove_reference_t<wrapped_first_t>,
                                                                            std::remove_reference_t<wrapped_second_t>,
                                                                            config_t>::type;
            using alignment_result_t = alignment_result<alignment_result_value_t>;
            using callback_on_result_t = std::function<void(alignment_result_t)>;
            // Define the function wrapper type.
            using function_wrapper_t = std::function<void(indexed_sequence_pair_chunk_t, callback_on_result_t)>;

            // Capture the alignment result type.
            auto config_with_result_type = cfg | align_cfg::alignment_result_capture<alignment_result_t>;

            // ----------------------------------------------------------------------------
            // Test some basic preconditions
            // ----------------------------------------------------------------------------

            using alignment_contract_t = alignment_contract<sequences_t, config_t>;

            static_assert(alignment_contract_t::expects_alignment_configuration(),
                          "Alignment configuration error: "
                          "The alignment can only be configured with alignment configurations.");

            static_assert(alignment_contract_t::expects_tuple_like_value_type(),
                          "Alignment configuration error: "
                          "The value type of the sequence ranges must model the seqan3::tuple_like "
                          "and must contain exactly 2 elements.");

            static_assert(alignment_contract_t::expects_valid_scoring_scheme(),
                          "Alignment configuration error: "
                          "Either the scoring scheme was not configured or the given scoring scheme cannot be invoked with "
                          "the value types of the passed sequences.");

            // ----------------------------------------------------------------------------
            // Configure the algorithm
            // ----------------------------------------------------------------------------

            // Use default edit distance if gaps are not set.
            auto const & gaps = config_with_result_type.template value_or<align_cfg::gap>(gap_scheme{gap_score{-1}});
            auto const & scoring_scheme = get<align_cfg::scoring>(cfg).value;
            auto align_ends_cfg = config_with_result_type.template value_or<align_cfg::aligned_ends>(free_ends_none);

            if constexpr (config_t::template exists<align_cfg::mode<detail::global_alignment_type>>())
            {
                // Only use edit distance if ...
                if (gaps.get_gap_open_score() == 0 &&  // gap open score is not set,
                    !(align_ends_cfg[2] || align_ends_cfg[3]) && // none of the free end gaps are set for second seq,
                    align_ends_cfg[0] == align_ends_cfg[1]) // free ends for leading and trailing gaps are equal in first seq.
                {
                    // TODO: Instead of relying on nucleotide scoring schemes we need to be able to determine the edit distance
                    //       option via the scheme.
                    if constexpr (is_type_specialisation_of_v<remove_cvref_t<decltype(scoring_scheme)>,
                                                              nucleotide_scoring_scheme>)
                    {
                        if ((scoring_scheme.score('A'_dna15, 'A'_dna15) == 0) &&
                            (scoring_scheme.score('A'_dna15, 'C'_dna15)) == -1)
                        {
                            return std::pair{configure_edit_distance<function_wrapper_t>(config_with_result_type),
                                             config_with_result_type};
                        }
                    }
                }
            }

            // ----------------------------------------------------------------------------
            // Check if invalid configuration was used.
            // ----------------------------------------------------------------------------

            // Do not allow max error configuration for alignments not computing the edit distance.
            if (config_t::template exists<align_cfg::max_error>())
                throw invalid_alignment_configuration{"The align_cfg::max_error configuration is only allowed for "
                                                      "the specific edit distance computation."};
            // Configure the alignment algorithm.
            return std::pair{configure_scoring_scheme<function_wrapper_t>(config_with_result_type),
                             config_with_result_type};
        }
    }

private:
    /*!\brief Configures the edit distance algorithm.
     * \tparam function_wrapper_t The invocable alignment function type-erased via std::function.
     * \tparam config_t           The alignment configuration type.
     * \param[in] cfg             The passed configuration object.
     */
    template <typename function_wrapper_t, typename config_t>
    static constexpr function_wrapper_t configure_edit_distance(config_t const & cfg)
    {
        using traits_t = alignment_configuration_traits<config_t>;

        // ----------------------------------------------------------------------------
        // Unsupported configurations
        // ----------------------------------------------------------------------------

        if constexpr (traits_t::is_banded)
            throw invalid_alignment_configuration{"Banded alignments are yet not supported."};

        // ----------------------------------------------------------------------------
        // Configure semi-global alignment
        // ----------------------------------------------------------------------------

        // Get the value for the sequence ends configuration.
        auto align_ends_cfg = cfg.template value_or<align_cfg::aligned_ends>(free_ends_none);
        using align_ends_cfg_t = remove_cvref_t<decltype(align_ends_cfg)>;

        auto configure_edit_traits = [&] (auto is_semi_global)
        {
            struct edit_traits_type
            {
                using is_semi_global_type [[maybe_unused]] = remove_cvref_t<decltype(is_semi_global)>;
            };

            edit_distance_algorithm<remove_cvref_t<config_t>, edit_traits_type> algorithm{cfg};
            return function_wrapper_t{std::move(algorithm)};
        };

        // Check if it has free ends set for the first sequence trailing gaps.
        auto has_free_ends_trailing = [&] (auto first) constexpr
        {
            if constexpr (!decltype(first)::value)
            {
                return configure_edit_traits(std::false_type{});
            }
            else if constexpr (align_ends_cfg_t::template is_static<1>())
            {
#if defined(__GNUC__) && __GNUC__ >= 9
                constexpr bool free_ends_trailing = align_ends_cfg_t::template get_static<1>();
                return configure_edit_traits(std::integral_constant<bool, free_ends_trailing>{});
#else // ^^^ workaround / no workaround vvv
                return configure_edit_traits(std::integral_constant<bool, align_ends_cfg_t::template get_static<1>()>{});
#endif // defined(__GNUC__) && __GNUC__ >= 9
            }
            else // Resolve correct property at runtime.
            {
                if (align_ends_cfg[1])
                    return configure_edit_traits(std::true_type{});
                else
                    return configure_edit_traits(std::false_type{});
            }
        };

        // Check if it has free ends set for the first sequence leading gaps.
        // If possible use static information.
        if constexpr (align_ends_cfg_t::template is_static<0>())
            return has_free_ends_trailing(std::integral_constant<bool, align_ends_cfg_t::template get_static<0>()>{});
        else // Resolve correct property at runtime.
        {
            if (align_ends_cfg[0])
                return has_free_ends_trailing(std::true_type{});
            else
                return has_free_ends_trailing(std::false_type{});
        }
    }

    /*!\brief Configures the dynamic programming matrix initialisation accoring to seqan3::align_cfg::aligned_ends
     *        settings.
     *
     * \tparam function_wrapper_t The invocable alignment function type-erased via std::function.
     * \tparam policies_t A template parameter pack for the already configured policy types.
     * \tparam config_t The alignment configuration type.
     *
     * \param[in] cfg The passed configuration object.
     *
     * \returns the configured alignment algorithm.
     *
     * \details
     *
     * The matrix initialisation depends on the settings for the leading gaps for the first and the second sequence
     * within the seqan3::align_cfg::aligned_ends configuration element.
     */
    template <typename function_wrapper_t, typename ...policies_t, typename config_t>
    static constexpr function_wrapper_t configure_free_ends_initialisation(config_t const & cfg);

    /*!\brief Configures the search space for the alignment algorithm according to seqan3::align_cfg::aligned_ends
     *        settings.
     *
     * \tparam function_wrapper_t The invocable alignment function type-erased via std::function.
     * \tparam policies_t A template parameter pack for the already configured policy types.
     * \tparam config_t The alignment configuration type.
     *
     * \param[in] cfg The passed configuration object.
     *
     * \returns the configured alignment algorithm.
     *
     * \details
     *
     * This option is configured in the seqan3::align_cfg::aligned_ends configuration element according to
     * the settings for the trailing gaps of the first and the second sequence.
     */
    template <typename function_wrapper_t, typename ...policies_t, typename config_t>
    static constexpr function_wrapper_t configure_free_ends_optimum_search(config_t const & cfg);

    /*!\brief Configures the scoring scheme to use for the alignment computation.
     *
     * \tparam function_wrapper_t The invocable alignment function type-erased via std::function.
     * \tparam config_t The alignment configuration type.
     *
     * \param[in] cfg The passed configuration object.
     *
     * \returns the configured alignment algorithm.
     *
     * \details
     *
     * The correct scoring scheme is selected based on the vectorisation mode. If no vectorisation is enabled, the
     * scoring scheme is the one configured in seqan3::align_config::scoring. If vectorisation is enabled, then the
     * appropriate scoring scheme for the vectorised alignment algorithm is selected. This involves checking whether the
     * passed scoring scheme is a matrix or a simple scoring scheme, which has only mismatch and match costs.
     */
    template <typename function_wrapper_t, typename config_t>
    static constexpr function_wrapper_t configure_scoring_scheme(config_t const & cfg);

    /*!\brief Constructs the actual alignment algorithm wrapped in the passed std::function object.
     *
     * \tparam function_wrapper_t The invocable alignment function type-erased via std::function.
     * \tparam policies_t A template parameter pack for the already configured policy types.
     * \tparam config_t The alignment configuration type.
     *
     * \param[in] cfg The passed configuration object.
     *
     * \returns the configured alignment algorithm.
     *
     * \details
     *
     * Configures the matrix and the gap policy and constructs the algorithm with the configured policies.
     */
    template <typename function_wrapper_t, typename ...policies_t, typename config_t>
    static constexpr function_wrapper_t make_algorithm(config_t const & cfg)
    {
        using traits_t = alignment_configuration_traits<config_t>;

        // Temporarily we will use the new and the old alignment implementation in order to
        // refactor step-by-step to the new implementation. The new implementation will be tested in
        // macrobenchmarks to show that it maintains a high performance.

        // Use old alignment implementation if...
        if constexpr (traits_t::is_local ||            // it is a local alignment,
                      traits_t::is_aligned_ends ||     // it has aligned ends configured,
                      traits_t::is_vectorised ||       // it is vectorised,
                      traits_t::is_banded ||           // it is banded,
                      traits_t::is_debug ||            // it runs in debug mode,
                      traits_t::result_type_rank > 0)  // it computes more than the score.
        {
            using matrix_policy_t = typename select_matrix_policy<traits_t>::type;
            using gap_policy_t = typename select_gap_policy<traits_t>::type;

            return alignment_algorithm<config_t, matrix_policy_t, gap_policy_t, policies_t...>{cfg};
        }
        else  // Use new alignment algorithm implementation.
        {
            using optimum_tracker_policy_t = policy_optimum_tracker<config_t>;
            using gap_cost_policy_t = policy_affine_gap_recursion<config_t>;

            return pairwise_alignment_algorithm<config_t, gap_cost_policy_t, optimum_tracker_policy_t>{cfg};
        }
    }
};

//!\cond
template <typename function_wrapper_t, typename config_t>
constexpr function_wrapper_t alignment_configurator::configure_scoring_scheme(config_t const & cfg)
{
    using traits_t = alignment_configuration_traits<config_t>;

    using alignment_scoring_scheme_t =
        lazy_conditional_t<traits_t::is_vectorised,
                           lazy<simd_match_mismatch_scoring_scheme,
                                typename traits_t::score_type,
                                typename traits_t::scoring_scheme_alphabet_type,
                                typename traits_t::alignment_mode_type>,
                                typename traits_t::scoring_scheme_type>;

    using scoring_scheme_policy_t = deferred_crtp_base<scoring_scheme_policy, alignment_scoring_scheme_t>;
    return configure_free_ends_initialisation<function_wrapper_t, scoring_scheme_policy_t>(cfg);
}

// This function returns a std::function object which can capture runtime dependent alignment algorithm types through
// a fixed invocation interface which is already defined by the caller of this function.
template <typename function_wrapper_t, typename ...policies_t, typename config_t>
constexpr function_wrapper_t alignment_configurator::configure_free_ends_initialisation(config_t const & cfg)
{
    using traits_t = alignment_configuration_traits<config_t>;
    // Get the value for the sequence ends configuration.
    auto align_ends_cfg = cfg.template value_or<align_cfg::aligned_ends>(free_ends_none);
    using align_ends_cfg_t = decltype(align_ends_cfg);

    // This lambda augments the initialisation policy of the alignment algorithm
    // with the aligned_ends configuration from before.
    auto configure_leading_both = [&] (auto first_seq, auto second_seq) constexpr
    {
        // Define the trait for the initialisation policy
        struct policy_trait_type
        {
            using free_first_leading_t  [[maybe_unused]] = decltype(first_seq);
            using free_second_leading_t [[maybe_unused]] = decltype(second_seq);
        };

        // Make initialisation policy a deferred CRTP base and delegate to configure the find optimum policy.
        using gap_init_policy_t = deferred_crtp_base<affine_gap_init_policy, policy_trait_type>;
        return configure_free_ends_optimum_search<function_wrapper_t, policies_t..., gap_init_policy_t>(cfg);
    };

    if constexpr (traits_t::is_local)
    {
        return configure_leading_both(std::true_type{}, std::true_type{});
    }
    else
    {
        // This lambda determines the initialisation configuration for the second sequence given
        // the leading gap property for it.
        auto configure_leading_second = [&] (auto first) constexpr
        {
            // If possible use static information.
            if constexpr (align_ends_cfg_t::template is_static<2>())
            {
                using second_t = std::integral_constant<bool, align_ends_cfg_t::template get_static<2>()>;
                return configure_leading_both(first, second_t{});
            }
            else
            {   // Resolve correct property at runtime.
                if (align_ends_cfg[2])
                    return configure_leading_both(first, std::true_type{});
                else
                    return configure_leading_both(first, std::false_type{});
            }
        };

        // Here the initialisation configuration for the first sequence is determined given
        // the leading gap property for it.
        // If possible use static information.
        if constexpr (align_ends_cfg_t::template is_static<0>())
        {
            using first_t = std::integral_constant<bool, align_ends_cfg_t::template get_static<0>()>;
            return configure_leading_second(first_t{});
        }
        else
        {  // Resolve correct property at runtime.
            if (align_ends_cfg[0])
                return configure_leading_second(std::true_type{});
            else
                return configure_leading_second(std::false_type{});
        }
    }
}
//!\endcond

//!\cond
// This function returns a std::function object which can capture runtime dependent alignment algorithm types through
// a fixed invocation interface which is already defined by the caller of this function.
template <typename function_wrapper_t, typename ...policies_t, typename config_t>
constexpr function_wrapper_t alignment_configurator::configure_free_ends_optimum_search(config_t const & cfg)
{
    using traits_t = alignment_configuration_traits<config_t>;

    // Get the value for the sequence ends configuration.
    auto align_ends_cfg = cfg.template value_or<align_cfg::aligned_ends>(free_ends_none);
    using align_ends_cfg_t = decltype(align_ends_cfg);

    // This lambda augments the find optimum policy of the alignment algorithm with the
    // respective aligned_ends configuration.
    auto configure_trailing_both = [&] (auto first_seq, auto second_seq) constexpr
    {
        struct policy_trait_type
        {
            using find_in_every_cell_type  [[maybe_unused]] = std::bool_constant<traits_t::is_local>;
            using find_in_last_row_type    [[maybe_unused]] = decltype(first_seq);
            using find_in_last_column_type [[maybe_unused]] = decltype(second_seq);
        };

        // We need to select the correct policy based on the configuration traits.
        using find_optimum_t = typename select_find_optimum_policy<traits_t, policy_trait_type>::type;
        return make_algorithm<function_wrapper_t, policies_t..., find_optimum_t>(cfg);
    };

    if constexpr (traits_t::is_local)
    {
        return configure_trailing_both(std::true_type{}, std::true_type{});
    }
    else
    {
        // This lambda determines the lookup configuration for the second sequence given
        // the trailing gap property for it.
        auto configure_trailing_second = [&] (auto first) constexpr
        {
            // If possible use static information.
            if constexpr (align_ends_cfg_t::template is_static<3>())
            {
                using second_t = std::integral_constant<bool, align_ends_cfg_t::template get_static<3>()>;
                return configure_trailing_both(first, second_t{});
            }
            else
            { // Resolve correct property at runtime.
                if (align_ends_cfg[3])
                    return configure_trailing_both(first, std::true_type{});
                else
                    return configure_trailing_both(first, std::false_type{});
            }
        };

        // Here the lookup configuration for the first sequence is determined given
        // the trailing gap property for it.
        // If possible use static information.
        if constexpr (align_ends_cfg_t::template is_static<1>())
        {
            using first_t = std::integral_constant<bool, align_ends_cfg_t::template get_static<1>()>;
            return configure_trailing_second(first_t{});
        }
        else
        { // Resolve correct property at runtime.
            if (align_ends_cfg[1])
                return configure_trailing_second(std::true_type{});
            else
                return configure_trailing_second(std::false_type{});
        }
    }
}
//!\endcond

} // namespace seqan3::detail
