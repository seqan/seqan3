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

#include <seqan3/alignment/configuration/align_config_output.hpp>
#include <seqan3/alignment/configuration/align_config_result_type.hpp>
#include <seqan3/alignment/matrix/detail/alignment_score_matrix_one_column.hpp>
#include <seqan3/alignment/matrix/detail/alignment_score_matrix_one_column_banded.hpp>
#include <seqan3/alignment/matrix/detail/alignment_trace_matrix_full.hpp>
#include <seqan3/alignment/matrix/detail/alignment_trace_matrix_full_banded.hpp>
#include <seqan3/alignment/matrix/detail/combined_score_and_trace_matrix.hpp>
#include <seqan3/alignment/matrix/detail/score_matrix_single_column.hpp>
#include <seqan3/alignment/matrix/detail/trace_matrix_full.hpp>
#include <seqan3/alignment/pairwise/detail/concept.hpp>
#include <seqan3/alignment/pairwise/detail/pairwise_alignment_algorithm.hpp>
#include <seqan3/alignment/pairwise/detail/pairwise_alignment_algorithm_banded.hpp>
#include <seqan3/alignment/pairwise/detail/policy_alignment_matrix.hpp>
#include <seqan3/alignment/pairwise/detail/policy_alignment_result_builder.hpp>
#include <seqan3/alignment/pairwise/detail/policy_affine_gap_recursion.hpp>
#include <seqan3/alignment/pairwise/detail/policy_affine_gap_recursion_banded.hpp>
#include <seqan3/alignment/pairwise/detail/policy_affine_gap_with_trace_recursion.hpp>
#include <seqan3/alignment/pairwise/detail/policy_affine_gap_with_trace_recursion_banded.hpp>
#include <seqan3/alignment/pairwise/detail/policy_optimum_tracker_simd.hpp>
#include <seqan3/alignment/pairwise/detail/policy_optimum_tracker.hpp>
#include <seqan3/alignment/pairwise/detail/policy_scoring_scheme.hpp>
#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/alignment/pairwise/edit_distance_algorithm.hpp>
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
#include <seqan3/alignment/scoring/detail/simd_match_mismatch_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/detail/simd_matrix_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/core/detail/deferred_crtp_base.hpp>
#include <seqan3/core/detail/template_inspection.hpp>
#include <seqan3/utility/simd/simd.hpp>
#include <seqan3/utility/tuple/concept.hpp>
#include <seqan3/utility/type_traits/lazy_conditional.hpp>
#include <seqan3/utility/views/type_reduce.hpp>
#include <seqan3/utility/views/zip.hpp>

namespace seqan3::detail
{

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
        if constexpr (alignment_config_type::template exists<align_cfg::scoring_scheme>())
        {
            using scoring_type =
                std::remove_reference_t<
                    decltype(get<align_cfg::scoring_scheme>(std::declval<alignment_config_type>()).scheme)>;
            return static_cast<bool>(scoring_scheme_for<scoring_type,
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
        const bool is_global = alignment_config_type::template exists<seqan3::align_cfg::method_global>();
        const bool is_local = alignment_config_type::template exists<seqan3::align_cfg::method_local>();

        return (is_global || is_local);
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
        static constexpr bool only_coordinates = !(traits_t::compute_begin_positions ||
                                                   traits_t::compute_sequence_alignment);

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
    template <typename traits_t>
    struct select_find_optimum_policy
    {
    private:
        //!\brief The score type for the alignment computation.
        using score_t = typename traits_t::score_type;

    public:
        //!\brief The find optimum policy for either scalar or vectorised alignment.
        using type = std::conditional_t<traits_t::is_vectorised,
                                        deferred_crtp_base<simd_find_optimum_policy, score_t>,
                                        deferred_crtp_base<find_optimum_policy>>;
    };

    //!\brief Selects either the banded or the unbanded alignment algorithm based on the given traits type.
    template <typename traits_t, typename ...args_t>
    using select_alignment_algorithm_t = lazy_conditional_t<traits_t::is_banded,
                                                            lazy<pairwise_alignment_algorithm_banded, args_t...>,
                                                            lazy<pairwise_alignment_algorithm, args_t...>>;

    /*!\brief Selects the gap recursion policy.
     * \tparam config_t The alignment configuration type.
     */
    template <typename config_t>
    struct select_gap_recursion_policy
    {
    private:
        //!\brief The traits type.
        using traits_type = alignment_configuration_traits<config_t>;
        //!\brief A flag indicating if trace is required.
        static constexpr bool with_trace = traits_type::requires_trace_information;

        //!\brief The gap recursion policy.
        using gap_recursion_policy_type = std::conditional_t<with_trace,
                                                             policy_affine_gap_with_trace_recursion<config_t>,
                                                             policy_affine_gap_recursion<config_t>>;
        //!\brief The banded gap recursion policy.
        using banded_gap_recursion_policy_type =
                    std::conditional_t<with_trace,
                                       policy_affine_gap_with_trace_recursion_banded<config_t>,
                                       policy_affine_gap_recursion_banded<config_t>>;
    public:
        //!\brief The configured recursion policy.
        using type = std::conditional_t<traits_type::is_banded,
                                        banded_gap_recursion_policy_type,
                                        gap_recursion_policy_type>;
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
     * In case of a missing configuration that has a default, e.g. the \ref seqan3_align_cfg_output_configurations
     * "seqan3::align_cfg::output_*" options, the default version of this configuration element is added to the passed
     * configuration object.
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
        auto config_with_output = maybe_default_output(cfg);
        using config_with_output_t = decltype(config_with_output);

        // ----------------------------------------------------------------------------
        // Configure the type-erased alignment function.
        // ----------------------------------------------------------------------------

        using first_seq_t = std::tuple_element_t<0, std::ranges::range_value_t<sequences_t>>;
        using second_seq_t = std::tuple_element_t<1, std::ranges::range_value_t<sequences_t>>;

        using wrapped_first_t  = type_reduce_t<first_seq_t &>;
        using wrapped_second_t = type_reduce_t<second_seq_t &>;

        // The alignment executor passes a chunk over an indexed sequence pair range to the alignment algorithm.
        using indexed_sequence_pair_range_t = typename chunked_indexed_sequence_pairs<sequences_t>::type;
        using indexed_sequence_pair_chunk_t = std::ranges::range_value_t<indexed_sequence_pair_range_t>;

        // Select the result type based on the sequences and the configuration.
        using alignment_result_value_t = typename align_result_selector<std::remove_reference_t<wrapped_first_t>,
                                                                        std::remove_reference_t<wrapped_second_t>,
                                                                        config_with_output_t>::type;
        using alignment_result_t = alignment_result<alignment_result_value_t>;
        using callback_on_result_t = std::function<void(alignment_result_t)>;
        // Define the function wrapper type.
        using function_wrapper_t = std::function<void(indexed_sequence_pair_chunk_t, callback_on_result_t)>;

        // Capture the alignment result type.
        auto config_with_result_type = config_with_output | align_cfg::detail::result_type<alignment_result_t>{};

        // ----------------------------------------------------------------------------
        // Test some basic preconditions
        // ----------------------------------------------------------------------------

        using alignment_contract_t = alignment_contract<sequences_t, config_with_output_t>;

        static_assert(alignment_contract_t::expects_alignment_configuration(),
                      "Alignment configuration error: "
                      "The alignment can only be configured with alignment configurations.");

        static_assert(alignment_contract_t::expects_tuple_like_value_type(),
                      "Alignment configuration error: "
                      "The value type of the sequence ranges must model the seqan3::tuple_like and must contain "
                      "exactly 2 elements.");

        static_assert(alignment_contract_t::expects_valid_scoring_scheme(),
                      "Alignment configuration error: "
                      "Either the scoring scheme was not configured or the given scoring scheme cannot be invoked with "
                      "the value types of the passed sequences.");

        // ----------------------------------------------------------------------------
        // Configure the algorithm
        // ----------------------------------------------------------------------------

        // Use default edit distance if gaps are not set.
        align_cfg::gap_cost_affine edit_gap_cost{};
        auto const & gap_cost = config_with_result_type.get_or(edit_gap_cost);
        auto const & scoring_scheme = get<align_cfg::scoring_scheme>(cfg).scheme;

        if constexpr (config_t::template exists<seqan3::align_cfg::method_global>())
        {
            // Only use edit distance if ...
            auto method_global_cfg = get<seqan3::align_cfg::method_global>(config_with_result_type);
            // Only use edit distance if ...
            if (gap_cost.open_score == 0 &&  // gap open score is not set,
                !(method_global_cfg.free_end_gaps_sequence2_leading ||
                  method_global_cfg.free_end_gaps_sequence2_trailing) && // none of the free end gaps are set for second seq,
                (method_global_cfg.free_end_gaps_sequence1_leading ==
                 method_global_cfg.free_end_gaps_sequence1_trailing)) // free ends for leading and trailing gaps are equal in first seq.
            {
                // TODO: Instead of relying on nucleotide scoring schemes we need to be able to determine the edit distance
                //       option via the scheme.
                if constexpr (is_type_specialisation_of_v<std::remove_cvref_t<decltype(scoring_scheme)>,
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

        // Do not allow min score configuration for alignments not computing the edit distance.
        if (config_t::template exists<align_cfg::min_score>())
            throw invalid_alignment_configuration{"The align_cfg::min_score configuration is only allowed for the "
                                                  "specific edit distance computation."};
        // Configure the alignment algorithm.
        return std::pair{configure_scoring_scheme<function_wrapper_t>(config_with_result_type),
                         config_with_result_type};
    }

private:
    /*!\brief Adds maybe the default output arguments if the user did not provide any.
     *
     * \tparam config_t The original type of the alignment configuration.
     *
     * \param[in] config The original user configuration to check.
     *
     * \returns Either the original config if the user specified any output configuration or a new config with all
     *          output options enabled.
     */
    template <typename config_t>
    static constexpr auto maybe_default_output(config_t const & config) noexcept
    {
        using traits_t = alignment_configuration_traits<config_t>;

        if constexpr (traits_t::has_output_configuration)
            return config;
        else
            return config | align_cfg::output_score{} |
                            align_cfg::output_begin_position{} |
                            align_cfg::output_end_position{} |
                            align_cfg::output_alignment{} |
                            align_cfg::output_sequence1_id{} |
                            align_cfg::output_sequence2_id{};
    }

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
        auto method_global_cfg = cfg.get_or(align_cfg::method_global{});

        auto configure_edit_traits = [&] (auto is_semi_global)
        {
            struct edit_traits_type
            {
                using is_semi_global_type [[maybe_unused]] = std::remove_cvref_t<decltype(is_semi_global)>;
            };

            edit_distance_algorithm<std::remove_cvref_t<config_t>, edit_traits_type> algorithm{cfg};
            return function_wrapper_t{std::move(algorithm)};
        };

        // Check if it has free ends set for the first sequence trailing gaps.
        auto has_free_ends_trailing = [&] (auto first) constexpr
        {
            if constexpr (!decltype(first)::value)
            {
                return configure_edit_traits(std::false_type{});
            }
            else // Resolve correct property at runtime.
            {
                if (method_global_cfg.free_end_gaps_sequence1_trailing)
                    return configure_edit_traits(std::true_type{});
                else
                    return configure_edit_traits(std::false_type{});
            }
        };

        // Check if it has free ends set for the first sequence leading gaps.
        if (method_global_cfg.free_end_gaps_sequence1_leading)
            return has_free_ends_trailing(std::true_type{});
        else
            return has_free_ends_trailing(std::false_type{});
    }

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
     * scoring scheme is the one configured in seqan3::align_cfg::scoring. If vectorisation is enabled, then the
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
        if constexpr (traits_t::is_local ||                                          // it is a local alignment,
                      traits_t::is_debug ||                                          // it runs in debug mode,
                      traits_t::compute_sequence_alignment ||                        // it computes more than the begin position.
                     (traits_t::is_banded && traits_t::compute_begin_positions) ||   // banded && more than end positions.
                     (traits_t::is_vectorised && traits_t::compute_end_positions))   // simd and more than the score.
        {
            using matrix_policy_t = typename select_matrix_policy<traits_t>::type;
            using gap_policy_t = typename select_gap_policy<traits_t>::type;
            using find_optimum_t = typename select_find_optimum_policy<traits_t>::type;
            using gap_init_policy_t = deferred_crtp_base<affine_gap_init_policy>;

            return alignment_algorithm<config_t, matrix_policy_t, gap_policy_t, find_optimum_t, gap_init_policy_t, policies_t...>{cfg};
        }
        else  // Use new alignment algorithm implementation.
        {
            //----------------------------------------------------------------------------------------------------------
            // Configure the optimum tracker policy.
            //----------------------------------------------------------------------------------------------------------

            using scalar_optimum_updater_t = std::conditional_t<traits_t::is_banded,
                                                                max_score_banded_updater,
                                                                max_score_updater>;

            using optimum_tracker_policy_t =
                lazy_conditional_t<traits_t::is_vectorised,
                                   lazy<policy_optimum_tracker_simd, config_t, max_score_updater_simd_global>,
                                   lazy<policy_optimum_tracker, config_t, scalar_optimum_updater_t>>;

            //----------------------------------------------------------------------------------------------------------
            // Configure the gap scheme policy.
            //----------------------------------------------------------------------------------------------------------

            using gap_cost_policy_t = typename select_gap_recursion_policy<config_t>::type;

            //----------------------------------------------------------------------------------------------------------
            // Configure the result builder policy.
            //----------------------------------------------------------------------------------------------------------

            using result_builder_policy_t = policy_alignment_result_builder<config_t>;

            //----------------------------------------------------------------------------------------------------------
            // Configure the scoring scheme policy.
            //----------------------------------------------------------------------------------------------------------

            using alignment_method_t = std::conditional_t<traits_t::is_global,
                                                          seqan3::align_cfg::method_global,
                                                          seqan3::align_cfg::method_local>;

            using score_t = typename traits_t::score_type;
            using scoring_scheme_t = typename traits_t::scoring_scheme_type;
            constexpr bool is_aminoacid_scheme = is_type_specialisation_of_v<scoring_scheme_t, aminoacid_scoring_scheme>;

            using simple_simd_scheme_t = lazy_conditional_t<traits_t::is_vectorised,
                                                            lazy<simd_match_mismatch_scoring_scheme,
                                                                 score_t,
                                                                 typename traits_t::scoring_scheme_alphabet_type,
                                                                 alignment_method_t>,
                                                            void>;
            using matrix_simd_scheme_t = lazy_conditional_t<traits_t::is_vectorised,
                                                            lazy<simd_matrix_scoring_scheme,
                                                                 score_t,
                                                                 typename traits_t::scoring_scheme_alphabet_type,
                                                                 alignment_method_t>,
                                                            void>;

            using alignment_scoring_scheme_t = std::conditional_t<traits_t::is_vectorised,
                                                                  std::conditional_t<is_aminoacid_scheme,
                                                                                     matrix_simd_scheme_t,
                                                                                     simple_simd_scheme_t>,
                                                                  scoring_scheme_t>;

            using scoring_scheme_policy_t = policy_scoring_scheme<config_t, alignment_scoring_scheme_t>;

            //----------------------------------------------------------------------------------------------------------
            // Configure the alignment matrix policy.
            //----------------------------------------------------------------------------------------------------------

            using score_matrix_t = score_matrix_single_column<score_t>;
            using trace_matrix_t = trace_matrix_full<trace_directions>;

            using alignment_matrix_t = std::conditional_t<traits_t::requires_trace_information,
                                                          combined_score_and_trace_matrix<score_matrix_t,
                                                                                          trace_matrix_t>,
                                                          score_matrix_t>;
            using alignment_matrix_policy_t = policy_alignment_matrix<traits_t, alignment_matrix_t>;

            //----------------------------------------------------------------------------------------------------------
            // Configure the final alignment algorithm.
            //----------------------------------------------------------------------------------------------------------

            using algorithm_t = select_alignment_algorithm_t<traits_t,
                                                             config_t,
                                                             gap_cost_policy_t,
                                                             optimum_tracker_policy_t,
                                                             result_builder_policy_t,
                                                             scoring_scheme_policy_t,
                                                             alignment_matrix_policy_t>;
            return algorithm_t{cfg};
        }
    }
};

//!\cond
template <typename function_wrapper_t, typename config_t>
constexpr function_wrapper_t alignment_configurator::configure_scoring_scheme(config_t const & cfg)
{
    using traits_t = alignment_configuration_traits<config_t>;

    using scoring_scheme_t = typename traits_t::scoring_scheme_type;
    constexpr bool is_aminoacid_scheme = is_type_specialisation_of_v<scoring_scheme_t, aminoacid_scoring_scheme>;
    using alignment_type_t = typename std::conditional_t<traits_t::is_global,
                                                         seqan3::align_cfg::method_global,
                                                         seqan3::align_cfg::method_local>;

    using simple_simd_scheme_t = lazy_conditional_t<traits_t::is_vectorised,
                                                    lazy<simd_match_mismatch_scoring_scheme,
                                                         typename traits_t::score_type,
                                                         typename traits_t::scoring_scheme_alphabet_type,
                                                         alignment_type_t>,
                                                    void>;
    using matrix_simd_scheme_t = lazy_conditional_t<traits_t::is_vectorised,
                                                    lazy<simd_matrix_scoring_scheme,
                                                         typename traits_t::score_type,
                                                         typename traits_t::scoring_scheme_alphabet_type,
                                                         alignment_type_t>,
                                                    void>;

    using alignment_scoring_scheme_t = std::conditional_t<traits_t::is_vectorised,
                                                          std::conditional_t<is_aminoacid_scheme, matrix_simd_scheme_t, simple_simd_scheme_t>,
                                                          scoring_scheme_t>;

    using scoring_scheme_policy_t = deferred_crtp_base<scoring_scheme_policy, alignment_scoring_scheme_t>;
    return make_algorithm<function_wrapper_t, scoring_scheme_policy_t>(cfg);
}
//!\endcond
} // namespace seqan3::detail
