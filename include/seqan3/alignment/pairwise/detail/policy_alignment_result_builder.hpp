// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::policy_alignment_result_builder.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/matrix/detail/aligned_sequence_builder.hpp>
#include <seqan3/alignment/pairwise/detail/policy_alignment_algorithm_logger.hpp>
#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/core/detail/template_inspection.hpp>

namespace seqan3::detail
{

/*!\brief Implements the alignment result builder.
 * \ingroup pairwise_alignment
 *
 * \tparam alignment_configuration_t The type of the alignment configuration; must be a type specialisation of
 *                                   seqan3::configuration.
 *
 * \details
 *
 * Implements the interfaces to build the alignment result based on the previously selected output configurations.
 */
template <typename alignment_configuration_t>
#if !SEQAN3_WORKAROUND_GCC_93467
//!\cond
    requires is_type_specialisation_of_v<alignment_configuration_t, configuration>
//!\endcond
#endif // !SEQAN3_WORKAROUND_GCC_93467
class policy_alignment_result_builder
{
protected:
    //!\brief The configuration traits type.
    using traits_type = alignment_configuration_traits<alignment_configuration_t>;
    //!\brief The alignment result type.
    using result_type = typename traits_type::alignment_result_type;

    static_assert(!std::same_as<result_type, empty_type>, "The alignment result type was not configured.");


    /*!\name Constructors, destructor and assignment
     * \{
     */
    policy_alignment_result_builder() = default; //!< Defaulted.
    policy_alignment_result_builder(policy_alignment_result_builder const &) = default; //!< Defaulted.
    policy_alignment_result_builder(policy_alignment_result_builder &&) = default; //!< Defaulted.
    policy_alignment_result_builder & operator=(policy_alignment_result_builder const &) = default; //!< Defaulted.
    policy_alignment_result_builder & operator=(policy_alignment_result_builder &&) = default; //!< Defaulted.
    ~policy_alignment_result_builder() = default; //!< Defaulted.

    /*!\brief Construction and initialisation using the alignment configuration.
     * \param[in] config The alignment configuration [not used in this context].
     */
    policy_alignment_result_builder(alignment_configuration_t const & SEQAN3_DOXYGEN_ONLY(config))
    {}
    //!\}

    /*!\brief Builds the seqan3::alignment_result based on the given alignment result type and then invokes the
     *        given callable with the result.
     *
     * \tparam sequence_pair_t The type of the sequence pair.
     * \tparam id_t The type of the id.
     * \tparam score_t The type of the score.
     * \tparam matrix_coordinate_t The type of the matrix coordinate.
     * \tparam alignment_matrix_t The type of the alignment matrix.
     * \tparam callback_t The type of the callback to invoke.
     * \tparam transfer_log_fn_t The type of the function which transfers the logged alignment matrix.
     *
     * \param[in] sequence_pair The indexed sequence pair.
     * \param[in] id The associated id.
     * \param[in] score The best alignment score.
     * \param[in] end_positions The matrix coordinate of the best alignment score.
     * \param[in] alignment_matrix The alignment matrix to obtain the trace back from.
     * \param[in] callback The callback to invoke with the generated result.
     * \param[in] transfer_log The invocable to transfer the logged alignment matrix.
     *
     * \details
     *
     * Generates a seqan3::alignment_result object with the results computed during the alignment. Depending on the
     * \ref seqan3_align_cfg_output_configurations "seqan3::align_cfg::output_*" configuration only the requested values
     * are stored. In some cases some additional work is done to generate the requested result. For example computing
     * the associated alignment from the traceback matrix.
     *
     * The first parameter is the logger and is only available if the alignment was run in debug mode. It stores the
     * debug score and, if applicable, the debug trace matrix, which are then stored in the created alignment result.
     */
    template <typename sequence_pair_t,
              typename index_t,
              typename score_t,
              typename matrix_coordinate_t,
              typename alignment_matrix_t,
              typename callback_t,
              typename transfer_log_fn_t>
    //!\cond
        requires std::invocable<callback_t, result_type>
    //!\endcond
    void make_result_and_invoke([[maybe_unused]] sequence_pair_t && sequence_pair,
                                [[maybe_unused]] index_t && id,
                                [[maybe_unused]] score_t score,
                                [[maybe_unused]] matrix_coordinate_t end_positions,
                                [[maybe_unused]] alignment_matrix_t const & alignment_matrix,
                                callback_t && callback,
                                [[maybe_unused]] transfer_log_fn_t && transfer_log)
    {
        using std::get;
        using invalid_t = std::nullopt_t *;

        result_type result{};

        if constexpr (traits_type::output_sequence1_id)
            result.data.sequence1_id = id;

        if constexpr (traits_type::output_sequence2_id)
            result.data.sequence2_id = id;

        if constexpr (traits_type::compute_score)
        {
            static_assert(!std::same_as<decltype(result.data.score), invalid_t>,
                          "Invalid configuration. Expected result with score!");
            result.data.score = std::move(score);
        }

        if constexpr (traits_type::compute_end_positions)
        {
            static_assert(!std::same_as<decltype(result.data.end_positions), invalid_t>,
                          "Invalid configuration. Expected result with end positions!");

            result.data.end_positions.first = end_positions.col;
            result.data.end_positions.second = end_positions.row;
        }

        if constexpr (traits_type::requires_trace_information)
        {
            aligned_sequence_builder builder{get<0>(sequence_pair), get<1>(sequence_pair)};
            auto aligned_sequence_result = builder(alignment_matrix.trace_path(end_positions));

            if constexpr (traits_type::compute_begin_positions)
            {
                result.data.begin_positions.first = aligned_sequence_result.first_sequence_slice_positions.first;
                result.data.begin_positions.second = aligned_sequence_result.second_sequence_slice_positions.first;
            }

            if constexpr (traits_type::compute_sequence_alignment)
                result.data.alignment = std::move(aligned_sequence_result.alignment);
        }

        // In case we run in debug mode, we need to store the score and possibly trace matrix.
        if constexpr (traits_type::is_debug)
            transfer_log(result.data);

        callback(std::move(result));
    }

};
} // namespace seqan3::detail
