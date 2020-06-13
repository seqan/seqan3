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

#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>

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
     * \tparam callback_t The type of the callback to invoke.
     *
     * \param[in] sequence_pair The indexed sequence pair.
     * \param[in] id The associated id.
     * \param[in] score The best alignment score.
     * \param[in] callback The callback to invoke with the generated result.
     *
     * \details
     *
     * Generates a seqan3::alignment_result object with the results computed during the alignment. Depending on the
     * seqan3::align_cfg::result configuration only the requested values are stored. In some cases some additional
     * work is done to generate the requested result. For example computing the associated alignment from the traceback
     * matrix.
     */
    template <typename sequence_pair_t,
              typename index_t,
              typename score_t,
              typename coordinate_t,
              typename alignment_matrix_t,
              typename callback_t>
    //!\cond
        requires std::invocable<callback_t, result_type>
    //!\endcond
    void make_result_and_invoke([[maybe_unused]] sequence_pair_t && sequence_pair,
                                [[maybe_unused]] index_t && id,
                                [[maybe_unused]] score_t score,
                                [[maybe_unused]] coordinate_t end_coordinate,
                                [[maybe_unused]] alignment_matrix_t alignment_matrix,
                                callback_t && callback)
    {
        using std::get;
        using invalid_t = std::nullopt_t *;

        result_type result{};

        static_assert(!std::same_as<decltype(result.data.id), invalid_t>,
                      "Invalid configuration. Expected result with id!");
        result.data.id = std::move(id);

        if constexpr (traits_type::compute_score)
        {
            static_assert(!std::same_as<decltype(result.data.score), invalid_t>,
                          "Invalid configuration. Expected result with score!");
            result.data.score = std::move(score);
        }

        if constexpr (traits_type::compute_back_coordinate)
        {
            result.data.back_coordinate.first = end_coordinate.col;
            result.data.back_coordinate.second = end_coordinate.row;
        }

        if constexpr (traits_type::compute_front_coordinate)
        {
            using std::get;
            // Get a aligned sequence builder for banded or un-banded case.
            aligned_sequence_builder builder{get<0>(sequence_pair), get<1>(sequence_pair)};

            auto trace_res = builder(alignment_matrix.trace_path(end_coordinate));
            result.data.front_coordinate.first = trace_res.first_sequence_slice_positions.first;
            result.data.front_coordinate.second = trace_res.second_sequence_slice_positions.first;

            if constexpr (traits_type::compute_sequence_alignment)
                result.data.alignment = std::move(trace_res.alignment);
        }

        callback(std::move(result));
    }
};
} // namespace seqan3::detail
