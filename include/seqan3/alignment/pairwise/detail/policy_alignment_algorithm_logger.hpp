// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::policy_alignment_algorithm_logger.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/concepts>

#include <seqan3/alignment/matrix/matrix_concept.hpp>
#include <seqan3/alignment/matrix/detail/two_dimensional_matrix.hpp>
#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/core/detail/template_inspection.hpp>

namespace seqan3::detail
{

/*!\brief Implements a logger to debug the alignment algorithms.
 * \ingroup pairwise_alignment
 *
 * \tparam debug_score_matrix_t The type of debug score matrix; must model seqan3::detail::matrix.
 * \tparam debug_trace_matrix_t The type of debug trace matrix; must model seqan3::detail::matrix or must be the same as
 *                              seqan3::detail::empty_type.
 *
 * \details
 *
 * The logger offers interfaces to store each column in a debug alignment (score and trace) matrix.
 * It needs to be initialised before computing the alignment. During the alignment computation each
 * alignment column must be logged. During the result construction a copy of the debug matrices will be stored inside
 * of the generated alignment result object. This logging mechanism is only available if the
 * seqan3::align_cfg::detail::debug configuration was enabled.
 *
 * If only the score matrix is available, then the type of the debug_trace_matrix_t can be set to
 * seqan3::detail::empty_type to disable the logging of the trace matrix entirely.
 */
template <matrix debug_score_matrix_t, typename debug_trace_matrix_t>
//!\cond
    requires (matrix<debug_trace_matrix_t> || std::same_as<debug_trace_matrix_t, empty_type>)
//!\endcond
class policy_alignment_algorithm_logger
{
private:
    //!\brief Whether the trace information is required.
    static constexpr bool with_trace = !std::same_as<debug_trace_matrix_t, empty_type>;

protected:
    //!\brief The debug score matrix.
    debug_score_matrix_t debug_score_matrix{};
    //!\brief The debug trace matrix.
    debug_trace_matrix_t debug_trace_matrix{};

    template <typename>
    friend class alignment_log_transfer;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    policy_alignment_algorithm_logger() = default; //!< Default.
    policy_alignment_algorithm_logger(policy_alignment_algorithm_logger const &) = default; //!< Default.
    policy_alignment_algorithm_logger(policy_alignment_algorithm_logger &&) = default; //!< Default.
    policy_alignment_algorithm_logger & operator=(policy_alignment_algorithm_logger const &) = default; //!< Default.
    policy_alignment_algorithm_logger & operator=(policy_alignment_algorithm_logger &&) = default; //!< Default.
    ~policy_alignment_algorithm_logger() = default; //!< Default.

    /*!\brief Construction and initialisation using the alignment configuration.
     * \param[in] config The alignment configuration [not used in this context].
     */
    template <typename alignment_configuration_t>
    //!\cond
        requires is_type_specialisation_of_v<alignment_configuration_t, configuration>
    //!\endcond
    policy_alignment_algorithm_logger(alignment_configuration_t const & SEQAN3_DOXYGEN_ONLY(config))
    {}
    //!\}

    /*!\brief Initialises the local debug matrices.
     *
     * \param[in] sequence1_size The size of the first sequence.
     * \param[in] sequence2_size The size of the second sequence.
     *
     * \details
     *
     * Resizes the debug score, and if requested, the trace matrix, to the given matrix dimensions.
     *
     * ### Exception
     *
     * Might throw std::bad_alloc if the requested matrix size exceeds the available memory.
     */
    void initialise_debug_matrices(size_t const sequence1_size, size_t const sequence2_size)
    {
        assert(sequence1_size < static_cast<uint64_t>(std::numeric_limits<int64_t>::max()));
        assert(sequence2_size < static_cast<uint64_t>(std::numeric_limits<int64_t>::max()));

        size_t const column_count = sequence1_size + 1;
        size_t const row_count = sequence2_size + 1;

        debug_score_matrix.resize(number_rows{row_count}, number_cols{column_count});

        if constexpr (with_trace)
            debug_trace_matrix.resize(number_rows{row_count}, number_cols{column_count});
    }

    /*!\brief Logs the current alignment column.
     * \tparam coordinate_column_t The type of the coordinate matrix column; must model std::ranges::input_range and
     *                             seqan3::detail::matrix_offset is std::constructible_from the range reference type.
     * \tparam alignment_column_t The type of the alignment matrix column; must model std::ranges::input_range.
     *
     * param[in] coordinate_column The current column over the coordinate matrix.
     * param[in] alignment_column The current column over the alignment matrix.
     *
     * \details
     *
     * Logs the current alignment column in the locally stored debug matrices for the score and, if applicable,
     * the trace matrix.
     * The coordinate column is used to determine the correct column offset at which the information is stored.
     * In the banded matrix the offset can be different depending on the valid area of the band within the global
     * alignment matrix.
     */
    template <std::ranges::input_range coordinate_column_t, std::ranges::input_range alignment_column_t>
    //!\cond
        requires std::constructible_from<matrix_offset, std::ranges::range_reference_t<coordinate_column_t>>
    //!\endcond
    void log_alignment_matrix_column(coordinate_column_t && coordinate_column,
                                     alignment_column_t && alignment_column)
    {
        matrix_offset column_coordinate_begin{*coordinate_column.begin()};

        assert(static_cast<size_t>(column_coordinate_begin.col) < debug_score_matrix.cols());
        assert(static_cast<size_t>(column_coordinate_begin.row) < debug_score_matrix.rows());

        std::ranges::copy(alignment_column | std::views::transform([] (auto && cell) { return cell.best_score(); }),
                          debug_score_matrix.begin() + column_coordinate_begin);

        if constexpr (with_trace)
        {
            assert(static_cast<size_t>(column_coordinate_begin.col) < debug_trace_matrix.cols());
            assert(static_cast<size_t>(column_coordinate_begin.row) < debug_trace_matrix.rows());

            std::ranges::copy(alignment_column | std::views::transform([] (auto && cell) { return cell.best_trace(); }),
                              debug_trace_matrix.begin() + column_coordinate_begin);
        }
    }
};

/*!\brief Helper function object to transfer the logging information.
 * \ingroup pairwise_alignment
 *
 * \tparam result_data_t The seqan3::detail::alignment_result_value_type.
 */
template <typename result_data_t>
struct alignment_log_transfer
{
    result_data_t & result_data; //!< The result data to receive the logged alignment matrix.

    /*!\name Constructors, destructor and assignment
     * \{
     */
    alignment_log_transfer() = delete; //!< Delete.
    alignment_log_transfer(alignment_log_transfer const &) = default; //!< Default.
    alignment_log_transfer(alignment_log_transfer &&) = default; //!< Default.
    alignment_log_transfer & operator=(alignment_log_transfer const &) = default; //!< Default.
    alignment_log_transfer & operator=(alignment_log_transfer &&) = default; //!< Default.
    ~alignment_log_transfer() = default; //!< Default.

    //!\brief Initialises the result data reference.
    alignment_log_transfer(result_data_t & result) : result_data{result}
    {}
    //!\}

    //!\brief Transfers the logged matrix information to the linked alignment result data.
    template <typename score_matrix_t, typename trace_matrix_t>
    void operator()(policy_alignment_algorithm_logger<score_matrix_t, trace_matrix_t> & logger)
    {
        using std::swap;
        swap(result_data.score_debug_matrix, logger.debug_score_matrix);

        if constexpr (!std::same_as<trace_matrix_t, empty_type>)
            swap(result_data.trace_debug_matrix, logger.debug_trace_matrix);
    }
};

} // namespace seqan3::detail
