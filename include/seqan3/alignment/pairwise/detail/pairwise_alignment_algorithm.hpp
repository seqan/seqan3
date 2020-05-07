// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::pairwise_alignment_algorithm.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

#include <seqan3/alignment/configuration/align_config_gap.hpp>
#include <seqan3/alignment/configuration/align_config_scoring.hpp>
#include <seqan3/alignment/matrix/detail/score_matrix_single_column.hpp>
#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/alignment/scoring/gap_scheme.hpp>
#include <seqan3/core/detail/empty_type.hpp>
#include <seqan3/core/detail/type_inspection.hpp>
#include <seqan3/range/views/drop.hpp>
#include <seqan3/range/views/zip.hpp>

namespace seqan3::detail
{

/*!\brief The alignment algorithm type to compute standard pairwise alignment using dynamic programming.
 * \implements std::invocable
 * \ingroup pairwise_alignment
 *
 * \tparam alignment_configuration_t The configuration type; must be of type seqan3::configuration.
 * \tparam policies_t Variadic template argument for the different policies of this alignment algorithm.
 *
 * \details
 *
 * ### Configuration
 *
 * The first template argument is the type of the alignment configuration. The alignment configuration was used to
 * configure the `alignment algorithm type` within the seqan3::detail::alignment_configurator.
 * The algorithm computes a column based dynamic programming matrix given two sequences.
 * After the computation a user defined callback function is invoked with the computed seqan3::alignment_result.
 */
template <typename alignment_configuration_t, typename ...policies_t>
//!\cond
    requires is_type_specialisation_of_v<alignment_configuration_t, configuration>
//!\endcond
class pairwise_alignment_algorithm : protected policies_t...
{
private:
    //!\brief The alignment configuration traits type with auxiliary information extracted from the configuration type.
    using traits_type = alignment_configuration_traits<alignment_configuration_t>;
    //!\brief The type of the scoring scheme.
    using scoring_scheme_type =  typename traits_type::scoring_scheme_type;
    //!\brief The configured alignment result type.
    using alignment_result_type = typename traits_type::alignment_result_type;

    static_assert(!std::same_as<alignment_result_type, empty_type>, "Alignment result type was not configured.");

    //!\brief The configured scoring scheme.
    scoring_scheme_type m_scoring_scheme{};

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    pairwise_alignment_algorithm() = default; //!< Defaulted.
    pairwise_alignment_algorithm(pairwise_alignment_algorithm const &) = default; //!< Defaulted.
    pairwise_alignment_algorithm(pairwise_alignment_algorithm &&) = default; //!< Defaulted.
    pairwise_alignment_algorithm & operator=(pairwise_alignment_algorithm const &) = default; //!< Defaulted.
    pairwise_alignment_algorithm & operator=(pairwise_alignment_algorithm &&) = default; //!< Defaulted.
    ~pairwise_alignment_algorithm() = default; //!< Defaulted.

    /*!\brief Constructs and initialises the algorithm using the alignment configuration.
     * \param config The configuration passed into the algorithm.
     *
     * \details
     *
     * Initialises the algorithm given the user settings from the alignment configuration object.
     */
    pairwise_alignment_algorithm(alignment_configuration_t const & config) : policies_t(config)...
    {
        m_scoring_scheme = seqan3::get<align_cfg::scoring>(config).value;
    }
    //!\}

    /*!\name Invocation
     * \{
     */
    /*!\brief Computes the pairwise sequence alignment for the given range over indexed sequence pairs.
     * \tparam indexed_sequence_pairs_t The type of indexed_sequence_pairs; must model
     *                                  seqan3::detail::indexed_sequence_pair_range.
     * \tparam callback_t The type of the callback function that is called with the alignment result; must model
    *                    std::invocable with seqan3::alignment_result as argument.
     *
     * \param[in] indexed_sequence_pairs A range over indexed sequence pairs to be aligned.
     * \param[in] callback The callback function to be invoked with each computed alignment result.
     *
     * \throws std::bad_alloc during allocation of the alignment matrices or
     *         seqan3::invalid_alignment_configuration if an invalid configuration for the given sequences is detected.
     *
     * \details
     *
     * Uses the standard dynamic programming algorithm to compute the pairwise sequence alignment for each
     * sequence pair. The space and runtime complexities depend on the selected configurations (see below).
     * For every computed alignment the given callback is invoked with the respective alignment result.
     *
     * ### Exception
     *
     * Strong exception guarantee. Might throw std::bad_alloc or seqan3::invalid_alignment_configuration.
     *
     * ### Thread-safety
     *
     * Calls to this functions in a concurrent environment are not thread safe. Instead use a copy of the alignment
     * algorithm type.
     *
     * ### Complexity
     *
     * The following table lists the runtime and space complexities for the banded and unbanded algorithm dependent
     * on the configured seqan3::align_cfg::result per sequence pair.
     * Let `n` be the length of the first sequence, `m` be the length of the second sequence and `k` be the size of
     * the band.
     *
     * |                        | unbanded         | banded            |
     * |:----------------------:|:----------------:|:-----------------:|
     * |runtime                 |\f$ O(n*m) \f$    |\f$ O(n*k) \f$     |
     * |space (score only)      |\f$ O(m) \f$      |\f$ O(k) \f$       |
     * |space (end positions)   |\f$ O(m) \f$      |\f$ O(k) \f$       |
     * |space (begin positions) |\f$ O(n*m) \f$    |\f$ O(n*k) \f$     |
     * |space (alignment)       |\f$ O(n*m) \f$    |\f$ O(n*k) \f$     |
     */
    template <indexed_sequence_pair_range indexed_sequence_pairs_t, typename callback_t>
    //!\cond
        requires std::invocable<callback_t, alignment_result_type>
    //!\endcond
    void operator()(indexed_sequence_pairs_t && indexed_sequence_pairs, callback_t && callback)
    {
        using result_value_t = typename alignment_result_value_type_accessor<alignment_result_type>::type;
        using std::get;

        for (auto && [sequence_pair, idx] : indexed_sequence_pairs)
        {
            result_value_t res{};
            res.id = idx;
            res.score = compute_matrix(get<0>(sequence_pair), get<1>(sequence_pair));
            callback(alignment_result_type{res});
        }
    }

protected:
    /*!\brief Compute the actual alignment.
     * \tparam sequence1_t The type of the first sequence; must model std::ranges::forward_range.
     * \tparam sequence2_t The type of the second sequence; must model std::ranges::forward_range.
     *
     * \param[in] sequence1 The first sequence to compute the alignment for.
     * \param[in] sequence2 The second sequence to compute the alignment for.
     */
    template <std::ranges::forward_range sequence1_t, std::ranges::forward_range sequence2_t>
    int32_t compute_matrix(sequence1_t & sequence1, sequence2_t & sequence2)
    {
        // ---------------------------------------------------------------------
        // Initialisation phase: allocate memory and initialise first column.
        // ---------------------------------------------------------------------

        thread_local score_matrix_single_column<int32_t> local_matrix{};
        local_matrix.reset_matrix(std::forward<sequence1_t>(sequence1), std::forward<sequence2_t>(sequence2));
        auto matrix_iter = local_matrix.begin();
        auto alignment_column = *matrix_iter;

        // Initialise the first column.
        *alignment_column.begin() = this->initialise_origin_cell();
        for (auto cell : alignment_column | seqan3::views::drop(1))
            cell = this->initialise_first_column_cell(cell);

        // ---------------------------------------------------------------------
        // Recursion phase: compute column-wise the alignment matrix.
        // ---------------------------------------------------------------------

        // Compute remaining matrix
        ++matrix_iter; // Move to next column.
        for (auto sequence1_iter = sequence1.begin();
             sequence1_iter != sequence1.end();
             ++sequence1_iter, ++matrix_iter)
        {
            alignment_column = *matrix_iter;
            auto alignment_column_iter = alignment_column.begin();

            // Initialise first cell of current column.
            auto cell = *alignment_column_iter;
            typename traits_type::score_type diagonal = cell.optimal_score();
            *alignment_column_iter = this->initialise_first_row_cell(cell);

            ++alignment_column_iter;  // Move to next cell in column.
            for (auto sequence2_iter = sequence2.begin();
                 sequence2_iter != sequence2.end();
                 ++sequence2_iter, ++alignment_column_iter)
            {
                auto cell = *alignment_column_iter;
                typename traits_type::score_type next_diagonal = cell.optimal_score();
                *alignment_column_iter = this->compute_inner_cell(diagonal,
                                                                  cell,
                                                                  m_scoring_scheme.score(*sequence1_iter,
                                                                                         *sequence2_iter));
                diagonal = next_diagonal;
            }
        }

        // ---------------------------------------------------------------------
        // Wrap up phase: track score of last column
        // ---------------------------------------------------------------------

        return alignment_column[std::ranges::distance(sequence2)].optimal_score();
    }
};

} // namespace seqan3::detail
