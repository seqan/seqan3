// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::detail::policy_alignment_matrix.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <tuple>

#include <seqan3/alignment/exception.hpp>
#include <seqan3/alignment/matrix/detail/coordinate_matrix.hpp>
#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/core/detail/template_inspection.hpp>

namespace seqan3::detail
{

/*!\brief A policy that provides a common interface to acquire the correct alignment matrices.
 * \ingroup alignment_pairwise
 *
 * \tparam traits_t The alignment configuration traits type; must be an instance of
 *                  seqan3::detail::alignment_configuration_traits.
 * \tparam alignment_matrix_t The type of the alignment matrix for this alignment configuration
 *                            [see requirements below].
 *
 * \details
 *
 * The alignment matrix must be a matrix type that is compatible with the configured alignment algorithm. It must offer
 * a resize member function that takes a seqan3::detail::column_index_type and seqan3::detail::row_index_type and an
 * additional parameter to initialise the allocated matrix memory.
 */
template <typename traits_t, typename alignment_matrix_t>
    requires (is_type_specialisation_of_v<traits_t, alignment_configuration_traits> &&
              requires (alignment_matrix_t & matrix, typename traits_t::score_type const initial_score)
              {
                  { matrix.resize(column_index_type{size_t{}}, row_index_type{size_t{}}, initial_score) };
              })
class policy_alignment_matrix
{
protected:
    //!\brief The configured score type.
    using score_type = typename traits_t::score_type;
    //!\brief The configured matrix index type to store the coordinates.
    using matrix_index_type = typename traits_t::matrix_index_type;

    //!\brief The selected lower diagonal.
    int32_t lower_diagonal{};
    //!\brief The selected upper diagonal.
    int32_t upper_diagonal{};
    //!\brief A flag indicating whether the final gaps in the last column are free.
    bool last_column_is_free{};
    //!\brief A flag indicating whether the final gaps in the last row are free.
    bool last_row_is_free{};

    /*!\name Constructors, destructor and assignment
     * \{
     */
    policy_alignment_matrix() = default;                                            //!< Defaulted.
    policy_alignment_matrix(policy_alignment_matrix const &) = default;             //!< Defaulted.
    policy_alignment_matrix(policy_alignment_matrix &&) = default;                  //!< Defaulted.
    policy_alignment_matrix & operator=(policy_alignment_matrix const &) = default; //!< Defaulted.
    policy_alignment_matrix & operator=(policy_alignment_matrix &&) = default;      //!< Defaulted.
    ~policy_alignment_matrix() = default;                                           //!< Defaulted.

    /*!\brief Constructs and initialises the algorithm using the alignment configuration.
     * \tparam alignment_configuration_t The type of the alignment configuration; must be an instance of
     *                                   seqan3::configuration.
     *
     * \param[in] config The configuration passed into the algorithm.
     *
     * \details
     *
     * Initialises the members for the lower and upper diagonal. These members are only used if the banded alignment
     * is computed.
     *
     * \throws seqan3::invalid_alignment_configuration if the given band settings are invalid.
     */
    template <typename alignment_configuration_t>
        requires (is_type_specialisation_of_v<alignment_configuration_t, configuration>)
    policy_alignment_matrix(alignment_configuration_t const & config)
    {
        using seqan3::get;

        auto band = config.get_or(seqan3::align_cfg::band_fixed_size{});

        lower_diagonal = band.lower_diagonal;
        upper_diagonal = band.upper_diagonal;

        bool invalid_band = upper_diagonal < lower_diagonal;
        std::string error_cause = (invalid_band) ? " The upper diagonal is smaller than the lower diagonal." : "";

        if constexpr (traits_t::is_global)
        {
            auto method_global_config = get<seqan3::align_cfg::method_global>(config);

            bool first_row_is_free = method_global_config.free_end_gaps_sequence1_leading;
            bool first_column_is_free = method_global_config.free_end_gaps_sequence2_leading;

            last_row_is_free = method_global_config.free_end_gaps_sequence1_trailing;
            last_column_is_free = method_global_config.free_end_gaps_sequence2_trailing;
            // band starts in first column without free gaps or band starts in first row without free gaps.
            invalid_band |= (upper_diagonal < 0 && !first_column_is_free) || (lower_diagonal > 0 && !first_row_is_free);
            error_cause += " The band starts in a region without free gaps.";
        }

        if (invalid_band)
            throw invalid_alignment_configuration{"The selected band [" + std::to_string(lower_diagonal) + ":"
                                                  + std::to_string(upper_diagonal)
                                                  + "] cannot be used with the current "
                                                    "alignment configuration:"
                                                  + error_cause};
    }
    //!\}

    /*!\brief Acquires a new thread local alignment and index matrix for the given sequence sizes.
     *
     * \param[in] sequence1_size The size of the first sequence.
     * \param[in] sequence2_size The size of the second sequence.
     * \param[in] initial_score The initial score used for the acquired alignment matrix.
     *
     * \returns A std::tuple storing lvalue references to the thread local alignment and index matrix.
     *
     * \details
     *
     * Acquires a thread local alignment and index matrix. Initialises the matrices with the given
     * sequence sizes and the initial score value. In the banded alignment, the alignment matrix is reduced to
     * the column count times the band size.
     *
     * ### Exception
     *
     * Might throw std::bad_alloc if the requested matrix size exceeds the available memory or
     * seqan3::invalid_alignment_configuration if the band does not allow a valid computation of the
     * configured alignment.
     *
     * \throws std::bad_alloc or seqan3::invalid_alignment_configuration
     */
    auto acquire_matrices(size_t const sequence1_size,
                          size_t const sequence2_size,
                          score_type initial_score = score_type{}) const
    {
        assert(sequence1_size < static_cast<uint64_t>(std::numeric_limits<int64_t>::max()));
        assert(sequence2_size < static_cast<uint64_t>(std::numeric_limits<int64_t>::max()));

        if constexpr (traits_t::is_banded)
            check_valid_band_configuration(sequence1_size, sequence2_size);

        static thread_local alignment_matrix_t alignment_matrix{};
        static thread_local coordinate_matrix<matrix_index_type> index_matrix{};

        // Increase dimension by one for the initialisation of the matrix.
        size_t const column_count = sequence1_size + 1;
        size_t row_count = sequence2_size + 1;

        index_matrix.resize(column_index_type{column_count}, row_index_type{row_count});

        if constexpr (traits_t::is_banded)
        {
            assert(upper_diagonal - lower_diagonal + 1 > 0); // Band size is a positive integer.
            // Allocate one more cell to compute the last cell of the band with standard recursion function.
            row_count = std::min<int64_t>(upper_diagonal - lower_diagonal + 2, row_count);
        }

        alignment_matrix.resize(column_index_type{column_count}, row_index_type{row_count}, initial_score);

        return std::tie(alignment_matrix, index_matrix);
    }

    /*!\brief Checks whether the band is valid for the given sequence sizes.
     *
     * \param[in] sequence1_size The size of the first sequence.
     * \param[in] sequence2_size The size of the second sequence.
     *
     * \throws seqan3::invalid_alignment_configuration if the band is invalid for the given sequence sizes and the
     *         alignment configuration.
     */
    void check_valid_band_configuration(size_t const sequence1_size, size_t const sequence2_size) const
    {
        bool const upper_diagonal_ends_before_last_cell = (upper_diagonal + sequence2_size) < sequence1_size;
        bool const lower_diagonal_ends_behind_last_cell = (-lower_diagonal + sequence1_size) < sequence2_size;

        bool invalid_band = false;
        std::string error_cause{};

        if constexpr (traits_t::is_global)
        {
            // band ends in last column without free gaps or band ends in last row without free gaps.
            invalid_band |= (lower_diagonal_ends_behind_last_cell && !last_column_is_free)
                         || (upper_diagonal_ends_before_last_cell && !last_row_is_free);
            error_cause = "The band ends in a region without free gaps.";
        }

        if (invalid_band)
            throw invalid_alignment_configuration{"The selected band [" + std::to_string(lower_diagonal) + ":"
                                                  + std::to_string(upper_diagonal)
                                                  + "] cannot be used with the current "
                                                    "alignment configuration: "
                                                  + error_cause};
    }
};
} // namespace seqan3::detail
