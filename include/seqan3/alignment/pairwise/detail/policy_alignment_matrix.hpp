// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::policy_alignment_matrix.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <tuple>

#include <seqan3/alignment/exception.hpp>
#include <seqan3/alignment/matrix/detail/coordinate_matrix.hpp>
#include <seqan3/alignment/pairwise/detail/type_traits.hpp>
#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>

namespace seqan3::detail
{

/*!\brief A policy that provides a common interface to acquire the correct alignment matrices.
 * \ingroup pairwise_alignment
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
//!\cond
    requires (is_type_specialisation_of_v<traits_t, alignment_configuration_traits> &&
              requires (alignment_matrix_t & matrix, typename traits_t::score_type const initial_score)
              {
                  { matrix.resize(column_index_type{size_t{}}, row_index_type{size_t{}}, initial_score) };
              })
//!\endcond
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

    /*!\name Constructors, destructor and assignment
     * \{
     */
    policy_alignment_matrix() = default; //!< Defaulted.
    policy_alignment_matrix(policy_alignment_matrix const &) = default; //!< Defaulted.
    policy_alignment_matrix(policy_alignment_matrix &&) = default; //!< Defaulted.
    policy_alignment_matrix & operator=(policy_alignment_matrix const &) = default; //!< Defaulted.
    policy_alignment_matrix & operator=(policy_alignment_matrix &&) = default; //!< Defaulted.
    ~policy_alignment_matrix() = default; //!< Defaulted.

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
    //!\cond
        requires (is_type_specialisation_of_v<alignment_configuration_t, configuration>)
    //!\endcond
    policy_alignment_matrix(alignment_configuration_t const & config)
    {
        using seqan3::get;

        auto band = config.get_or(seqan3::align_cfg::band_fixed_size{});

        lower_diagonal = band.lower_diagonal.get();
        upper_diagonal = band.upper_diagonal.get();

        if (upper_diagonal < lower_diagonal || upper_diagonal < 0 || lower_diagonal > 0)
            throw invalid_alignment_configuration{"The selected band [" + std::to_string(lower_diagonal) + ":" +
                                                  std::to_string(upper_diagonal) + "] is not supported because it "
                                                  "does not completely encloses both sequences."};
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
     * Might throw std::bad_alloc if the requested matrix size exceeds the available memory.
     *
     * \throws std::bad_alloc
     */
    auto acquire_matrices(size_t const sequence1_size,
                          size_t const sequence2_size,
                          score_type initial_score = score_type{}) const
    {
        assert(sequence1_size < static_cast<uint64_t>(std::numeric_limits<int64_t>::max()));
        assert(sequence2_size < static_cast<uint64_t>(std::numeric_limits<int64_t>::max()));

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
};
} // namespace seqan3::detail
