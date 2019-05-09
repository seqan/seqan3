// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Contains seqan3::detail::edit_distance_score_matrix_full.
 */

#pragma once

#include <bitset>

#include <seqan3/alignment/pairwise/edit_distance_fwd.hpp>
#include <seqan3/core/bit_manipulation.hpp>

namespace seqan3::detail
{

template <typename word_t, typename score_t, bool is_semi_global, bool use_max_errors>
class edit_distance_score_matrix_full
{
public:
    template <std::ranges::ViewableRange database_t,
              std::ranges::ViewableRange query_t,
              typename align_config_t,
              EditDistanceTrait traits_t>
    friend class pairwise_alignment_edit_distance_unbanded;

    /*!\name Constructors, destructor and assignment
     * \{
     */
     edit_distance_score_matrix_full() = default; //!< Defaulted
     edit_distance_score_matrix_full(edit_distance_score_matrix_full const &) = default; //!< Defaulted
     edit_distance_score_matrix_full(edit_distance_score_matrix_full &&) = default; //!< Defaulted
     edit_distance_score_matrix_full & operator=(edit_distance_score_matrix_full const &) = default; //!< Defaulted
     edit_distance_score_matrix_full & operator=(edit_distance_score_matrix_full &&) = default; //!< Defaulted
     ~edit_distance_score_matrix_full() = default; //!< Defaulted
    //!\}

    using word_type = word_t;

    static constexpr auto word_size = sizeof_bits<word_type>;

    using score_type = score_t;
    using entry_type = std::conditional_t<use_max_errors, std::optional<score_type>, score_type>;

    static constexpr std::optional<score_type> inf = std::nullopt;

protected:
    edit_distance_score_matrix_full(size_t const rows_size)
        : rows_size{rows_size}, columns{}
    {}

    void reserve(size_t const columns_size)
    {
        columns.reserve(columns_size);
    }

    template <typename score_type>
    static size_t max_rows(word_type const score_mask, unsigned const last_block,
                           score_type const score, score_type const max_errors)
    {
        size_t const offset = score_mask == 0u ? 0u : bit_scan_reverse(score_mask) + 1u;
        size_t const active_row = word_size * last_block + offset;
        return active_row + (score <= max_errors);
    }

    static score_type score_delta_within_word(word_type const & vp, word_type const & vn,
                                              uint8_t const offset) noexcept
    {
        using bitset = std::bitset<word_size>;
        assert(offset < word_size);

        score_type const p = bitset(vp)[offset] ? 1u : 0u;
        score_type const n = bitset(vn)[offset] ? 1u : 0u;

        return p - n;
    }

    static score_type score_delta_of_word(word_type const & vp, word_type const & vn) noexcept
    {
        using bitset = std::bitset<word_size>;

        score_type const p = bitset(vp).count();
        score_type const n = bitset(vn).count();

        return p - n;
    }

public:
    //!\copydoc seqan3::detail::Matrix::at
    entry_type at(size_t const row, size_t const col) const noexcept
    {
        assert(row < rows());
        assert(col < cols());

        column_type const & column = columns[col];
        if constexpr(use_max_errors)
            if (!(row < column.max_rows))
                return inf;

        score_type score = is_semi_global ? 0u : static_cast<score_type>(col);

        size_t current_row = 1u;
        size_t word_idx = 0u;

        for (; current_row + word_size <= row; ++word_idx, current_row += word_size)
            score += score_delta_of_word(column.vp[word_idx], column.vn[word_idx]);

        for (size_t offset = 0u; current_row + offset <= row; ++offset)
            score += score_delta_within_word(column.vp[word_idx], column.vn[word_idx], offset);

        return -score;
    }

    //!\copydoc seqan3::detail::Matrix::rows
    size_t rows() const noexcept
    {
        return rows_size;
    }

    //!\copydoc seqan3::detail::Matrix::cols
    size_t cols() const noexcept
    {
        return columns.size();
    }

protected:
    //!\brief If use_max_errors is true store these additional
    //!state information in state_type.
    struct max_errors_state
    {
        size_t max_rows;
    };

    //!\brief This information is needed to infer the score matrix.
    struct score_matrix_state
    {
        //!\copydoc pairwise_alignment_edit_distance_unbanded::vp
        std::vector<word_type> vp{};
        //!\copydoc pairwise_alignment_edit_distance_unbanded::vn
        std::vector<word_type> vn{};
    };

    //!\brief The state of one computation step.
    struct column_type :
        enable_state_t<true, score_matrix_state>,
        enable_state_t<use_max_errors, max_errors_state>
    {};

    void add_column(std::vector<word_type> vp, std::vector<word_type> vn)
        requires !use_max_errors
    {
        column_type column{};
        column.vp = std::move(vp);
        column.vn = std::move(vn);

        columns.push_back(std::move(column));
    }

    void add_column(std::vector<word_type> vp, std::vector<word_type> vn, size_t max_rows)
        requires use_max_errors
    {
        column_type column{};
        column.vp = std::move(vp);
        column.vn = std::move(vn);
        column.max_rows = max_rows;

        columns.push_back(std::move(column));
    }

private:
    size_t rows_size{};
    std::vector<column_type> columns{};
};

} // namespace seqan3::detail
