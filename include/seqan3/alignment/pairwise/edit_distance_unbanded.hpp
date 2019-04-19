// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides a pairwise alignment algorithm for edit distance but without band.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <bitset>
#include <utility>

#include <range/v3/algorithm/copy.hpp>

#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/matrix/alignment_coordinate.hpp>
#include <seqan3/alignment/matrix/alignment_trace_algorithms.hpp>
#include <seqan3/alignment/matrix/edit_distance_score_matrix_full.hpp>
#include <seqan3/alignment/matrix/edit_distance_trace_matrix_full.hpp>
#include <seqan3/alignment/pairwise/align_result_selector.hpp>
#include <seqan3/alignment/pairwise/alignment_result.hpp>
#include <seqan3/alignment/pairwise/edit_distance_fwd.hpp>
#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/bit_manipulation.hpp>
#include <seqan3/range/shortcuts.hpp>
#include <seqan3/std/ranges>

namespace seqan3::detail
{

/*!\brief This calculates an alignment using the edit distance and without a band.
 * \ingroup pairwise_alignment
 * \tparam database_t     \copydoc pairwise_alignment_edit_distance_unbanded::database_type
 * \tparam query_t        \copydoc pairwise_alignment_edit_distance_unbanded::query_type
 * \tparam align_config_t The configuration type; must be of type seqan3::configuration.
 */
template <std::ranges::ViewableRange database_t,
          std::ranges::ViewableRange query_t,
          typename align_config_t,
          EditDistanceTrait traits_t>
class pairwise_alignment_edit_distance_unbanded
{
    //!\brief The horizontal/database sequence.
    database_t database;
    //!\brief The vertical/query sequence.
    query_t query;
    //!\brief The configuration.
    align_config_t config;

public:
    //!\brief The type of one machine word.
    using word_type = typename std::remove_reference_t<traits_t>::word_type;
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
private:
    //!\brief The type of an iterator of the database sequence.
    using database_iterator = std::ranges::iterator_t<database_type>;
    //!\brief The alphabet type of the query sequence.
    using query_alphabet_type = std::remove_reference_t<reference_t<query_type>>;
    //!\brief The intermediate result type of the execution of this function object.
    using result_value_type = typename align_result_selector<database_type, query_type, align_config_type>::type;

    //!\brief When true the computation will use the ukkonen trick with the last active cell and bounds the error to config.max_errors.
    static constexpr bool use_max_errors = align_config_type::template exists<align_cfg::max_error>();
    //!\brief Whether the alignment is a semi-global alignment or not.
    static constexpr bool is_semi_global = traits_t::is_semi_global_type::value;
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
    static constexpr bool compute_score_matrix = compute_front_coordinate || compute_sequence_alignment;
    //!\brief Whether the alignment configuration indicates to compute and/or store the trace matrix.
    static constexpr bool compute_trace_matrix = compute_front_coordinate || compute_sequence_alignment;
    //!\brief Whether the alignment configuration indicates to compute and/or store the score or trace matrix.
    static constexpr bool compute_matrix = compute_score_matrix || compute_trace_matrix;

    //!\brief How to pre-initialise hp.
    static constexpr word_type hp0 = is_global ? 1u : 0u;
    //!\brief How to pre-initialise hn.
    static constexpr word_type hn0 = 0u;
    //!\brief How to pre-initialise vp.
    static constexpr word_type vp0 = ~word_type{0u};
    //!\brief How to pre-initialise vn.
    static constexpr word_type vn0 = 0u;

    //!\brief The score of the current column.
    score_type _score{};
    //!\brief The mask with a bit set at the position where the score change.
    //!\details If #use_max_errors is true this corresponds to the last active cell.
    word_type score_mask{0};
    //!\brief The machine words which stores the positive vertical differences.
    std::vector<word_type> vp{};
    //!\brief The machine words which stores the negative vertical differences.
    std::vector<word_type> vn{};
    //!\brief The machine words which stores the positive horizontal differences.
    std::vector<word_type> hp{};
    //!\brief The machine words which stores if trace_directions::diagonal is true.
    std::vector<word_type> db{};
    //!\brief The machine words which translate a letter of the query into a bit mask.
    //!\details Each bit position which is true (= 1) corresponds to a match of a letter in the query at this position.
    std::vector<word_type> bit_masks{};
    /*!\brief The best score of the alignment in the last row
     * (if is_semi_global = true) or the last entry in the
     * score matrix (if is_global = true).
     */
    score_type _best_score{};
    /*!\brief In which column the best score of the alignment is
     * located. Will only be tracked if is_semi_global is true.
     *
     * \details
     *
     * If is_global is true this is always at the last entry in
     * the score matrix, i.e. at position (`|query|`,
     * `|database|`).
     */
    database_iterator _best_score_col{};

    /*!\name Only used when use_max_errors is true
     * \{
     */
    //!\brief Which score value is considered as a hit?
    score_type max_errors{255};
    //!\brief The block containing the last active cell.
    size_t last_block{0};
    //!\brief A mask with a bit set on the position of the last row.
    word_type last_score_mask{};
    //!\}

    //!\brief The current position in the database.
    database_iterator database_it{};
    //!\brief The end position of the database.
    database_iterator database_it_end{};

    using score_matrix_type = edit_distance_score_matrix_full<word_type, score_type, is_semi_global, use_max_errors>;
    score_matrix_type _score_matrix{};

    using trace_matrix_type = edit_distance_trace_matrix_full<word_type, is_semi_global, use_max_errors>;
    trace_matrix_type _trace_matrix{};

    //!\brief Add a computation step
    void add_state()
    {
        if constexpr(!use_max_errors)
        {
            _score_matrix.add_column(vp, vn);
            _trace_matrix.add_column(hp, db, vp);
        }

        if constexpr(use_max_errors)
        {
            auto max_rows = _score_matrix.max_rows(score_mask, last_block, _score, max_errors);
            _score_matrix.add_column(vp, vn, max_rows);
            _trace_matrix.add_column(hp, db, vp, max_rows);
        }
    }

public:

    /*!\name Constructors, destructor and assignment
     * \{
     */
     //!\brief The class template parameter may resolve to an lvalue reference which prohibits default constructibility.
     pairwise_alignment_edit_distance_unbanded() = delete;
     //!\brief Defaulted
     pairwise_alignment_edit_distance_unbanded(pairwise_alignment_edit_distance_unbanded const &) = default;
     pairwise_alignment_edit_distance_unbanded(pairwise_alignment_edit_distance_unbanded &&) = default; //!< Defaulted
     //!\brief Defaulted
     pairwise_alignment_edit_distance_unbanded & operator=(pairwise_alignment_edit_distance_unbanded const &) = default;
     //!\brief Defaulted
     pairwise_alignment_edit_distance_unbanded & operator=(pairwise_alignment_edit_distance_unbanded &&) = default;
     ~pairwise_alignment_edit_distance_unbanded() = default;                                           //!< Defaulted

    /*!\brief Constructor
     * \param[in] _database \copydoc database
     * \param[in] _query    \copydoc query
     * \param[in] _config   \copydoc config
     * \param[in] _traits   The traits object. Only the type information will be used.
     */
    pairwise_alignment_edit_distance_unbanded(database_t && _database,
                                              query_t && _query,
                                              align_config_t _config,
                                              traits_t const & SEQAN3_DOXYGEN_ONLY(_traits) = traits_t{}) :
        database{std::forward<database_t>(_database)},
        query{std::forward<query_t>(_query)},
        config{std::forward<align_config_t>(_config)},
        _score{static_cast<score_type>(query.size())},
        _best_score{static_cast<score_type>(query.size())},
        _best_score_col{ranges::end(database)},
        database_it{ranges::begin(database)},
        database_it_end{ranges::end(database)},
        _score_matrix{query.size() + 1u},
        _trace_matrix{query.size() + 1u}
    {
        static constexpr size_t alphabet_size_ = alphabet_size<query_alphabet_type>;
        _score_matrix.reserve(database.size() + 1u);
        _trace_matrix.reserve(database.size() + 1u);

        if constexpr (use_max_errors)
        {
            max_errors = get<align_cfg::max_error>(config).value;
            assert(max_errors >= score_type{0});
        }

        size_t block_count = (query.size() - 1 + word_size) / word_size;
        score_mask = (word_type)1 << ((query.size() - 1 + word_size) % word_size);
        last_score_mask = score_mask;
        last_block = block_count - 1;

        if constexpr(use_max_errors)
        {
            // local_max_errors either stores the maximal number of _score (me.max_errors) or the needle size minus one.
            // It is used for the mask computation and setting the initial score (the minus one is there because of the Ukkonen trick).
            size_t local_max_errors = std::min<size_t>(max_errors, query.size() - 1);
            score_mask = (word_type)1 << (local_max_errors % word_size);
            last_block = std::min(local_max_errors / word_size, block_count - 1);
            _score = local_max_errors + 1;
        }

        vp.resize(block_count, vp0);
        vn.resize(block_count, vn0);
        hp.resize(block_count, 0u);
        db.resize(block_count, 0u);
        bit_masks.resize((alphabet_size_ + 1) * block_count, 0);

        // encoding the letters as bit-vectors
        for (size_t j = 0; j < query.size(); j++)
        {
            size_t i = block_count * to_rank(query[j]) + j / word_size;
            bit_masks[i] |= (word_type)1 << (j % word_size);
        }

        add_state();
    }
    //!\}

private:

    //!\brief One compute step in one column.
    template <typename carry_type>
    void compute_step(word_type b, word_type & d0, word_type & hp, word_type & hn, word_type & vp, word_type & vn,
                      carry_type carry_d0, carry_type carry_hp, carry_type carry_hn)
    {
        word_type x, t;
        assert(carry_d0 <= 1u);
        assert(carry_hp <= 1u);
        assert(carry_hn <= 1u);

        x = b | vn;
        t = vp + (x & vp) + carry_d0;

        d0 = (t ^ vp) | x;
        hn = vp & d0;
        hp = vn | ~(vp | d0);

        if constexpr(std::Same<carry_type, word_type &>)
            carry_d0 = (carry_d0 != 0u) ? t <= vp : t < vp;

        x = (hp << 1u) | carry_hp;
        vn = x & d0;
        vp = (hn << 1u) | ~(x | d0) | carry_hn;

        if constexpr(std::Same<carry_type, word_type &>)
        {
            carry_hp = hp >> (word_size - 1u);
            carry_hn = hn >> (word_size - 1u);
        }
    }

    //!\brief Increase or decrease the score.
    void advance_score(word_type P, word_type N, word_type mask)
    {
        if ((P & mask) != (word_type)0)
            _score++;
        else if ((N & mask) != (word_type)0)
            _score--;
    }

    //!\brief Returns true if the current active cell is within the last row.
    bool is_last_active_cell_within_last_row()
    {
        return (score_mask == last_score_mask) && (last_block == vp.size() - 1);
    }

    //!\brief Update the current best known score if the current score is better.
    void update_best_score()
    {
        if constexpr(is_global)
            _best_score = _score;

        if constexpr(is_semi_global)
        {
            // we have to make sure that update_best_score is only called after a score update within the
            // last row.
            assert(is_last_active_cell_within_last_row());

            _best_score_col = (_score <= _best_score) ? database_it : _best_score_col;
            _best_score     = (_score <= _best_score) ? _score : _best_score;
        }
    }

    //!\brief Decrement the last active cell position.
    bool prev_last_active_cell()
    {
        score_mask >>= 1;
        if (score_mask != (word_type)0)
            return true;

        if constexpr (is_global)
            if (last_block == 0u)  // [[unlikely]]
                return false;

        last_block--;

        score_mask = (word_type)1 << (word_size - 1);
        return true;
    }

    //!\brief Increment the last active cell position.
    void next_last_active_cell()
    {
        score_mask <<= 1;
        if (score_mask)
            return;

        score_mask = 1;
        last_block++;
    }

    //!\brief Use the ukkonen trick and update the last active cell.
    bool update_last_active_cell()
    {
        // updating the last active cell
        while (!(_score <= max_errors))
        {
            advance_score(vn[last_block], vp[last_block], score_mask);
            if (!prev_last_active_cell())
                break;
        }

        if (is_last_active_cell_within_last_row())
            return on_hit();
        else
        {
            next_last_active_cell();
            advance_score(vp[last_block], vn[last_block], score_mask);
        }

        return false;
    }

    //!\brief Will be called if a hit was found (e.g., score <= max_errors).
    bool on_hit()
    {
        assert(_score <= max_errors);

        if constexpr(is_semi_global)
            update_best_score();

        // TODO: call external on_hit functor

        return false;
    }

    //!\brief Pattern is small enough that it fits into one machine word. Use
    //!faster computation with less overhead.
    inline bool small_patterns();

    //!\brief Pattern is larger than one machine word. Use overflow aware computation.
    inline bool large_patterns();

    //!\brief Compute the alignment.
    void _compute()
    {
        // limit search width for prefix search
        if constexpr(use_max_errors && is_global)
        {
            size_t max_length = query.size() + max_errors + 1;
            size_t haystack_length = std::min(database.size(), max_length);
            database_it_end -= database.size() - haystack_length;
        }

        // distinguish between the version for needles not longer than
        // one machine word and the version for longer needles
        if (vp.size() <= 1)
            small_patterns();
        else
            large_patterns();

        if constexpr(is_global)
            update_best_score();
    }

public:

    /*!\brief Generic invocable interface.
     * \param[in]     idx The index of the currently processed sequence pair.
     * \returns A reference to the filled alignment result.
     */
    alignment_result<result_value_type> operator()(size_t const idx)
    {
        _compute();
        result_value_type res_vt{};
        res_vt.id = idx;
        if constexpr (compute_score)
        {
            res_vt.score = score().value_or(matrix_inf<score_type>);
        }

        if constexpr (compute_back_coordinate)
        {
            res_vt.back_coordinate = back_coordinate();
        }

        if constexpr (compute_front_coordinate)
        {
            if (is_valid())
                res_vt.front_coordinate = alignment_front_coordinate(trace_matrix(), res_vt.back_coordinate);
            else
                res_vt.front_coordinate = invalid_coordinate();
        }

        if constexpr (compute_sequence_alignment)
        {
            if (is_valid())
            {
                using alignment_t = decltype(res_vt.alignment);
                res_vt.alignment = alignment_trace<alignment_t>(database,
                                                                query,
                                                                trace_matrix(),
                                                                res_vt.back_coordinate,
                                                                res_vt.front_coordinate);
            }
        }
        return alignment_result<result_value_type>{std::move(res_vt)};
    }

    //!\brief Return the score of the alignment.
    std::optional<score_type> score() const noexcept
    {
        static_assert(compute_score, "score() can only be computed if you specify the result type within "
                                     "your alignment config.");
        if (!is_valid())
            return std::nullopt;

        return -_best_score;
    }

    //!\brief Return the score matrix of the alignment.
    score_matrix_type const & score_matrix() const noexcept
    {
        static_assert(compute_score_matrix, "score_matrix() can only be computed if you specify the result type within "
                                            "your alignment config.");
        return _score_matrix;
    }

    //!\brief Return the trace matrix of the alignment.
    trace_matrix_type const & trace_matrix() const noexcept
    {
        static_assert(compute_trace_matrix, "trace_matrix() can only be computed if you specify the result type within "
                                            "your alignment config.");
        return _trace_matrix;
    }

    //!\brief Return the begin position of the alignment
    alignment_coordinate front_coordinate() const noexcept
    {
        static_assert(compute_front_coordinate, "front_coordinate() can only be computed if you specify the result type "
                                                "within your alignment config.");
        if (!is_valid())
            return invalid_coordinate();

        alignment_coordinate back = back_coordinate();
        return alignment_front_coordinate(trace_matrix(), back);
    }

    //!\brief Return the end position of the alignment
    alignment_coordinate back_coordinate() const noexcept
    {
        static_assert(compute_back_coordinate, "back_coordinate() can only be computed if you specify the result type "
                                               "within your alignment config.");
        if (!is_valid())
            return invalid_coordinate();

        size_t col = database.size() - 1;
        if constexpr(is_semi_global)
            col = std::ranges::distance(begin(database), _best_score_col);

        return {column_index_type{col}, row_index_type{query.size() - 1}};
    }

    //!\brief Returns true if the computation produced a valid alignment.
    bool is_valid() const noexcept
    {
        // This condition is obvious, because _best_score_col will only be set to a new position if the last active cell
        // is within the last row AND score <= max_errors. And _best_score_col can't reach database_it_end again.
        if constexpr(use_max_errors && is_semi_global)
            return _best_score_col != database_it_end;

        // This condition uses the observation that after each computation of a column, _score has either the initial
        // value of the first row (i.e. the entire column consist of INF's), has the value _score = max_errors + 1 (there
        // exists a cell within the column that has value <= max_errors, but is not on the last row) or _score <=
        // max_errors (the score of the last active cell is <= max_errors)
        if constexpr(use_max_errors && is_global)
            return _best_score <= max_errors;

        // When not using max_errors there is always a valid alignment, because the last row will always be updated and
        // with it the score.
        return true;
    }

    //!\brief Returns an invalid_coordinate for this alignment.
    alignment_coordinate invalid_coordinate() const noexcept
    {
        return {column_index_type{database.size()}, row_index_type{query.size()}};
    }

    //!\brief Return the alignment, i.e. the actual base pair matching.
    auto alignment() const noexcept
    {
        static_assert(compute_sequence_alignment, "alignment() can only be computed if you specify the result type "
                                                  "within your alignment config.");
        using alignment_t = decltype(result_value_type{}.alignment);

        if (!is_valid())
            return alignment_t{};

        return alignment_trace<alignment_t>(database, query, trace_matrix(), back_coordinate(), front_coordinate());
    }
};

template <typename database_t, typename query_t, typename align_config_t, typename traits_t>
bool pairwise_alignment_edit_distance_unbanded<database_t, query_t, align_config_t, traits_t>::small_patterns()
{
    // computing the blocks
    while (database_it != database_it_end)
    {
        word_type d0, hn;

        word_type const b = bit_masks[to_rank((query_alphabet_type) *database_it)];
        compute_step<word_type>(b, d0, hp[0], hn, vp[0], vn[0], 0u, hp0, 0u);
        advance_score(hp[0], hn, score_mask);
        db[0] = ~(b ^ d0);

        // semi-global without max_errors guarantees that the score stays within the last row
        if constexpr(is_semi_global && !use_max_errors)
            update_best_score();

        if constexpr(use_max_errors)
        {
            // updating the last active cell
            if (update_last_active_cell())
            {
                add_state();
                ++database_it;
                return true;
            }
        }

        add_state();
        ++database_it;
    }

    return false;
}

template <typename database_t, typename query_t, typename align_config_t, typename traits_t>
bool pairwise_alignment_edit_distance_unbanded<database_t, query_t, align_config_t, traits_t>::large_patterns()
{
    while (database_it != database_it_end)
    {
        word_type d0, hn;
        word_type carry_d0{0u}, carry_hp{hp0}, carry_hn{0u};
        size_t block_offset = vp.size() * to_rank((query_alphabet_type) *database_it);

        // computing the necessary blocks, carries between blocks following one another are stored
        for (size_t current_block = 0; current_block <= last_block; current_block++)
        {
            word_type const b = bit_masks[block_offset + current_block];
            compute_step<word_type &>(b, d0, hp[current_block], hn, vp[current_block], vn[current_block],
                                      carry_d0, carry_hp, carry_hn);
            db[current_block] = ~(b ^ d0);
        }
        advance_score(hp[last_block], hn, score_mask);

        // semi-global without max_errors guarantees that the score stays within the last row
        if constexpr(is_semi_global && !use_max_errors)
            update_best_score();

        if constexpr(use_max_errors)
        {
            // if the active cell is the last of it's block, one additional block has to be calculated
            bool additional_block = score_mask >> (word_size - 1);
            if (last_block + 1 == vp.size())
                additional_block = false;

            if (additional_block)
            {
                size_t const current_block = last_block + 1u;
                word_type const b = bit_masks[block_offset + current_block];
                // this might not be necessary, but carry_d0 = 1u might have an influence on the result of vn and vp.
                vp[current_block] = vp0;
                vn[current_block] = vn0;
                compute_step<word_type>(b, d0, hp[current_block], hn, vp[current_block], vn[current_block],
                                        carry_d0, carry_hp, carry_hn);
                db[current_block] = ~(b ^ d0);
            }

            // updating the last active cell
            if (update_last_active_cell())
            {
                add_state();
                ++database_it;
                return true;
            }
        }

        add_state();
        ++database_it;
    }

    return false;
}

/*!\name Type deduction guides
 * \relates seqan3::detail::pairwise_alignment_edit_distance_unbanded
 * \{
 */

//!\brief Deduce the type from the provided arguments.
template<typename database_t, typename query_t, typename config_t>
pairwise_alignment_edit_distance_unbanded(database_t && database, query_t && query, config_t config)
    -> pairwise_alignment_edit_distance_unbanded<database_t, query_t, config_t>;

//!\brief Deduce the type from the provided arguments.
template<typename database_t, typename query_t, typename config_t, typename traits_t>
pairwise_alignment_edit_distance_unbanded(database_t && database, query_t && query, config_t config, traits_t)
    -> pairwise_alignment_edit_distance_unbanded<database_t, query_t, config_t, traits_t>;
//!\}

} // namespace seqan3::detail
