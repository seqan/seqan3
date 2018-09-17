// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

/*!\file
 * \brief Contains a pairwise alignment algorithm for edit distance but without band.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <bitset>
#include <utility>

#include <range/v3/algorithm/copy.hpp>
#include <range/v3/utility/iterator.hpp>

#include <seqan3/alignment/configuration/all.hpp>
#include <seqan3/alignment/matrix/alignment_score_matrix.hpp>
#include <seqan3/alignment/matrix/alignment_trace_matrix.hpp>
#include <seqan3/alignment/pairwise/align_result.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/add_enum_bitwise_operators.hpp>
#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/range/shortcuts.hpp>
#include <seqan3/std/ranges>

namespace seqan3::detail
{
//!\cond
template <typename align_config_t>
concept semi_global_config_concept = requires (align_config_t & cfg)
{
    requires get<align_cfg::id::sequence_ends>(cfg) == free_ends_at::seq1;
};

template <typename align_config_t>
concept global_config_concept = has_align_cfg_v<align_cfg::id::global, std::remove_reference_t<align_config_t>>;

template <typename align_config_t>
concept max_errors_concept = has_align_cfg_v<align_cfg::id::max_error, std::remove_reference_t<align_config_t>>;
//!\endcond

/*!\todo Document me
 * \ingroup pairwise
 */
template <typename traits_type>
concept edit_distance_trait_concept = requires
{
    typename std::remove_reference_t<traits_type>::word_type;
};

/*!\brief The default traits type for the edit distance algorithm.
 * \ingroup pairwise
 */
struct default_edit_distance_trait_type
{
    //!\brief The default word type.
    using word_type = uint64_t;
};

/*!\brief This calculates an alignment using the edit distance and without a band.
 * \ingroup pairwise
 * \tparam database_t     \copydoc pairwise_alignment_edit_distance_unbanded::database_type
 * \tparam query_t        \copydoc pairwise_alignment_edit_distance_unbanded::query_type
 * \tparam align_config_t The type of the alignment config.
 */
template <std::ranges::ViewableRange database_t,
          std::ranges::ViewableRange query_t,
          typename align_config_t,
          edit_distance_trait_concept traits_t = default_edit_distance_trait_type>
class pairwise_alignment_edit_distance_unbanded
{
    /*!\name Befriended classes
     * \{
     */
    //!\brief Befriend seqan3::detail::alignment_score_matrix<pairwise_alignment_edit_distance_unbanded>
    friend alignment_score_matrix<pairwise_alignment_edit_distance_unbanded>;
    //!\brief Befriend seqan3::detail::alignment_trace_matrix<pairwise_alignment_edit_distance_unbanded>
    friend alignment_trace_matrix<pairwise_alignment_edit_distance_unbanded>;
    //!\}

    //!\brief The horizontal/database sequence.
    database_t database;
    //!\brief The vertical/query sequence.
    query_t query;
    //!\brief The configuration.
    align_config_t config;

public:
    //!\brief The type of one machine word.
    using word_type = typename std::remove_reference_t<traits_t>::word_type;
    //!\brief The type of the score.
    using score_type = int;
    //!\brief The type of the database sequence.
    using database_type = std::remove_reference_t<database_t>;
    //!\brief The type of the query sequence.
    using query_type = std::remove_reference_t<query_t>;
    //!\brief The type of the score matrix.
    using score_matrix_type = detail::alignment_score_matrix<pairwise_alignment_edit_distance_unbanded>;
    //!\brief The type of the trace matrix.
    using trace_matrix_type = detail::alignment_trace_matrix<pairwise_alignment_edit_distance_unbanded>;

    //!\brief The size of one machine word.
    static constexpr uint8_t word_size = sizeof(word_type) * 8;

private:
    //!\brief The type of an iterator of the database sequence.
    using database_iterator = std::ranges::iterator_t<database_type>;
    //!\brief The alphabet type of the query sequence.
    using query_alphabet_type = std::remove_reference_t<reference_t<query_type>>;

    //TODO Make it dynamic.
    // using result_type = align_result<type_list<uint32_t, int>>;

    //!\brief When true the computation will use the ukkonen trick with the last active cell and bounds the error to config.max_errors.
    static constexpr bool use_max_errors = detail::max_errors_concept<align_config_t>;
    //!\brief Whether the alignment is a semi-global alignment or not.
    static constexpr bool is_semi_global = detail::semi_global_config_concept<align_config_t>;
    //!\brief Whether the alignment is a global alignment or not.
    static constexpr bool is_global = detail::global_config_concept<align_config_t> && !is_semi_global;

    //!\brief How to pre-initialize hp.
    static constexpr word_type hp0 = is_global ? 1 : 0;

    static_assert(8 * sizeof(word_type) <= 64, "we assume at most uint64_t as word_type");
    static_assert((is_global && !is_semi_global) || (!is_global && is_semi_global), "Either set global or semi-global");

    //!\brief The score of the current column.
    score_type _score{};
    //!\brief The mask with a bit set at the position where the score change.
    //!\details If #use_max_errors is true this corresponds to the last active cell.
    word_type score_mask{0};
    //!\brief The machine words which stores the positive vertical differences.
    std::vector<word_type> vp{};
    //!\brief The machine words which stores the negative vertical differences.
    std::vector<word_type> vn{};
    //!\brief The machine words which translate a letter of the query into a bit mask.
    //!\details Each bit position which is true (=1) corresponds to a match of a letter in the query at this position.
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
    size_t max_errors{255};
    //!\brief The block containing the last active cell.
    size_t last_block{0};
    //!\brief A mask with a bit set on the position of the last row.
    word_type last_score_mask{};
    //!\}

    //!\brief The current position in the database.
    database_iterator database_it{};
    //!\brief The end position of the database.
    database_iterator database_it_end{};

    //!\brief The state of one column computation.
    struct state_type
    {
        //!\copydoc pairwise_alignment_edit_distance_unbanded::vp
        std::vector<word_type> vp{};
        //!\copydoc pairwise_alignment_edit_distance_unbanded::vn
        std::vector<word_type> vn{};
    };

    //!\brief The collection of each computation step.
    std::vector<state_type> states{};

    //!\brief Add a computation step
    void add_state()
    {
        states.push_back(state_type{vp, vn});
    }

public:

    /*!\name Constructors, destructor and assignment
     * \{
     */
     pairwise_alignment_edit_distance_unbanded() = delete;
     pairwise_alignment_edit_distance_unbanded(pairwise_alignment_edit_distance_unbanded const &) = default;
     pairwise_alignment_edit_distance_unbanded(pairwise_alignment_edit_distance_unbanded &&) = default;
     pairwise_alignment_edit_distance_unbanded & operator=(pairwise_alignment_edit_distance_unbanded const &) = default;
     pairwise_alignment_edit_distance_unbanded & operator=(pairwise_alignment_edit_distance_unbanded &&) = default;

    /*!\brief Constructor
     * \param[in] _database \copydoc database
     * \param[in] _query    \copydoc query
     * \param[in] _config   \copydoc config
     */
    pairwise_alignment_edit_distance_unbanded(database_t && _database, query_t && _query, align_config_t _config) :
        database{std::forward<database_t>(_database)},
        query{std::forward<query_t>(_query)},
        config{std::forward<align_config_t>(_config)},
        _score{query.size()},
        _best_score{query.size()},
        _best_score_col{ranges::begin(database)},
        database_it{ranges::begin(database)},
        database_it_end{ranges::end(database)}
    {
        static constexpr std::size_t alphabet_size = alphabet_size_v<query_alphabet_type>;

        if constexpr(use_max_errors)
            max_errors = get<align_cfg::id::max_error>(config);

        size_t block_count = (query.size() - 1 + word_size) / word_size;
        score_mask = (word_type)1 << ((query.size() - 1 + word_size) % word_size);
        last_score_mask = score_mask;
        last_block = block_count - 1;

        if constexpr(use_max_errors)
        {
            // localMaxErrors either stores the maximal number of _score (me.max_errors) or the needle size minus one.
            // It is used for the mask computation and setting the initial score (the minus one is there because of the Ukkonen trick).
            size_t localMaxErrors = std::min<size_t>(max_errors, query.size() - 1);
            score_mask = (word_type)1 << (localMaxErrors % word_size);
            last_block = std::min(localMaxErrors / word_size, block_count - 1);
            _score = localMaxErrors + 1;
            _best_score = _score;
        }

        word_type vp0{~static_cast<word_type>(0)};
        word_type vn0{0};

        vp.resize(block_count, vp0);
        vn.resize(block_count, vn0);
        bit_masks.resize((alphabet_size + 1) * block_count, 0);

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
    template <bool with_overflow_check>
    void compute_step(word_type b, word_type & hp, word_type & hn, word_type & vp, word_type & vn, word_type & carry_d0, word_type & carry_hp, word_type & carry_hn)
    {
        word_type x, d0, t;

        x = b | vn;
        t = vp + (x & vp) + (with_overflow_check ? carry_d0 : 0);

        d0 = (t ^ vp) | x;
        hn = vp & d0;
        hp = vn | ~(vp | d0);

        if constexpr(with_overflow_check)
            carry_d0 = (carry_d0 != (word_type)0) ? t <= vp : t < vp;

        x = (hp << 1) | (with_overflow_check ? carry_hp : hp0);
        vn = x & d0;
        vp = (hn << 1) | ~(x | d0) | (with_overflow_check ? carry_hn : 0);

        if constexpr(with_overflow_check)
        {
            carry_hp = hp >> (word_size - 1);
            carry_hn = hn >> (word_size - 1);
        }
    }

    //!\brief Increase or decrease the score.
    void advance_score(word_type P, word_type N, word_type mask)
    {
        if ((P & mask) != (word_type)0)
            _score++;
        else if ((N & mask) != (word_type)0)
            _score--;

        if constexpr(is_semi_global)
        {
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

        last_block--;
        if (is_global && last_block == (size_t)-1)
            return false;

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

        if ((score_mask == last_score_mask) && (last_block == vp.size() - 1))
            return on_hit();
        else
        {
            next_last_active_cell();
            advance_score(vp[last_block], vn[last_block], score_mask);
        }

        return false;
    }

    //!\brief Will be called if a hit was found (e.g., score < max_errors).
    bool on_hit()
    {
        // _setFinderEnd(finder);
        //
        // if constexpr(is_global)
        //     _setFinderLength(finder, endPosition());

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
            std::size_t max_length = query.size() + max_errors + 1;
            std::size_t haystack_length = std::min(database.size(), max_length);
            database_it_end -= database.size() - haystack_length;
        }

        // distinguish between the version for needles not longer than
        // one machine word and the version for longer needles
        if (vp.size() <= 1)
            small_patterns();
        else
            large_patterns();

        if constexpr(is_global)
            _best_score = _score;
    }

public:

    /*!\brief Generic invocable interface.
     * \param[in,out] res The alignment result to fill.
     * \returns A reference to the filled alignment result.
     */
    template <typename result_type>
    result_type & operator()(result_type & res)
    {
        _compute();
        if constexpr (std::tuple_size_v<result_type> >= 2)
        {
            get<align_result_key::score>(res) = score();
        }

        if constexpr (std::tuple_size_v<result_type> >= 3)
        {
            get<align_result_key::end>(res) = end_coordinate();
        }

        [[maybe_unused]] alignment_trace_matrix matrix = trace_matrix();
        if constexpr (std::tuple_size_v<result_type> >= 4)
        {
            get<align_result_key::begin>(res) = alignment_begin_coordinate(matrix, get<align_result_key::end>(res));
        }

        if constexpr (std::tuple_size_v<result_type> >= 5)
        {
            get<align_result_key::trace>(res) = alignment_trace(database, query, matrix, get<align_result_key::end>(res));
        }
        return res;
    }

    //!\brief Return the score of the alignment.
    score_type score() const noexcept
    {
        return -_best_score;
    }

    //!\brief Return the score matrix of the alignment.
    score_matrix_type score_matrix() const noexcept
    {
        return score_matrix_type{*this};
    }

    //!\brief Return the trace matrix of the alignment.
    trace_matrix_type trace_matrix() const noexcept
    {
        return trace_matrix_type{*this};
    }

    //!\brief Return the begin position of the alignment
    alignment_coordinate begin_coordinate() const noexcept
    {
        alignment_coordinate end = end_coordinate();
        return alignment_begin_coordinate(trace_matrix(), end);
    }

    //!\brief Return the end position of the alignment
    alignment_coordinate end_coordinate() const noexcept
    {
        size_t col = database.size() - 1;
        if constexpr(is_semi_global)
            col = std::distance(begin(database), _best_score_col);

        return {col, query.size() - 1};
    }

    //!\brief Return the trace of the alignment
    auto trace() const noexcept
    {
        return alignment_trace(database, query, trace_matrix(), end_coordinate());
    }
};

template <typename database_t, typename query_t, typename align_config_t, typename traits_t>
bool pairwise_alignment_edit_distance_unbanded<database_t, query_t, align_config_t, traits_t>::small_patterns()
{
    // computing the blocks
    while (database_it != database_it_end)
    {
        word_type hn, hp, _;

        word_type b = bit_masks[to_rank((query_alphabet_type) *database_it)];
        compute_step<false>(b, hp, hn, vp[0], vn[0], _, _, _);
        advance_score(hp, hn, score_mask);

        if constexpr(use_max_errors)
            if (_score <= max_errors && on_hit())
            {
                add_state();
                ++database_it;
                return true;
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
        word_type hn, hp;
        word_type carry_d0{0}, carry_hp{hp0}, carry_hn{0};
        size_t block_offset = vp.size() * to_rank((query_alphabet_type) *database_it);

        // computing the necessary blocks, carries between blocks following one another are stored
        for (size_t current_block = 0; current_block <= last_block; current_block++)
        {
            word_type b = bit_masks[block_offset + current_block];
            compute_step<true>(b, hp, hn, vp[current_block], vn[current_block], carry_d0, carry_hp, carry_hn);
        }
        advance_score(hp, hn, score_mask);

        if constexpr(use_max_errors)
        {
            // if the active cell is the last of it's block, one additional block has to be calculated
            bool additional_block = score_mask >> (word_size - 1);
            if (last_block+1 == vp.size())
                additional_block = false;

            if (additional_block)
            {
                size_t current_block = last_block + 1;
                word_type b = bit_masks[block_offset + current_block];
                compute_step<false>(b, hp, hn, vp[current_block], vn[current_block], carry_d0, carry_hp, carry_hn);
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
template<typename database_t, typename query_t, typename config_t>
pairwise_alignment_edit_distance_unbanded(database_t && database, query_t && query, config_t config)
    -> pairwise_alignment_edit_distance_unbanded<database_t, query_t, config_t>;

template<typename database_t, typename query_t, typename config_t, typename traits_t>
pairwise_alignment_edit_distance_unbanded(database_t && database, query_t && query, config_t config, traits_t)
    -> pairwise_alignment_edit_distance_unbanded<database_t, query_t, config_t, traits_t>;
//!\}

//!\cond
template<typename database_t, typename query_t, typename align_config_t, typename traits_t>
struct alignment_score_matrix<pairwise_alignment_edit_distance_unbanded<database_t, query_t, align_config_t, traits_t>>
    : public alignment_score_matrix<std::vector<typename pairwise_alignment_edit_distance_unbanded<database_t, query_t, align_config_t, traits_t>::score_type>>
{
    using alignment_type = pairwise_alignment_edit_distance_unbanded<database_t, query_t, align_config_t, traits_t>;
    using score_type = typename alignment_type::score_type;
    using base_score_matrix_type = alignment_score_matrix<std::vector<score_type>>;
    using word_type = typename alignment_type::word_type;

    static constexpr size_t word_size = sizeof(word_type)*8;

    alignment_score_matrix(alignment_type const & alignment) :
        base_score_matrix_type
        {
            [&]{
                size_t _cols = alignment.database.size() + 1;
                size_t _rows = alignment.query.size() + 1;
                std::vector<score_type> scores{};
                scores.reserve(_cols * _rows);

                // init first row with 0, 1, 2, 3, ...
                for (size_t col=0; col < _cols; ++col)
                    scores[col] = alignment_type::is_global ? col : 0;

                auto deltas = [&](size_t col)
                {
                    return [state = alignment.states[col]](size_t row)
                    {
                        using bitset = std::bitset<word_size>;

                        size_t chunk = row / word_size;
                        size_t row_in_chunk = row % word_size;
                        word_type vp = state.vp[chunk];
                        word_type vn = state.vn[chunk];

                        int8_t p = bitset(vp)[row_in_chunk] ? 1 : 0;
                        int8_t n = bitset(vn)[row_in_chunk] ? 1 : 0;
                        return p - n;
                    };
                };

                for (size_t col=0; col < _cols; ++col)
                {
                    auto delta = deltas(col);
                    for (size_t row=1; row < _rows; ++row)
                        scores[row * _cols + col] = scores[(row-1) * _cols + col] + delta(row-1);
                }

                return scores;
            }(),
            alignment.query.size()+1,
            alignment.database.size()+1
        }
    {
    }
};

template<typename database_t, typename query_t, typename align_config_t, typename traits_t>
struct alignment_trace_matrix<pairwise_alignment_edit_distance_unbanded<database_t, query_t, align_config_t, traits_t>>
    : public alignment_trace_matrix<database_t const &, query_t const &, align_config_t, alignment_score_matrix<pairwise_alignment_edit_distance_unbanded<database_t, query_t, align_config_t, traits_t>>>
{
    using alignment_type = pairwise_alignment_edit_distance_unbanded<database_t, query_t, align_config_t, traits_t>;
    using score_matrix_type = alignment_score_matrix<alignment_type>;
    using base_trace_matrix_type = alignment_trace_matrix<database_t const &, query_t const &, align_config_t, score_matrix_type>;

    alignment_trace_matrix(alignment_type const & alignment) :
        base_trace_matrix_type{alignment.database, alignment.query, alignment.config, score_matrix_type{alignment}}
    {
    }
};

//!\endcond

} // namespace seqan3::detail
