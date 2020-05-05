// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
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

#include <seqan3/alignment/matrix/alignment_coordinate.hpp>
#include <seqan3/alignment/matrix/alignment_trace_algorithms.hpp>
#include <seqan3/alignment/matrix/edit_distance_score_matrix_full.hpp>
#include <seqan3/alignment/matrix/edit_distance_trace_matrix_full.hpp>
#include <seqan3/alignment/pairwise/alignment_result.hpp>
#include <seqan3/alignment/pairwise/edit_distance_fwd.hpp>
#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/std/ranges>

namespace seqan3::detail
{
/*!\brief Only available when default_edit_distance_trait_type::use_max_errors is true.
 * \extends default_edit_distance_trait_type
 */
template <typename derived_t, typename edit_traits>
class edit_distance_unbanded_max_errors_policy :
//!\cond
    edit_traits
//!\endcond
{
protected:
    static_assert(edit_traits::use_max_errors, "This policy assumes that edit_traits::use_max_errors is true.");

    //!\brief Befriends the derived type.
    friend derived_t;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    edit_distance_unbanded_max_errors_policy() noexcept = default;  //!< Defaulted.
    edit_distance_unbanded_max_errors_policy(edit_distance_unbanded_max_errors_policy const &) noexcept
        = default; //!< Defaulted.
    edit_distance_unbanded_max_errors_policy(edit_distance_unbanded_max_errors_policy &&) noexcept
        = default; //!< Defaulted.
    edit_distance_unbanded_max_errors_policy & operator=(edit_distance_unbanded_max_errors_policy const &) noexcept
        = default; //!< Defaulted.
    edit_distance_unbanded_max_errors_policy & operator=(edit_distance_unbanded_max_errors_policy &&) noexcept
        = default; //!< Defaulted.
    ~edit_distance_unbanded_max_errors_policy() noexcept = default; //!< Defaulted.
    //!\}

    using typename edit_traits::word_type;
    using typename edit_traits::score_type;

    /*!\name Max Error Policy: Protected Attributes
     * \copydoc edit_distance_unbanded_max_errors_policy
     * \{
     */
    //!\brief Which score value is considered as a hit?
    score_type max_errors{255};
    //!\brief The block containing the last active cell.
    size_t last_block{0u};
    //!\brief A mask with a bit set on the position of the last row.
    word_type last_score_mask{};
    //!\}

    /*!\name Max Error Policy: Protected Member Functions
     * \copydoc edit_distance_unbanded_max_errors_policy
     * \{
     */
    //!\brief Initialises max_errors policy.
    void max_errors_init(size_t block_count) noexcept
    {
        derived_t * self = static_cast<derived_t *>(this);

        max_errors = get<align_cfg::max_error>(self->config).value;
        assert(max_errors >= score_type{0});

        if (std::ranges::empty(self->query)) // [[unlikely]]
        {
            last_block = 0u;
            self->score_mask = 0u;
            last_score_mask = self->score_mask;
            return;
        }

        last_block = block_count - 1u;
        last_score_mask = self->score_mask;

        // local_max_errors either stores the maximal number of _score (me.max_errors) or the needle size minus one. It
        // is used for the mask computation and setting the initial score (the minus one is there because of the Ukkonen
        // trick).
        size_t const local_max_errors = std::min<size_t>(max_errors, std::ranges::size(self->query) - 1u);
        self->score_mask = word_type{1u} << (local_max_errors % self->word_size);
        last_block = std::min(local_max_errors / self->word_size, last_block);
        self->_score = local_max_errors + 1u;
    }

    //!\brief Returns true if the current active cell is within the last row.
    bool is_last_active_cell_within_last_row() const noexcept
    {
        derived_t const * self = static_cast<derived_t const *>(this);
        return (self->score_mask == this->last_score_mask) && (this->last_block == self->vp.size() - 1u);
    }

    //!\brief Decrement the last active cell position.
    bool prev_last_active_cell() noexcept
    {
        derived_t * self = static_cast<derived_t *>(this);
        self->score_mask >>= 1u;
        if (self->score_mask != 0u)
            return true;

        if constexpr (edit_traits::is_global)
        {
            if (last_block == 0u) // [[unlikely]]
                return false;
        }

        last_block--;

        self->score_mask = word_type{1u} << (edit_traits::word_size - 1u);
        return true;
    }

    //!\brief Increment the last active cell position.
    void next_last_active_cell() noexcept
    {
        derived_t * self = static_cast<derived_t *>(this);
        self->score_mask <<= 1u;
        if (self->score_mask)
            return;

        self->score_mask = 1u;
        last_block++;
    }

    /*!\brief Use the ukkonen trick and update the last active cell.
     * \returns `true` if computation should be aborted, `false` if computation should continue.
     */
    bool update_last_active_cell() noexcept
    {
        derived_t * self = static_cast<derived_t *>(this);
        // update the last active cell
        while (!(self->_score <= max_errors))
        {
            self->advance_score(self->vn[last_block], self->vp[last_block], self->score_mask);
            if (!prev_last_active_cell())
            {
                // prev_last_active_cell = false can only happen for global alignments
                assert(edit_traits::is_global);
                // we abort here if we don't need to compute a matrix, because the continued
                // computation can't produce an alignment.
                return !edit_traits::compute_matrix;
            }
        }

        if (is_last_active_cell_within_last_row())
        {
            assert(self->_score <= max_errors);

            if constexpr(edit_traits::is_semi_global)
                self->update_best_score();

            return self->on_hit();
        }
        else
        {
            next_last_active_cell();
            self->advance_score(self->vp[last_block], self->vn[last_block], self->score_mask);
        }

        return false;
    }
    //!\}

    //!\copydoc edit_distance_score_matrix_full::max_rows
    static size_t max_rows(word_type const score_mask, unsigned const last_block,
                           score_type const score, score_type const max_errors) noexcept
    {
        using score_matrix_type = typename edit_traits::score_matrix_type;
        return score_matrix_type::max_rows(score_mask,
                                           last_block,
                                           score,
                                           max_errors);
    }
};

/*!\brief Only available when default_edit_distance_trait_type::is_global is true.
 * \extends default_edit_distance_trait_type
 */
template <typename derived_t, typename edit_traits>
class edit_distance_unbanded_global_policy :
//!\cond
    edit_traits
//!\endcond
{
protected:
    static_assert(edit_traits::is_global || edit_traits::is_semi_global,
                  "This policy assumes that edit_traits::is_global or edit_traits::is_semi_global is true.");

    //!\brief Befriends the derived type.
    friend derived_t;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    edit_distance_unbanded_global_policy() noexcept = default;  //!< Defaulted.
    edit_distance_unbanded_global_policy(edit_distance_unbanded_global_policy const &) noexcept
        = default; //!< Defaulted.
    edit_distance_unbanded_global_policy(edit_distance_unbanded_global_policy &&) noexcept
        = default; //!< Defaulted.
    edit_distance_unbanded_global_policy & operator=(edit_distance_unbanded_global_policy const &) noexcept
        = default; //!< Defaulted.
    edit_distance_unbanded_global_policy & operator=(edit_distance_unbanded_global_policy &&) noexcept
        = default; //!< Defaulted.
    ~edit_distance_unbanded_global_policy() noexcept = default; //!< Defaulted.
    //!\}

    //!\copydoc default_edit_distance_trait_type::score_type
    using typename edit_traits::score_type;

    /*!\name Global Policy: Protected Attributes
     * \copydoc edit_distance_unbanded_global_policy
     * \{
     */
    /*!\brief The best score of the alignment in the last row
     * (if is_semi_global = true) or the last entry in the
     * score matrix (if is_global = true).
     */
    score_type _best_score{};
    //!\}

    /*!\name Global Policy: Protected Member Functions
     * \copydoc edit_distance_unbanded_global_policy
     * \{
     */
    //!\brief Initialises global policy.
    void score_init() noexcept
    {
        derived_t const * self = static_cast<derived_t const *>(this);
        _best_score = self->_score;
    }

    //!\brief Returns true if the computation produced a valid alignment.
    bool is_valid() const noexcept
    {
        [[maybe_unused]] derived_t const * self = static_cast<derived_t const *>(this);
        // This condition uses the observation that after each computation of a column, _score has either the initial
        // value of the first row (i.e. the entire column consist of INF's), has the value _score = max_errors + 1
        // (there exists a cell within the column that has value <= max_errors, but is not on the last row) or _score <=
        // max_errors (the score of the last active cell is <= max_errors)
        if constexpr(edit_traits::use_max_errors)
            return _best_score <= self->max_errors;

        // When not using max_errors there is always a valid alignment, because the last row will always be updated and
        // with it the score.
        return true;
    }

    //!\brief Returns an invalid_coordinate for this alignment.
    alignment_coordinate invalid_coordinate() const noexcept
    {
        derived_t const * self = static_cast<derived_t const *>(this);
        return {column_index_type{std::ranges::size(self->database)}, row_index_type{std::ranges::size(self->query)}};
    }

    //!\brief Update the current best known score if the current score is better.
    void update_best_score() noexcept
    {
        derived_t const * self = static_cast<derived_t const *>(this);
        _best_score = self->_score;
    }

    //!\brief Returns the first component of the #back_coordinate.
    size_t back_coordinate_first() const noexcept
    {
        derived_t const * self = static_cast<derived_t const *>(this);
        return std::ranges::size(self->database);
    }
    //!\}

public:
    /*!\name Global Policy: Public Member Functions
     * \copydoc edit_distance_unbanded_global_policy
     * \{
     */
    //!\brief Return the score of the alignment.
    //!       Only available if default_edit_distance_trait_type::compute_score is true.
    std::optional<score_type> score() const noexcept
    {
        derived_t const * self = static_cast<derived_t const *>(this);
        static_assert(edit_traits::compute_score, "score() can only be computed if you specify the result type within "
                                                  "your alignment config.");
        if (!self->is_valid())
            return std::nullopt;

        return -_best_score;
    }

    //!\brief Return the end position of the alignment
    //!       Only available if default_edit_distance_trait_type::compute_back_coordinate is true.
    alignment_coordinate back_coordinate() const noexcept
    {
        derived_t const * self = static_cast<derived_t const *>(this);
        static_assert(edit_traits::compute_back_coordinate, "back_coordinate() can only be computed if you specify the"
                                                            "result type within your alignment config.");
        if (!self->is_valid())
            return self->invalid_coordinate();

        column_index_type const first{self->back_coordinate_first()};
        row_index_type const second{std::ranges::size(self->query)};
        return {first, second};
    }
    //!\}
};

//!\brief Only available when default_edit_distance_trait_type::is_semi_global is true.
template <typename derived_t, typename edit_traits>
class edit_distance_unbanded_semi_global_policy :
    public edit_distance_unbanded_global_policy<derived_t, edit_traits>
{
protected:
    static_assert(edit_traits::is_semi_global, "This policy assumes that edit_traits::is_semi_global is true.");

    //!\brief Befriends the derived type.
    friend derived_t;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    edit_distance_unbanded_semi_global_policy() noexcept = default;  //!< Defaulted.
    edit_distance_unbanded_semi_global_policy(edit_distance_unbanded_semi_global_policy const &) noexcept
        = default; //!< Defaulted.
    edit_distance_unbanded_semi_global_policy(edit_distance_unbanded_semi_global_policy &&) noexcept
        = default; //!< Defaulted.
    edit_distance_unbanded_semi_global_policy & operator=(edit_distance_unbanded_semi_global_policy const &) noexcept
        = default; //!< Defaulted.
    edit_distance_unbanded_semi_global_policy & operator=(edit_distance_unbanded_semi_global_policy &&) noexcept
        = default; //!< Defaulted.
    ~edit_distance_unbanded_semi_global_policy() noexcept = default; //!< Defaulted.
    //!\}

    //!\brief The base policy of this policy.
    using base_t = edit_distance_unbanded_global_policy<derived_t, edit_traits>;
    //!\cond
    using database_iterator = typename edit_traits::database_iterator;
    using base_t::_best_score;
    //!\endcond
    /*!\name Semi-Global Policy: Protected Attributes
     * \copydoc edit_distance_unbanded_semi_global_policy
     * \{
     */
    /*!\brief In which column the best score of the alignment is located. Will only be tracked if is_semi_global is
     *        true.
     *
     * \details
     *
     * If is_global is true this is always at the last entry in the score matrix, i.e. at position (`|query|`,
     * `|database|`).
     */
    database_iterator _best_score_col{};
    //!\}

    /*!\name Semi-Global Policy: Protected Member Functions
     * \copydoc edit_distance_unbanded_semi_global_policy
     * \{
     */
    //!\brief Initialises semi-global policy.
    void score_init() noexcept
    {
        derived_t const * self = static_cast<derived_t const *>(this);
        base_t::score_init();
        _best_score_col = self->database_it_end;
    }

    //!\copydoc edit_distance_unbanded_global_policy::update_best_score
    void update_best_score() noexcept
    {
        derived_t const * self = static_cast<derived_t const *>(this);
        // we have to make sure that update_best_score is only called after a score update within the last row.
        if constexpr(edit_traits::use_max_errors)
        {
            assert(std::ranges::empty(self->query) || self->is_last_active_cell_within_last_row());
        }

        _best_score_col = (self->_score <= _best_score) ? self->database_it : _best_score_col;
        _best_score     = (self->_score <= _best_score) ? self->_score : _best_score;
    }

    //!\copydoc edit_distance_unbanded_global_policy::back_coordinate_first
    size_t back_coordinate_first() const noexcept
    {
        derived_t const * self = static_cast<derived_t const *>(this);
        // offset == 0u is a special case if database sequence is empty, because in this case the best column is zero.
        size_t offset = std::ranges::empty(self->database) ? 0u : 1u;
        return std::ranges::distance(std::ranges::begin(self->database), _best_score_col) + offset;
    }
    //!\}
};

/*!\brief Only available when default_edit_distance_trait_type::compute_score_matrix is true.
 * \extends default_edit_distance_trait_type
 */
template <typename derived_t, typename edit_traits>
class edit_distance_unbanded_score_matrix_policy :
//!\cond
    edit_traits
//!\endcond
{
protected:
    static_assert(edit_traits::compute_score_matrix,
                  "This policy assumes that edit_traits::compute_score_matrix is true.");

    //!\brief Befriends the derived type.
    friend derived_t;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    edit_distance_unbanded_score_matrix_policy() noexcept = default;  //!< Defaulted.
    edit_distance_unbanded_score_matrix_policy(edit_distance_unbanded_score_matrix_policy const &) noexcept
        = default; //!< Defaulted.
    edit_distance_unbanded_score_matrix_policy(edit_distance_unbanded_score_matrix_policy &&) noexcept
        = default; //!< Defaulted.
    edit_distance_unbanded_score_matrix_policy & operator=(edit_distance_unbanded_score_matrix_policy const &) noexcept
        = default; //!< Defaulted.
    edit_distance_unbanded_score_matrix_policy & operator=(edit_distance_unbanded_score_matrix_policy &&) noexcept
        = default; //!< Defaulted.
    ~edit_distance_unbanded_score_matrix_policy() noexcept = default; //!< Defaulted.
    //!\}

    using typename edit_traits::score_matrix_type;

    /*!\name Score Matrix Policy: Protected Attributes
     * \copydoc edit_distance_unbanded_score_matrix_policy
     * \{
     */
    //!\brief The score matrix of the edit distance alignment.
    score_matrix_type _score_matrix{};
    //!\}

    /*!\name Score Matrix Policy: Protected Member Functions
     * \copydoc edit_distance_unbanded_score_matrix_policy
     * \{
     */
    /*!\brief Initialises score-matrix policy.
     *
     * \details
     *
     * ### Exception
     *
     * Strong exception guarantee.
     */
    void score_matrix_init()
    {
        derived_t const * self = static_cast<derived_t const *>(this);

        _score_matrix = score_matrix_type{std::ranges::size(self->query) + 1u};
        _score_matrix.reserve(std::ranges::size(self->database) + 1u);
    }
    //!\}

public:
    /*!\name Score Matrix Policy: Public Member Functions
     * \copydoc edit_distance_unbanded_score_matrix_policy
     * \{
     */
    /*!\brief Return the score matrix of the alignment.
     * Only available if default_edit_distance_trait_type::compute_score_matrix is true.
     */
    score_matrix_type const & score_matrix() const noexcept
    {
        static_assert(edit_traits::compute_score_matrix, "score_matrix() can only be computed if you specify the "
                                                         "result type within your alignment config.");
        return _score_matrix;
    }
    //!\}
};

/*!\brief Only available when default_edit_distance_trait_type::compute_trace_matrix is true.
 * \extends default_edit_distance_trait_type
 */
template <typename derived_t, typename edit_traits>
class edit_distance_unbanded_trace_matrix_policy :
//!\cond
    edit_traits
//!\endcond
{
protected:
    static_assert(edit_traits::compute_trace_matrix,
                  "This policy assumes that edit_traits::compute_trace_matrix is true.");

    //!\brief Befriends the derived type.
    friend derived_t;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    edit_distance_unbanded_trace_matrix_policy() noexcept = default;  //!< Defaulted.
    edit_distance_unbanded_trace_matrix_policy(edit_distance_unbanded_trace_matrix_policy const &) noexcept
        = default; //!< Defaulted.
    edit_distance_unbanded_trace_matrix_policy(edit_distance_unbanded_trace_matrix_policy &&) noexcept
        = default; //!< Defaulted.
    edit_distance_unbanded_trace_matrix_policy & operator=(edit_distance_unbanded_trace_matrix_policy const &) noexcept
        = default; //!< Defaulted.
    edit_distance_unbanded_trace_matrix_policy & operator=(edit_distance_unbanded_trace_matrix_policy &&) noexcept
        = default; //!< Defaulted.
    ~edit_distance_unbanded_trace_matrix_policy() noexcept = default; //!< Defaulted.
    //!\}

    using typename edit_traits::word_type;
    using typename edit_traits::trace_matrix_type;
    using typename edit_traits::alignment_result_type;

    /*!\name Trace matrix Policy: Protected Attributes
     * \copydoc edit_distance_unbanded_trace_matrix_policy
     * \{
     */
    //!\copydoc edit_distance_unbanded::compute_state::hp
    std::vector<word_type> hp{};
    //!\copydoc edit_distance_unbanded::compute_state_trace_matrix::db
    std::vector<word_type> db{};

    //!\brief The trace matrix of the edit distance alignment.
    trace_matrix_type _trace_matrix{};
    //!\}

    /*!\name Trace matrix Policy: Protected Member Functions
     * \copydoc edit_distance_unbanded_trace_matrix_policy
     * \{
     */
    /*!\brief Initialises trace-matrix policy.
     *
     * \details
     *
     * ### Exception
     *
     * Strong exception guarantee.
     */
    void trace_matrix_init(size_t block_count)
    {
        derived_t const * self = static_cast<derived_t const *>(this);

        _trace_matrix = trace_matrix_type{std::ranges::size(self->query) + 1u};
        _trace_matrix.reserve(std::ranges::size(self->database) + 1u);

        hp.resize(block_count, 0u);
        db.resize(block_count, 0u);
    }
    //!\}

public:
    /*!\name Trace matrix Policy: Public Member Functions
     * \copydoc edit_distance_unbanded_trace_matrix_policy
     * \{
     */
    //!\brief Return the trace matrix of the alignment.
    trace_matrix_type const & trace_matrix() const noexcept
    {
        static_assert(edit_traits::compute_trace_matrix, "trace_matrix() can only be computed if you specify the "
                                                         "result type within your alignment config.");
        return _trace_matrix;
    }

    //!\brief Return the begin position of the alignment.
    //!       Only available if default_edit_distance_trait_type::compute_front_coordinate is true.
    alignment_coordinate front_coordinate() const noexcept
    {
        derived_t const * self = static_cast<derived_t const *>(this);
        static_assert(edit_traits::compute_front_coordinate, "front_coordinate() can only be computed if you specify "
                                                             "the result type within your alignment config.");
        if (!self->is_valid())
            return self->invalid_coordinate();

        alignment_coordinate const back = self->back_coordinate();
        return alignment_front_coordinate(trace_matrix(), back);
    }

    //!\brief Return the alignment, i.e. the actual base pair matching.
    //!       Only available if default_edit_distance_trait_type::compute_sequence_alignment is true.
    auto alignment() const noexcept
    {
        using alignment_t = remove_cvref_t<decltype(std::declval<alignment_result_type &>().alignment())>;

        derived_t const * self = static_cast<derived_t const *>(this);
        static_assert(edit_traits::compute_sequence_alignment, "alignment() can only be computed if you specify the "
                                                               "result type within your alignment config.");

        if (!self->is_valid())
            return alignment_t{};

        return alignment_trace<alignment_t>(self->database,
                                            self->query,
                                            trace_matrix(),
                                            self->back_coordinate(),
                                            front_coordinate());
    }
    //!\}
};

/*!\brief The same as `value_t &` but it is default constructible and is re-assignable.
 * \tparam value_t The value type of the reference.
 */
template <typename value_t>
class proxy_reference
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    proxy_reference() noexcept = default;                                    //!< Defaulted.
    proxy_reference(proxy_reference const &) noexcept = default;             //!< Defaulted.
    proxy_reference(proxy_reference &&) noexcept = default;                  //!< Defaulted.
    proxy_reference & operator=(proxy_reference const &) noexcept = default; //!< Defaulted.
    proxy_reference & operator=(proxy_reference &&) noexcept = default;      //!< Defaulted.
    ~proxy_reference() = default;                                            //!< Defaulted.

    //!\brief Use the lvalue `t` as the stored reference.
    proxy_reference(value_t & t) noexcept
        : ptr(std::addressof(t))
    {}

    proxy_reference(value_t &&) = delete; //!< Deleted.

    //!\brief Assign a value to the stored reference.
    template <typename other_value_t>
    //!\cond
        requires std::convertible_to<other_value_t, value_t>
    //!\endcond
    proxy_reference & operator=(other_value_t && u) noexcept
    {
        get() = std::forward<other_value_t>(u);
        return *this;
    }
    //!\}

    //!\brief Get the stored reference.
    value_t & get() const noexcept
    {
        assert(ptr != nullptr);
        return *ptr;
    }

    //!\brief Get the stored reference.
    operator value_t & () const noexcept
    {
        return get();
    }

private:
    //!\brief The stored reference.
    value_t * ptr{nullptr};
};

/*!\brief This calculates an alignment using the edit distance and without a band.
 * \ingroup pairwise_alignment
 * \tparam database_t     \copydoc default_edit_distance_trait_type::database_type
 * \tparam query_t        \copydoc default_edit_distance_trait_type::query_type
 * \tparam align_config_t The configuration type; must be of type seqan3::configuration.
 * \extends edit_distance_unbanded_global_policy
 * \extends edit_distance_unbanded_semi_global_policy
 * \extends edit_distance_unbanded_score_matrix_policy
 * \extends edit_distance_unbanded_trace_matrix_policy
 * \extends edit_distance_unbanded_max_errors_policy
 */
template <std::ranges::viewable_range database_t,
          std::ranges::viewable_range query_t,
          typename align_config_t,
          typename edit_traits>
class edit_distance_unbanded :
//!\cond
// Hide this section in doxygen, because it messes up the inheritance.
    public edit_distance_base<
        edit_traits::use_max_errors,
        edit_distance_unbanded_max_errors_policy,
        edit_traits,
        edit_distance_unbanded<database_t, query_t, align_config_t, edit_traits>>,
    public edit_distance_base<
        edit_traits::is_global,
        edit_distance_unbanded_global_policy,
        edit_traits,
        edit_distance_unbanded<database_t, query_t, align_config_t, edit_traits>>,
    public edit_distance_base<
        edit_traits::is_semi_global,
        edit_distance_unbanded_semi_global_policy,
        edit_traits,
        edit_distance_unbanded<database_t, query_t, align_config_t, edit_traits>>,
    public edit_distance_base<
        edit_traits::compute_score_matrix,
        edit_distance_unbanded_score_matrix_policy,
        edit_traits,
        edit_distance_unbanded<database_t, query_t, align_config_t, edit_traits>>,
    public edit_distance_base<
        edit_traits::compute_trace_matrix,
        edit_distance_unbanded_trace_matrix_policy,
        edit_traits,
        edit_distance_unbanded<database_t, query_t, align_config_t, edit_traits>>
//!\endcond
{
public:
    using typename edit_traits::word_type;
    using typename edit_traits::score_type;
    using typename edit_traits::database_type;
    using typename edit_traits::query_type;
    using typename edit_traits::align_config_type;
    using edit_traits::word_size;

private:
    //!\brief Allows seqan3::detail::edit_distance_unbanded_max_errors_policy to access this class.
    template <typename other_derived_t, typename other_edit_traits>
    friend class edit_distance_unbanded_max_errors_policy;
    //!\brief Allows seqan3::detail::edit_distance_unbanded_global_policy to access this class.
    template <typename other_derived_t, typename other_edit_traits>
    friend class edit_distance_unbanded_global_policy;
    //!\brief Allows seqan3::detail::edit_distance_unbanded_semi_global_policy to access this class.
    template <typename other_derived_t, typename other_edit_traits>
    friend class edit_distance_unbanded_semi_global_policy;
    //!\brief Allows seqan3::detail::edit_distance_unbanded_score_matrix_policy to access this class.
    template <typename other_derived_t, typename other_edit_traits>
    friend class edit_distance_unbanded_score_matrix_policy;
    //!\brief Allows seqan3::detail::edit_distance_unbanded_trace_matrix_policy to access this class.
    template <typename other_derived_t, typename other_edit_traits>
    friend class edit_distance_unbanded_trace_matrix_policy;

    using typename edit_traits::database_iterator;
    using typename edit_traits::query_alphabet_type;
    using typename edit_traits::alignment_result_type;
    using edit_traits::use_max_errors;
    using edit_traits::is_semi_global;
    using edit_traits::is_global;
    using edit_traits::compute_score;
    using edit_traits::compute_back_coordinate;
    using edit_traits::compute_front_coordinate;
    using edit_traits::compute_sequence_alignment;
    using edit_traits::compute_score_matrix;
    using edit_traits::compute_trace_matrix;
    using edit_traits::compute_matrix;
    using typename edit_traits::score_matrix_type;
    using typename edit_traits::trace_matrix_type;

    //!\brief The horizontal/database sequence.
    database_t database;
    //!\brief The vertical/query sequence.
    query_t query;
    //!\brief The configuration.
    align_config_t config;

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
    /*!\brief The mask with a bit set at the position where the score change.
     *
     * \details
     *
     * If #use_max_errors is true this corresponds to the last active cell.
     */
    word_type score_mask{0u};
    //!\copydoc edit_distance_unbanded::compute_state::vp
    std::vector<word_type> vp{};
    //!\copydoc edit_distance_unbanded::compute_state::vn
    std::vector<word_type> vn{};
    /*!\brief The machine words which translate a letter of the query into a bit mask.
     *
     * \details
     *
     * Each bit position which is true (= 1) corresponds to a match of a letter in the query at this position.
     */
    std::vector<word_type> bit_masks{};

    //!\brief The current position in the database.
    database_iterator database_it{};
    //!\brief The end position of the database.
    database_iterator database_it_end{};

    //!\brief The internal state needed to compute the trace matrix.
    struct compute_state_trace_matrix
    {
        //!\brief The machine word which stores if trace_directions::diagonal is true.
        proxy_reference<word_type> db{};
    };

    //!\brief The internal state needed to compute the alignment.
    struct compute_state : enable_state_t<compute_trace_matrix, compute_state_trace_matrix>
    {
        //!\brief The type of hp.
        using hp_type = std::conditional_t<compute_trace_matrix, proxy_reference<word_type>, word_type>;

        //!\brief The machine word which stores wether the current character matches.
        word_type b{};
        //!\brief The machine word which stores the diagonal differences.
        word_type d0{};
        //!\brief The machine word which stores the positive horizontal differences.
        hp_type hp{};
        //!\brief The machine word which stores the negative horizontal differences.
        word_type hn{};
        //!\brief The machine word which stores the positive vertical differences.
        proxy_reference<word_type> vp{};
        //!\brief The machine word which stores the negative vertical differences.
        proxy_reference<word_type> vn{};
        //!\brief The carry-bit of d0.
        word_type carry_d0{};
        //!\brief The carry-bit of hp.
        word_type carry_hp{hp0};
        //!\brief The carry-bit of hn.
        word_type carry_hn{};
    };

    //!\brief Add a computation step
    void add_state()
    {
        if constexpr(!use_max_errors && compute_score_matrix)
            this->_score_matrix.add_column(vp, vn);

        if constexpr(!use_max_errors && compute_trace_matrix)
            this->_trace_matrix.add_column(this->hp, this->db, vp);

        if constexpr(use_max_errors && compute_matrix)
        {
            size_t max_rows = this->max_rows(score_mask, this->last_block, _score, this->max_errors);
            if constexpr(compute_score_matrix)
                this->_score_matrix.add_column(vp, vn, max_rows);

            if constexpr(compute_trace_matrix)
                this->_trace_matrix.add_column(this->hp, this->db, vp, max_rows);
        }
    }

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
     //!\brief The class template parameter may resolve to an lvalue reference which prohibits default constructibility.
     edit_distance_unbanded() = delete;                                            //!< Defaulted.
     edit_distance_unbanded(edit_distance_unbanded const &) = default;             //!< Defaulted.
     edit_distance_unbanded(edit_distance_unbanded &&) = default;                  //!< Defaulted.
     edit_distance_unbanded & operator=(edit_distance_unbanded const &) = default; //!< Defaulted.
     edit_distance_unbanded & operator=(edit_distance_unbanded &&) = default;      //!< Defaulted.
     ~edit_distance_unbanded() = default;                                          //!< Defaulted.

    /*!\brief Constructor
     * \param[in] _database \copydoc database
     * \param[in] _query    \copydoc query
     * \param[in] _config   \copydoc config
     * \param[in] _traits   The traits object. Only the type information will be used.
     */
    edit_distance_unbanded(database_t _database,
                           query_t _query,
                           align_config_t _config,
                           edit_traits const & SEQAN3_DOXYGEN_ONLY(_traits)) :
        database{std::forward<database_t>(_database)},
        query{std::forward<query_t>(_query)},
        config{std::forward<align_config_t>(_config)},
        _score{static_cast<score_type>(std::ranges::size(query))},
        database_it{ranges::begin(database)},
        database_it_end{ranges::end(database)}
    {
        static constexpr size_t alphabet_size_ = alphabet_size<query_alphabet_type>;

        size_t const block_count = (std::ranges::size(query) - 1u + word_size) / word_size;
        score_mask = word_type{1u} << ((std::ranges::size(query) - 1u + word_size) % word_size);

        this->score_init();
        if constexpr(use_max_errors)
            this->max_errors_init(block_count);

        if constexpr(compute_score_matrix)
            this->score_matrix_init();

        if constexpr(compute_trace_matrix)
            this->trace_matrix_init(block_count);

        vp.resize(block_count, vp0);
        vn.resize(block_count, vn0);
        bit_masks.resize((alphabet_size_ + 1u) * block_count, 0u);

        // encoding the letters as bit-vectors
        for (size_t j = 0u; j < std::ranges::size(query); j++)
        {
            size_t const i = block_count * seqan3::to_rank(query[j]) + j / word_size;
            bit_masks[i] |= word_type{1u} << (j % word_size);
        }

        add_state();
    }
    //!\}

private:
    //!\brief A single compute step in the current column.
    template <bool with_carry>
    static void compute_step(compute_state & state) noexcept
    {
        word_type x, t;
        assert(state.carry_d0 <= 1u);
        assert(state.carry_hp <= 1u);
        assert(state.carry_hn <= 1u);

        x = state.b | state.vn;
        t = state.vp + (x & state.vp) + state.carry_d0;

        state.d0 = (t ^ state.vp) | x;
        state.hn = state.vp & state.d0;
        state.hp = state.vn | ~(state.vp | state.d0);

        if constexpr(with_carry)
            state.carry_d0 = (state.carry_d0 != 0u) ? t <= state.vp : t < state.vp;

        x = (state.hp << 1u) | state.carry_hp;
        state.vn = x & state.d0;
        state.vp = (state.hn << 1u) | ~(x | state.d0) | state.carry_hn;

        if constexpr(with_carry)
        {
            state.carry_hp = state.hp >> (word_size - 1u);
            state.carry_hn = state.hn >> (word_size - 1u);
        }
    }

    //!\brief A single compute step in the current column at a given position.
    template <bool with_carry>
    void compute_kernel(compute_state & state, size_t const block_offset, size_t const current_block) noexcept
    {
        state.vp = proxy_reference<word_type>{this->vp[current_block]};
        state.vn = proxy_reference<word_type>{this->vn[current_block]};
        if constexpr(compute_trace_matrix)
        {
            state.hp = proxy_reference<word_type>{this->hp[current_block]};
            state.db = proxy_reference<word_type>{this->db[current_block]};
        }
        state.b = bit_masks[block_offset + current_block];

        compute_step<with_carry>(state);
        if constexpr(compute_trace_matrix)
            state.db = ~(state.b ^ state.d0);
    }

    //!\brief Increase or decrease the score.
    void advance_score(word_type P, word_type N, word_type mask) noexcept
    {
        if ((P & mask) != word_type{0u})
            _score++;
        else if ((N & mask) != word_type{0u})
            _score--;
    }

    //!\brief Will be called if a hit was found (e.g., score <= max_errors).
    bool on_hit() noexcept
    {
        // TODO: call external on_hit functor
        return false;
    }

    //!\brief Pattern is small enough that it fits into one machine word. Use faster computation with less overhead.
    inline bool small_patterns();

    //!\brief Pattern is larger than one machine word. Use overflow aware computation.
    inline bool large_patterns();

    //!\brief Special case if query sequence is empty.
    inline void compute_empty_query_sequence()
    {
        assert(std::ranges::empty(query));

        bool abort_computation = false;

        for (; database_it != database_it_end; ++database_it)
        {
            if constexpr(is_global)
                ++_score;
            else // is_semi_global
                this->update_best_score();

            // call on_hit
            if constexpr(use_max_errors)
                abort_computation = on_hit();

            this->add_state();
            if (abort_computation)
                break;
        }
    }

    //!\brief Compute the alignment.
    void compute()
    {
        // limit search width for prefix search (if no matrix needs to be computed)
        if constexpr(use_max_errors && is_global && !compute_matrix)
        {
            // Note: For global alignments we know that the database can only be max_length long to have a score less
            // than or equal max_errors in the last cell.
            //
            // The following matrix shows a minimal value for each entry (because a diagonal must be +0 or +1, each
            // diagonal is at least the value of the initial value in the first row)
            // 0 1 2 3 4 5 6...
            // 1 0 1 2 3 4 5...
            // ...
            // m ... 3 2 1 0 1 2 3 4 5
            // Thus, after |query| + max_errors entries the score will always be higher than max_errors.
            size_t const max_length = std::ranges::size(query) + this->max_errors + 1u;
            size_t const haystack_length = std::min(std::ranges::size(database), max_length);
            database_it_end -= std::ranges::size(database) - haystack_length;
        }

        // distinguish between the version for needles not longer than
        // one machine word and the version for longer needles
        // A special cases is if the second sequence is empty (vp.size() == 0u).
        if (vp.size() == 0u) // [[unlikely]]
            compute_empty_query_sequence();
        else if (vp.size() == 1u)
            small_patterns();
        else
            large_patterns();

        if constexpr(is_global)
            this->update_best_score();
    }

public:
    /*!\brief Generic invocable interface.
     * \param[in] idx The index of the currently processed sequence pair.
     * \param[in] callback The callback function to be invoked with the alignment result.
     */
    template <typename callback_t>
    void operator()(size_t const idx, callback_t && callback)
    {
        using result_value_type = typename alignment_result_value_type_accessor<alignment_result_type>::type;

        compute();
        result_value_type res_vt{};
        res_vt.id = idx;
        if constexpr (compute_score)
        {
            res_vt.score = this->score().value_or(matrix_inf<score_type>);
        }

        if constexpr (compute_back_coordinate)
        {
            res_vt.back_coordinate = this->back_coordinate();
        }

        if constexpr (compute_front_coordinate)
        {
            if (this->is_valid())
                res_vt.front_coordinate = alignment_front_coordinate(this->trace_matrix(), res_vt.back_coordinate);
            else
                res_vt.front_coordinate = this->invalid_coordinate();
        }

        if constexpr (compute_sequence_alignment)
        {
            if (this->is_valid())
            {
                using alignment_t = decltype(res_vt.alignment);
                res_vt.alignment = alignment_trace<alignment_t>(database,
                                                                query,
                                                                this->trace_matrix(),
                                                                res_vt.back_coordinate,
                                                                res_vt.front_coordinate);
            }
        }
        callback(alignment_result_type{std::move(res_vt)});
    }
};

template <typename database_t, typename query_t, typename align_config_t, typename traits_t>
bool edit_distance_unbanded<database_t, query_t, align_config_t, traits_t>::small_patterns()
{
    bool abort_computation = false;

    // computing the blocks
    while (database_it != database_it_end)
    {
        compute_state state{};
        size_t const block_offset = seqan3::to_rank((query_alphabet_type) *database_it);

        compute_kernel<false>(state, block_offset, 0u);
        advance_score(state.hp, state.hn, score_mask);

        // semi-global without max_errors guarantees that the score stays within the last row
        if constexpr(is_semi_global && !use_max_errors)
            this->update_best_score();

        // updating the last active cell
        if constexpr(use_max_errors)
            abort_computation = this->update_last_active_cell();

        add_state();
        ++database_it;
        if (abort_computation)
            return true;
    }

    return false;
}

template <typename database_t, typename query_t, typename align_config_t, typename traits_t>
bool edit_distance_unbanded<database_t, query_t, align_config_t, traits_t>::large_patterns()
{
    bool abort_computation = false;

    while (database_it != database_it_end)
    {
        compute_state state{};
        size_t const block_offset = vp.size() * seqan3::to_rank((query_alphabet_type) *database_it);

        size_t block_count = vp.size();
        if constexpr(use_max_errors)
            block_count = this->last_block + 1;

        // compute each block in the current column; carries between blocks will be propagated.
        for (size_t current_block = 0u; current_block < block_count; current_block++)
            compute_kernel<true>(state, block_offset, current_block);

        advance_score(state.hp, state.hn, score_mask);

        // semi-global without max_errors guarantees that the score stays within the last row
        if constexpr(is_semi_global && !use_max_errors)
            this->update_best_score();

        if constexpr(use_max_errors)
        {
            // if the last active cell reached the end within the current block we have to compute the next block.
            bool additional_block = score_mask >> (word_size - 1u);
            bool reached_last_block = this->last_block + 1u == vp.size();
            // If there is no next block we skip the computation.
            if (reached_last_block)
                additional_block = false;

            if (additional_block)
            {
                size_t const current_block = this->last_block + 1u;
                // this might not be necessary, but carry_d0 = 1u might have an influence on the result of vn and vp.
                vp[current_block] = vp0;
                vn[current_block] = vn0;
                compute_kernel<false>(state, block_offset, current_block);
            }

            // updating the last active cell
            abort_computation = this->update_last_active_cell();
        }

        add_state();
        ++database_it;

        if (abort_computation)
            return true;
    }

    return false;
}

/*!\name Type deduction guides
 * \relates seqan3::detail::edit_distance_unbanded
 * \{
 */

//!\brief Deduce the type from the provided arguments.
template <typename database_t, typename query_t, typename config_t, typename traits_t>
edit_distance_unbanded(database_t && database, query_t && query, config_t config, traits_t)
    -> edit_distance_unbanded<database_t, query_t, config_t, traits_t>;
//!\}

} // namespace seqan3::detail
