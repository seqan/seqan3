// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Provides the declaration of seqan3::detail::debug_matrix.
 */

#pragma once

#include <iomanip>

#include <seqan3/alignment/matrix/matrix_concept.hpp>
#include <seqan3/alignment/matrix/row_wise_matrix.hpp>
#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/core/detail/debug_stream_alphabet.hpp>
#include <seqan3/core/detail/debug_stream_range.hpp>
#include <seqan3/core/detail/debug_stream_type.hpp>
#include <seqan3/core/detail/debug_stream_optional.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>

namespace seqan3::detail
{

/*!\brief A debug matrix to wrap alignment matrices and sequences and make them printable together.
 * \ingroup alignment_matrix
 * \implements seqan3::detail::matrix
 * \tparam matrix_t          An alignment matrix; Must model seqan3::detail::matrix.
 * \tparam first_sequence_t  The type of the first sequence; If no sequences are given this is std::nullopt_t.
 * \tparam second_sequence_t The type of the second sequence; If no sequences are given this is std::nullopt_t.
 *
 * \details
 *
 * This debug matrix allows you to print an alignment matrix (e.g. score or trace matrix) combined with two sequences.
 *
 * \cond DEV
 * This type is used internally
 *   * to print seqan3::detail::matrix's and,
 *   * to compare alignment matrices in in test cases.
 * \endcond
 *
 * # Score matrix example
 *
 * \include test/snippet/alignment/matrix/debug_matrix_score.cpp
 *
 * ### Output
 * \include test/snippet/alignment/matrix/debug_matrix_score.out
 *
 * # Trace matrix example
 *
 * \include test/snippet/alignment/matrix/debug_matrix_trace.cpp
 *
 * ### Output
 * \include test/snippet/alignment/matrix/debug_matrix_trace.out
 */
template <matrix matrix_t, typename first_sequence_t = std::nullopt_t, typename second_sequence_t = std::nullopt_t>
class debug_matrix
{
protected:
    //!\brief Whether the current debug_matrix was given a first_sequence.
    static constexpr bool has_first_sequence = !std::is_same_v<std::decay_t<first_sequence_t>, std::nullopt_t>;
    //!\brief Whether the current debug_matrix was given a second_sequence.
    static constexpr bool has_second_sequence = !std::is_same_v<std::decay_t<second_sequence_t>, std::nullopt_t>;
    //!\brief The entry type.
    using entry_t = typename std::remove_reference_t<matrix_t>::value_type;
    //!\brief Whether the value_type is trace_directions.
    static constexpr bool is_traceback_matrix = std::is_same_v<std::decay_t<entry_t>, trace_directions>;
    //!\brief Whether a score matrix already returns std::optional scores. (Where std::nullopt means
    //!       unset/invalid/infinite score)
    static constexpr bool is_optional_score = is_type_specialisation_of_v<entry_t, std::optional>;
public:

    /*!\name Associated types
     * \{
     */
    //!\copydoc seqan3::detail::matrix::value_type
    using value_type = std::conditional_t<is_traceback_matrix || is_optional_score,
                                          entry_t,
                                          std::optional<entry_t>>;
    //!\copydoc seqan3::detail::matrix::reference
    using reference = value_type;
    //!\brief The const reference type.
    using const_reference = reference;
    //!\brief The size type.
    using size_type = typename std::remove_reference_t<matrix_t>::size_type;
    //!\}

    /*!\name Constructors, destructor and assignment
     * \{
     */
    debug_matrix() = default; //!< Defaulted
    debug_matrix(debug_matrix const &) = default; //!< Defaulted
    debug_matrix(debug_matrix &&) = default; //!< Defaulted
    debug_matrix & operator=(debug_matrix const &) = default; //!< Defaulted
    debug_matrix & operator=(debug_matrix &&) = default; //!< Defaulted
    ~debug_matrix() = default;  //!< Defaulted

    /*!\brief Construct the matrix out of an existing matrix.
     * \param matrix An alignment matrix; Must model seqan3::detail::matrix.
     */
    debug_matrix(matrix_t matrix)
        : debug_matrix(std::forward<matrix_t>(matrix), std::nullopt, std::nullopt)
    {}

    /*!\brief Construct the matrix out of an existing matrix and two sequences.
     * \param matrix An alignment matrix; Must model seqan3::detail::matrix.
     * \param first_sequence The first sequence of the sequence alignment.
     * \param second_sequence The second sequence of the sequence alignment.
     */
    debug_matrix(matrix_t matrix, first_sequence_t first_sequence, second_sequence_t second_sequence)
        : _matrix{std::forward<matrix_t>(matrix)},
          _first_sequence{std::forward<first_sequence_t>(first_sequence)},
          _second_sequence{std::forward<second_sequence_t>(second_sequence)}
    {
        if constexpr(has_first_sequence)
        {
            assert(_matrix.cols() <= _first_sequence.size() + 1u);
        }

        if constexpr(has_second_sequence)
        {
            assert(_matrix.rows() <= _second_sequence.size() + 1u);
        }
    }
    //!\}

    //!\copydoc seqan3::detail::matrix::rows
    size_t rows() const noexcept
    {
        if (!_transpose)
            return _rows.value_or(_matrix.rows());
        else
            return _cols.value_or(_matrix.cols());
    }

    //!\copydoc seqan3::detail::matrix::cols
    size_t cols() const noexcept
    {
        if (!_transpose)
            return _cols.value_or(_matrix.cols());
        else
            return _rows.value_or(_matrix.rows());
    }

    //!\copydoc seqan3::detail::debug_matrix::_first_sequence
    first_sequence_t const & first_sequence() const noexcept
    {
        if (!_transpose)
            return _first_sequence;
        else
            return _second_sequence;
    }

    //!\copydoc seqan3::detail::debug_matrix::_second_sequence
    second_sequence_t const & second_sequence() const noexcept
    {
        if (!_transpose)
            return _second_sequence;
        else
            return _first_sequence;
    }

    //!\copydoc seqan3::detail::matrix::at
    const_reference at(matrix_coordinate const & coordinate) const noexcept
    {
        size_t row = coordinate.row;
        size_t col = coordinate.col;

        assert(row < rows() && col < cols());

        row_index_type const _row{!_transpose ? row : col};
        column_index_type const _col{!_transpose ? col : row};
        row_index_type const _mask_row{_transpose == _transpose_mask ? row : col};
        column_index_type const _mask_col{_transpose == _transpose_mask ? col : row};

        if (!_masking_matrix.has_value() || _masking_matrix.value().at({_mask_row, _mask_col}))
        {
            entry_t const & entry = _matrix.at({_row, _col});

            if (!is_traceback_matrix || !_transpose)
                return entry;

            if constexpr(is_traceback_matrix)
            {
                trace_directions reverse{};
                if ((entry & trace_directions::left) == trace_directions::left)
                    reverse |= trace_directions::up;
                if ((entry & trace_directions::up) == trace_directions::up)
                    reverse |= trace_directions::left;
                if ((entry & trace_directions::diagonal) == trace_directions::diagonal)
                    reverse |= trace_directions::diagonal;
                return reverse;
            }
        }

        if constexpr(is_traceback_matrix)
            return trace_directions::none;
        else
            return std::nullopt;
    }

    /*!\brief Masks entries out of the current matrix. This operations changes the way `this.at(i, j)` will operate.
     * If `masking_matrix.at(i,j)` returns true `this.at(i, j)` will operate as usual.
     * But, if false `this.at(i, j)` will return std::nullopt.
     * \param masking_matrix \copydoc _masking_matrix
     * \returns *this
     */
    debug_matrix & mask_matrix(row_wise_matrix<bool> masking_matrix) noexcept
    {
        assert(masking_matrix.rows() == rows());
        assert(masking_matrix.cols() == cols());
        _transpose_mask = _transpose;
        _masking_matrix = std::move(masking_matrix);
        return *this;
    }

    /*!\brief Creates the masking_matrix out of the given masking_vector and calls #mask_matrix(row_wise_matrix<bool>)
     * \param masking_vector The masking vector to construct the masking_matrix.
     * \returns *this
     */
    debug_matrix & mask_matrix(std::vector<bool> masking_vector) noexcept
    {
        return mask_matrix(row_wise_matrix<bool>{number_rows{rows()},
                                                 number_cols{cols()},
                                                 std::move(masking_vector)});
    }

    /*!\brief Limits the view port of the current matrix.
     * \param new_rows \copydoc _rows
     * \param new_cols \copydoc _cols
     * \returns *this
     */
    debug_matrix & sub_matrix(size_t const new_rows, size_t const new_cols) noexcept
    {
        assert(new_rows <= rows());
        assert(new_cols <= cols());
        if (!_transpose)
        {
            _rows = new_rows;
            _cols = new_cols;
        }
        else
        {
            _rows = new_cols;
            _cols = new_rows;
        }
        return *this;
    }

    /*!\brief Transposes the current matrix.
     * \returns *this
     */
    debug_matrix & transpose_matrix() noexcept
    {
        _transpose = !_transpose;
        return *this;
    }

protected:
    //!\cond
    struct format_type; // forward declaration
    //!\endcond

public:
    /*!\brief Prints this matrix into the given stream.
     * \param cout  The stream to print to.
     * \param flags Modify the way the matrix is printed.
     *
     * \details
     *
     * The matrix will be printed with unicode characters if seqan3::fmtflags2::utf8 is set in the flags. Ascii
     * otherwise.
     */
    template <typename ostream_t>
    void stream_matrix(ostream_t & cout, fmtflags2 const flags) const noexcept
    {
        format_type const & symbols = (flags & fmtflags2::utf8) == fmtflags2::utf8 ? unicode : csv;
        size_t const column_width = this->column_width.has_value() ?
                                    this->column_width.value() : auto_column_width(flags);

        auto char_first_sequence = [&]([[maybe_unused]] size_t const i) -> std::string
        {
            if constexpr(!has_first_sequence)
                return " ";
            else
                return as_string(first_sequence()[i], flags);
        };

        auto char_second_sequence = [&]([[maybe_unused]] size_t const i) -> std::string
        {
            if constexpr(!has_second_sequence)
                return " ";
            else
                return as_string(second_sequence()[i], flags);
        };

        auto print_cell = [&](std::string const & symbol)
        {
            // deal with unicode chars that mess up std::setw
            size_t const length_bytes = symbol.size();
            size_t const length = unicode_str_length(symbol);
            size_t const offset = length_bytes - length;

            cout << std::left
                 << std::setw(column_width + offset)
                 << symbol
                 << symbols.col_sep;
        };

        auto print_first_cell = [&](std::string const & symbol)
        {
            cout << symbol << symbols.col_sep;
        };

        // |_|d|a|t|a|b|a|s|e|
        auto print_first_row = [&]
        {
            print_first_cell(" ");
            print_cell(symbols.epsilon);

            for (size_t col = 0; col < cols() - 1; ++col)
                print_cell(char_first_sequence(col));

            cout << "\n";
        };

        // |-|-|-|-|-|-|-|-|-|
        auto print_divider = [&]
        {
            cout << " " << symbols.row_col_sep;
            for (size_t col = 0; col < cols(); ++col)
            {
                for (size_t i = 0; i < column_width; ++i)
                    cout << symbols.row_sep;

                cout << symbols.row_col_sep;
            }
            cout << "\n";
        };

        print_first_row();
        for (size_t row = 0; row < rows(); ++row)
        {
            if (symbols.row_sep[0] != '\0')
                print_divider();

            // one query letter + one row of scores / traces
            if (row == 0)
                print_first_cell(symbols.epsilon);
            else
                print_first_cell(char_second_sequence(row - 1));

            for (size_t col = 0; col < cols(); ++col)
                print_cell(entry_at({row_index_type{row}, column_index_type{col}}, flags));

            cout << "\n";
        }
    }

    //!\brief Determines the largest width of all entries in the matrix, e.g. `-152` has width 4.
    size_t auto_column_width(fmtflags2 const flags) const noexcept
    {
        size_t col_width = 1;
        for (size_t row = 0; row < rows(); ++row)
            for (size_t col = 0; col < cols(); ++col)
                col_width = std::max(col_width,
                                     unicode_str_length(entry_at({row_index_type{row}, column_index_type{col}}, flags)));

        return col_width;
    }

protected:
    //!\brief Same as at(*coordinate*), but as string.
    std::string entry_at(matrix_coordinate const coordinate, fmtflags2 flags) const noexcept
    {
        format_type const & symbols = (flags & fmtflags2::utf8) == fmtflags2::utf8 ? unicode : csv;

        value_type const & entry = at(coordinate);
        if (!is_traceback_matrix && entry == matrix_inf<value_type>)
            return symbols.inf;

        return as_string(entry, flags);
    }

    //!\brief Convert a value into a std::string.
    template <typename value_type>
    static std::string as_string(value_type && entry, fmtflags2 const flags) noexcept
    {
        std::stringstream strstream;
        debug_stream_type stream{strstream};
        stream << flags << entry;
        return strstream.str();
    }

    //!\brief The length of the *str* (traceback symbols are unicode aware).
    //!\sa https://en.wikipedia.org/wiki/UTF-8 for encoding details
    static size_t unicode_str_length(std::string const & str) noexcept
    {
        size_t length = 0u;
        for (auto it = str.cbegin(); it < str.cend(); ++it, ++length)
        {
            uint8_t v = *it;
            if ((v & 0b11100000) == 0b11000000)
                ++it;
            else if ((v & 0b11110000) == 0b11100000)
                it += 2;
            else if ((v & 0b11111000) == 0b11110000)
                it += 3;
        }
        return length;
    }

    //!\brief Format used by seqan3::detail::debug_matrix.
    struct format_type
    {
        //!\brief Epsilon symbol (a single symbol)
        char const * epsilon{};
        //!\brief Column separator symbol (a single symbol)
        char const * col_sep{};
        //!\brief Row separator symbol (a single symbol)
        char const * row_sep{};
        //!\brief Row separator symbol (a single symbol)
        char const * row_col_sep{};
        //!\brief Infinity symbol (a single symbol)
        char const * inf{};
    };

    //!\brief The format when printing to a ascii stream.
    static constexpr format_type csv{" ", ";", "", "", ""};
    //!\brief The format when printing to a unicode stream.
    static constexpr format_type unicode{"ε", "║", "═", "╬", "∞"};

public:
    //!\brief What is the width (number of chars) of an entry. Defaults to #auto_column_width.
    std::optional<size_t> column_width{std::nullopt};

protected:
    //!\brief The matrix
    matrix_t _matrix;
    //!\brief The first sequence of the sequence alignment.
    first_sequence_t _first_sequence;
    //!\brief The second sequence of the sequence alignment.
    second_sequence_t _second_sequence;
    //!\brief The number of rows the debug matrix should have. Must be at most the size of the original matrix.
    std::optional<size_t> _rows{};
    //!\brief The number of columns the debug matrix should have. Must be at most the size of the original matrix.
    std::optional<size_t> _cols{};
    //!\brief The masking matrix.
    std::optional<row_wise_matrix<bool>> _masking_matrix{};
    //!\brief Whether the current matrix should be transposed.
    bool _transpose{};
    //!\brief Whether the masking matrix should be transposed.
    bool _transpose_mask{};
};

/*!\name Type deduction guides
 * \relates seqan3::detail::debug_matrix
 * \{
 */
//!\brief The type deduction guide for the constructor seqan3::detail::debug_matrix(matrix_t)
template <matrix matrix_t>
debug_matrix(matrix_t &&)
   -> debug_matrix<matrix_t>;

//!\brief The type deduction guide for the constructor
//!       seqan3::detail::debug_matrix(matrix_t, first_sequence_t, second_sequence_t)
template <matrix matrix_t, typename first_sequence_t, typename second_sequence_t>
debug_matrix(matrix_t &&, first_sequence_t &&, second_sequence_t &&)
   -> debug_matrix<matrix_t, first_sequence_t, second_sequence_t>;
//!\}

} // namespace seqan3::detail

namespace seqan3
{
/*!\brief An alignment matrix can be printed to the seqan3::debug_stream.
 * \tparam alignment_matrix_t Type of the alignment matrix to be printed; must model seqan3::detail::matrix.
 * \param s The seqan3::debug_stream.
 * \param matrix The alignment matrix.
 * \relates seqan3::debug_stream_type
 *
 * \details
 *
 * This prints out an alignment matrix which can be a score matrix or a trace matrix.
 */
template <detail::matrix alignment_matrix_t, typename char_t>
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & s, alignment_matrix_t && matrix)
{
    detail::debug_matrix debug{std::forward<alignment_matrix_t>(matrix)};

    std::stringstream sstream{};
    debug.stream_matrix(sstream, s.flags2());
    s << sstream.str();
    return s;
}

//!\overload
template <std::ranges::input_range alignment_matrix_t, typename char_t>
//!\cond
    requires detail::debug_stream_range_guard<alignment_matrix_t> && detail::matrix<alignment_matrix_t>
//!\endcond
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & s, alignment_matrix_t && matrix)
{
    return s << detail::debug_matrix{std::forward<alignment_matrix_t>(matrix)};
}

} // namespace seqan3
