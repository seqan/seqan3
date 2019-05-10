// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Contains the declaration of seqan3::detail::debug_matrix.
 */

#pragma once

#include <iomanip>

#include <seqan3/alignment/matrix/matrix_concept.hpp>
#include <seqan3/alignment/matrix/row_wise_matrix.hpp>
#include <seqan3/alignment/matrix/trace_directions.hpp>
#include <seqan3/core/debug_stream.hpp>

namespace seqan3::detail
{

/*!\brief A debug matrix to wrap alignment matrices and sequences and make them printable together.
 * \ingroup alignment_matrix
 * \implements seqan3::detail::Matrix
 * \tparam matrix_t    An alignment matrix; Must model seqan3::detail::Matrix.
 * \tparam sequence1_t The type of the first sequence; If no sequences are given this is std::nullopt_t.
 * \tparam sequence2_t The type of the second sequence; If no sequences are given this is std::nullopt_t.
 *
 *
 * This debug matrix allows you to print an alignment matrix (e.g. score or trace matrix) combined with two sequences.
 *
 * \cond DEV
 * This type is used internally
 *   * to print seqan3::detail::Matrix's and,
 *   * to compare alignment matrices in in test cases.
 * \endcond
 *
 * # Score Matrix Example
 *
 * \snippet test/snippet/alignment/matrix/debug_matrix.cpp score_matrix
 *
 * ### Output
 * \snippet test/snippet/alignment/matrix/debug_matrix.out score_matrix::out
 *
 * # Trace Matrix Example
 *
 * \snippet test/snippet/alignment/matrix/debug_matrix.cpp trace_matrix
 *
 * ### Output
 * \snippet test/snippet/alignment/matrix/debug_matrix.out trace_matrix::out
 */
template <Matrix matrix_t, typename sequence1_t = std::nullopt_t, typename sequence2_t = std::nullopt_t>
class debug_matrix
{
public:
    //!\copydoc seqan3::detail::Matrix::entry_type
    using entry_type = typename std::remove_reference_t<matrix_t>::entry_type;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    debug_matrix() = default; //!< Defaulted
    debug_matrix(debug_matrix const &) = default; //!< Defaulted
    debug_matrix(debug_matrix &&) = default; //!< Defaulted
    debug_matrix & operator=(debug_matrix const &) = default; //!< Defaulted
    debug_matrix & operator=(debug_matrix &&) = default; //!< Defaulted
    ~debug_matrix() = default;  //!< Defaulted

    /*!\brief Construct the matrix out of the *entries*, the *rows*, and the *cols*.
     * \param entries The entry values as a flat std::vector <#entry_type>.
     * \param rows    The number of rows.
     * \param cols    The number of columns.
     */
    debug_matrix(std::vector<entry_type> entries, size_t const rows, size_t const cols)
        : _matrix{std::move(entries), rows, cols},
          _sequence1{std::nullopt},
          _sequence2{std::nullopt}
    {}

    /*!\brief Construct the matrix out of the *entries*, the *rows*, and the *cols* and the corresponding sequences.
     * \param entries The entry values as a flat std::vector <#entry_type>.
     * \param rows    The number of rows.
     * \param cols    The number of columns.
     * \param sequence1 The first sequence of the sequence alignment.
     * \param sequence2 The second sequence of the sequence alignment.
     */
    debug_matrix(std::vector<entry_type> entries, size_t const rows, size_t const cols,
                 sequence1_t sequence1, sequence2_t sequence2)
        : _matrix{std::move(entries), rows, cols},
          _sequence1{std::forward<sequence1_t>(sequence1)},
          _sequence2{std::forward<sequence2_t>(sequence2)}
    {
        assert(cols == _sequence1.size() + 1u);
        assert(rows == _sequence2.size() + 1u);
    }

    /*!\brief Construct the matrix out of an existing matrix.
     * \param matrix An alignment matrix; Must model seqan3::detail::Matrix.
     */
    debug_matrix(matrix_t matrix)
        : _matrix{std::forward<matrix_t>(matrix)},
          _sequence1{std::nullopt},
          _sequence2{std::nullopt}
    {}

    /*!\brief Construct the matrix out of an existing matrix and two sequences.
     * \param matrix An alignment matrix; Must model seqan3::detail::Matrix.
     * \param sequence1 The first sequence of the sequence alignment.
     * \param sequence2 The second sequence of the sequence alignment.
     */
    debug_matrix(matrix_t matrix, sequence1_t sequence1, sequence2_t sequence2)
        : _matrix{std::forward<matrix_t>(matrix)},
          _sequence1{std::forward<sequence1_t>(sequence1)},
          _sequence2{std::forward<sequence2_t>(sequence2)}
    {
        assert(_matrix.cols() == _sequence1.size() + 1u);
        assert(_matrix.rows() == _sequence2.size() + 1u);
    }
    //!\}

    //!\copydoc seqan3::detail::Matrix::rows
    size_t rows() const noexcept
    {
        return _matrix.rows();
    }

    //!\copydoc seqan3::detail::Matrix::cols
    size_t cols() const noexcept
    {
        return _matrix.cols();
    }

    //!\copydoc seqan3::detail::debug_matrix::_sequence1
    sequence1_t const & sequence1() const noexcept
    {
        return _sequence1;
    }

    //!\copydoc seqan3::detail::debug_matrix::_sequence2
    sequence2_t const & sequence2() const noexcept
    {
        return _sequence2;
    }

    //!\copydoc seqan3::detail::Matrix::at
    entry_type at(size_t const row, size_t const col) const noexcept
    {
        assert(row < rows() && col < cols());
        return _matrix.at(row, col);
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
     *
     * The matrix will be printed out with unicode characters if seqan3::fmtflags2::utf8 is set in the flags. Ascii
     * otherwise.
     */
    template <typename oostream_t>
    void print(oostream_t & cout, fmtflags2 const flags) const noexcept
    {
        format_type const & symbols = (flags & fmtflags2::utf8) == fmtflags2::utf8 ? unicode : csv;
        size_t const column_width = this->column_width.has_value() ?
                                    this->column_width.value() : auto_column_width(flags);

        auto char_sequence1 = [&]([[maybe_unused]] size_t const i) -> std::string
        {
            if constexpr(std::is_same_v<sequence1_t, std::nullopt_t>)
                return " ";
            else
                return as_string(sequence1()[i], flags);
        };

        auto char_sequence2 = [&]([[maybe_unused]] size_t const i) -> std::string
        {
            if constexpr(std::is_same_v<sequence2_t, std::nullopt_t>)
                return " ";
            else
                return as_string(sequence2()[i], flags);
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
                print_cell(char_sequence1(col));
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
                print_first_cell(char_sequence2(row - 1));

            for (size_t col = 0; col < cols(); ++col)
                print_cell(entry_at(row, col, flags));
            cout << "\n";
        }
    }

    //!\brief Determines the largest width of all entries in the matrix, e.g. `-152` has width 4.
    size_t auto_column_width(fmtflags2 const flags) const noexcept
    {
        size_t col_width = 1;
        for (size_t row = 0; row < rows(); ++row)
            for (size_t col = 0; col < cols(); ++col)
                col_width = std::max(col_width, unicode_str_length(entry_at(row, col, flags)));
        return col_width;
    }

protected:
    //!\brief Whether the entry_type is trace_directions.
    static constexpr bool is_traceback_matrix = std::is_same_v<std::decay_t<entry_type>, trace_directions>;

    //!\brief Same as at(*row*, *col*), but as string.
    std::string entry_at(size_t const row, size_t const col, fmtflags2 flags) const noexcept
    {
        format_type const & symbols = (flags & fmtflags2::utf8) == fmtflags2::utf8 ? unicode : csv;

        entry_type const & entry = at(row, col);
        if (!is_traceback_matrix && entry == matrix_inf<entry_type>)
            return symbols.inf;
        return as_string(entry, flags);
    }

    //!\brief Convert a value into a std::string.
    template <typename entry_type>
    static std::string as_string(entry_type && entry, fmtflags2 const flags) noexcept
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
        //!\brief epsilon symbol (a single symbol)
        char const * epsilon{};
        //!\brief column separator symbol (a single symbol)
        char const * col_sep{};
        //!\brief row separator symbol (a single symbol)
        char const * row_sep{};
        //!\brief row separator symbol (a single symbol)
        char const * row_col_sep{};
        //!\brief infinity symbol (a single symbol)
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
    sequence1_t _sequence1;
    //!\brief The second sequence of the sequence alignment.
    sequence2_t _sequence2;
};

/*!\name Type deduction guides
 * \relates seqan3::detail::debug_matrix
 * \{
 */
//!\brief The type deduction guide for the constructor
//!seqan3::detail::debug_matrix(std::vector<entry_type>, size_t const, size_t const)
template <typename entry_type>
debug_matrix(std::vector<entry_type>, size_t const, const size_t)
    -> debug_matrix<row_wise_matrix<entry_type>>;

//!\brief The type deduction guide for the constructor
//!seqan3::detail::debug_matrix(std::vector<entry_type>, size_t const, size_t const, sequence1_t, sequence2_t)
template <typename entry_type, typename sequence1_t, typename sequence2_t>
debug_matrix(std::vector<entry_type>, size_t const, const size_t, sequence1_t &&, sequence2_t &&)
   -> debug_matrix<row_wise_matrix<entry_type>, sequence1_t, sequence2_t>;


//!\brief The type deduction guide for the constructor seqan3::detail::debug_matrix(matrix_t)
template <Matrix matrix_t>
debug_matrix(matrix_t &&)
   -> debug_matrix<matrix_t>;

//!\brief The type deduction guide for the constructor seqan3::detail::debug_matrix(matrix_t, sequence1_t, sequence2_t)
template <Matrix matrix_t, typename sequence1_t, typename sequence2_t>
debug_matrix(matrix_t &&, sequence1_t &&, sequence2_t &&)
   -> debug_matrix<matrix_t, sequence1_t, sequence2_t>;
//!\}

} // namespace seqan3::detail
