// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Contains seqan3::detail::alignment_matrix_formatter and seqan3::detail::alignment_matrix_format.
 */

#pragma once

#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>

#include <seqan3/alignment/matrix/matrix_concept.hpp>
#include <seqan3/alignment/matrix/alignment_trace_matrix.hpp>
#include <seqan3/core/debug_stream.hpp>

namespace seqan3::detail
{

/*!\brief Format used by seqan3::detail::alignment_matrix_formatter
 * \ingroup alignment_matrix
 *
 * \details
 *
 * With seqan3::detail::alignment_matrix_format you can style:
 *  * the epsilon symbol (alignment_matrix_format::epsilon)
 *  * the column symbol that separates each cell in a row (alignment_matrix_format::col_sep)
 *  * the row symbol that divides each row (alignment_matrix_format::row_sep)
 *  * the column symbol that comes after the row symbol (alignment_matrix_format::row_col_sep)
 *  * the traces symbols of a traceback matrix (alignment_matrix_format::trace_dir)
 *
 * # Example
 * \snippet test/snippet/alignment/matrix/alignment_matrix_formatter.cpp alignment_matrix_format
 *
 * ### Output
 * \snippet test/snippet/alignment/matrix/alignment_matrix_formatter.out alignment_matrix_format::out
 */
struct alignment_matrix_format
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

    /*!\brief Eight symbols for each combination of directions a trace can have (each entry can have multiple symbols).
     *
     * * 1st bit: **D** = diagonal
     * * 2nd bit: **U** = up
     * * 3th bit: **L** = left
     *
     * | i     | 1st bit | 2nd bit | 3rd bit | trace_dir[i] |
     * | :---: |---------|---------|---------| :----------: |
     * | **0** | 0       | 0       | 0       | **No dir**   |
     * | **1** | 1       | 0       | 0       | **D**        |
     * | **2** | 0       | 1       | 0       | **U**        |
     * | **3** | 1       | 1       | 0       | **DU**       |
     * | **4** | 0       | 0       | 1       | **L**        |
     * | **5** | 1       | 0       | 1       | **DL**       |
     * | **6** | 0       | 1       | 1       | **UL**       |
     * | **7** | 1       | 1       | 1       | **DUL**      |
     *
     */
    char const * trace_dir[8]{};

    /*!\brief The CSV format that makes it easy to export the matrix
     * \sa https://en.wikipedia.org/wiki/Comma-separated_values
     *
     *
     * ### Score output
     * \snippet test/snippet/alignment/matrix/alignment_matrix_formatter.out alignment_matrix_format::csv#score
     *
     * ### Trace output
     * \snippet test/snippet/alignment/matrix/alignment_matrix_formatter.out alignment_matrix_format::csv#trace
     */
    static const alignment_matrix_format csv;

    /*!\brief A format that uses only ascii symbols.
     * \sa https://en.wikipedia.org/wiki/Ascii
     *
     *
     * ### Score output
     * \snippet test/snippet/alignment/matrix/alignment_matrix_formatter.out alignment_matrix_format::ascii#score
     *
     * ### Trace output
     * \snippet test/snippet/alignment/matrix/alignment_matrix_formatter.out alignment_matrix_format::ascii#trace
     */
    static const alignment_matrix_format ascii;

    /*!\brief A format that uses unicode block symbols.
     * \sa https://en.wikipedia.org/wiki/Unicode
     *
     *
     * ### Score output
     * \snippet test/snippet/alignment/matrix/alignment_matrix_formatter.out alignment_matrix_format::unicode_block#score
     *
     * ### Trace output
     * \snippet test/snippet/alignment/matrix/alignment_matrix_formatter.out alignment_matrix_format::unicode_block#trace
     */
    static const alignment_matrix_format unicode_block;

    /*!\brief A format that uses unicode braille symbols.
     * \sa https://en.wikipedia.org/wiki/Unicode
     *
     *
     * ### Score output
     * \snippet test/snippet/alignment/matrix/alignment_matrix_formatter.out alignment_matrix_format::unicode_braille#score
     *
     * ### Trace output
     * \snippet test/snippet/alignment/matrix/alignment_matrix_formatter.out alignment_matrix_format::unicode_braille#trace
     */
    static const alignment_matrix_format unicode_braille;

   /*!\brief A format that uses unicode arrow symbols.
    * \sa https://en.wikipedia.org/wiki/Unicode
    *
    *
    * ### Score output
    * \snippet test/snippet/alignment/matrix/alignment_matrix_formatter.out alignment_matrix_format::unicode_arrows#score
    *
    * ### Trace output
    * \snippet test/snippet/alignment/matrix/alignment_matrix_formatter.out alignment_matrix_format::unicode_arrows#trace
    */
    static const alignment_matrix_format unicode_arrows;
};

constexpr alignment_matrix_format alignment_matrix_format::csv
{
    " ", ";", "", "", "",
    {"N","D","U","DU","L","DL","UL","DUL"}
};

constexpr alignment_matrix_format alignment_matrix_format::ascii
{
    " ", "|", "-", "/", "INF",
    {" ","D","U","DU","L","DL","UL","DUL"}
};

constexpr alignment_matrix_format alignment_matrix_format::unicode_block
{
    "ε", "║", "═", "╬", "∞",
    {"█","▘","▝","▀","▖","▌","▞","▛"}
};

constexpr alignment_matrix_format alignment_matrix_format::unicode_braille
{
    "ε", "║", "═", "╬", "∞",
    {"⠀","⠁","⠈","⠉","⠄","⠅","⠌","⠍"}
};

constexpr alignment_matrix_format alignment_matrix_format::unicode_arrows
{
    "ε", "║", "═", "╬", "∞",
    {"↺","↖","↑","↖↑","←","↖←","↑←","↖↑←"}
};

/*!\brief Formats and prints trace and score matrices that satisfy the
 *        seqan3::detail::Matrix.
 * \ingroup alignment_matrix
 * \tparam alignment_matrix_t The matrix class which satisfies the seqan3::detail::Matrix
 *
 * \details
 *
 * # Example
 *
 * \snippet test/snippet/alignment/matrix/alignment_score_matrix.cpp code
 *
 * ### Output
 * \include test/snippet/alignment/matrix/alignment_score_matrix.out
 */
template <Matrix alignment_matrix_t>
class alignment_matrix_formatter
{
public:
    //!\brief The type of the #matrix
    using alignment_matrix_type = alignment_matrix_t;

private:
    //!\brief The matrix to format.
    alignment_matrix_type const & matrix;

public:
    //!\brief The actual format used by the formatter.
    alignment_matrix_format symbols;

    /*!\name Constructors, destructor and assignment
     * \{
     */
     alignment_matrix_formatter() = delete;                                                //!< Matrix taken as ref.
     alignment_matrix_formatter(alignment_matrix_formatter const &) = default;             //!< Defaulted
     alignment_matrix_formatter(alignment_matrix_formatter &&) = default;                  //!< Defaulted
     alignment_matrix_formatter & operator=(alignment_matrix_formatter const &) = default; //!< Defaulted
     alignment_matrix_formatter & operator=(alignment_matrix_formatter &&) = default;      //!< Defaulted
     ~alignment_matrix_formatter() = default;                                              //!< Defaulted

    /*!\brief Construct the matrix out of the *entries*, the *rows*,
     *        and the *cols*.
     * \param _matrix  \copydoc matrix
     * \param _symbols \copydoc symbols
     */
    alignment_matrix_formatter(alignment_matrix_type const & _matrix, alignment_matrix_format _symbols = alignment_matrix_format::unicode_arrows)
        : matrix{_matrix}, symbols{_symbols}
    {}
    //!\}

    //!\brief \copydoc seqan3::detail::Matrix::entry_type
    using entry_type = typename alignment_matrix_type::entry_type;

    //!\brief Whether #alignment_matrix_type is a traceback matrix.
    static constexpr bool is_traceback_matrix = std::is_same_v<entry_type, trace_directions>;

    //!\brief Determines the largest width of all entries in the #matrix,
    //!       e.g. `-152` has width 4.
    size_t auto_width() const noexcept
    {
        size_t col_width = 1;
        for (size_t row = 0; row < matrix.rows(); ++row)
            for (size_t col = 0; col < matrix.cols(); ++col)
                col_width = std::max(col_width, unicode_str_length(entry_at(row, col)));
        return col_width;
    }

    //!\brief Print the formatted #matrix to std::cout.
    //!\param[in]  database        the database sequence
    //!\param[in]  query           the query sequence
    //!\param[in]  column_width    width of each cell, std::nullopt defaults to auto_width()
    template <typename database_t, typename query_t>
    void format(database_t && database, query_t && query, std::optional<size_t> column_width = std::nullopt) const noexcept
    {
        format(std::forward<database_t>(database), std::forward<query_t>(query), std::cout, column_width);
    }

    /*!\brief Print the formatted #matrix to the *cout* stream.
     * \tparam         char_t          char type of std::basic_ostream
     * \tparam         traits_t        traits type of std::basic_ostream
     * \param[in]      database        the database sequence
     * \param[in]      query           the query sequence
     * \param[in,out]  cout            print formatted #matrix into this ostream
     * \param[in]      column_width    width of each cell, std::nullopt defaults to auto_width()
     */
    template <typename database_t, typename query_t, typename oostream_t>
    void format(database_t && database, query_t && query, oostream_t & cout, std::optional<size_t> const column_width) const noexcept
    {
        size_t const _column_width = column_width.has_value() ? column_width.value() : auto_width();

        auto print_cell = [&](std::string const & symbol)
        {
            // deal with unicode chars that mess up std::setw
            size_t const length_bytes = unicode_str_length_bytes(symbol);
            size_t const length = unicode_str_length(symbol);
            size_t const offset = length_bytes - length;

            cout << std::left
                 << std::setw(_column_width + offset)
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

            for (size_t col = 0; col < matrix.cols() - 1; ++col)
                print_cell(as_string(database[col]));
            cout << "\n";
        };

        // |-|-|-|-|-|-|-|-|-|
        auto print_divider = [&]
        {
            cout << " " << symbols.row_col_sep;
            for (size_t col = 0; col < matrix.cols(); ++col)
            {
                for (size_t i = 0; i < _column_width; ++i)
                    cout << symbols.row_sep;
                cout << symbols.row_col_sep;
            }
            cout << "\n";
        };

        print_first_row();
        for (size_t row = 0; row < matrix.rows(); ++row)
        {
            if (symbols.row_sep[0] != '\0')
                print_divider();

            // one query letter + one row of scores / traces
            if (row == 0)
                print_first_cell(symbols.epsilon);
            else
                print_first_cell(as_string(query[row - 1]));
            for (size_t col = 0; col < matrix.cols(); ++col)
                print_cell(entry_at(row, col));
            cout << "\n";
        }
    }

private:
    //!\brief Same as #matrix\.at(*row*, *col*), but converts the value to
    //!       a trace symbol (alignment_matrix_format::trace_dir) if the
    //!       #matrix is a traceback matrix.
    std::string entry_at(size_t const row, size_t const col) const noexcept
    {
        if constexpr(is_traceback_matrix)
        {
            trace_directions direction = matrix.at(row, col);
            return symbols.trace_dir[(size_t)(direction) % 8u];
        }
        else
        {
            entry_type entry = matrix.at(row, col);
            if (entry == matrix_inf<entry_type>)
                return as_string(symbols.inf);
            return as_string(entry);
        }
    }

    //!\brief Convert a matrix entry into a std::string
    template <typename entry_type>
    static std::string as_string(entry_type && entry) noexcept
    {
        std::stringstream strstream;
        debug_stream_type stream{strstream};
        stream << entry;
        return strstream.str();
    }

    //!\brief The length of the *str* (traceback symbols are unicode aware)
    //!\sa https://en.wikipedia.org/wiki/UTF-8 for encoding details
    static size_t unicode_str_length(std::string const & str) noexcept
    {
        size_t length = 0u;
        for (auto it = str.cbegin(), it_end = str.cend(); it < it_end; ++it, ++length)
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

    //!\brief The number of bytes the *str* uses
    static size_t unicode_str_length_bytes(std::string const & str) noexcept
    {
        return str.length();
    }

    //!\brief Befriend test case
    friend struct matrix_formatter_test;
};

/*!\name Type deduction guides
 * \relates seqan3::detail::alignment_matrix_formatter
 * \{
 */

//!\brief Deduce the matrix type from the provided argument.
template<typename alignment_matrix_t>
alignment_matrix_formatter(alignment_matrix_t const &) -> alignment_matrix_formatter<alignment_matrix_t>;

//!\brief Deduce the matrix type from the provided arguments.
template<typename alignment_matrix_t>
alignment_matrix_formatter(alignment_matrix_t const &, alignment_matrix_format) -> alignment_matrix_formatter<alignment_matrix_t>;
//!\}

} // namespace seqan3::detail
