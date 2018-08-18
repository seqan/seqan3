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
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Contains seqan3::alignment_matrix_formatter and seqan3::alignment_matrix_format.
 */

#pragma once

#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>

#include <seqan3/alignment/matrix/alignment_matrix_concept.hpp>
#include <seqan3/alignment/matrix/alignment_trace_matrix.hpp>
#include <seqan3/core/metafunction/template_inspection.hpp>

namespace seqan3
{

/*!\brief Format used by seqan3::alignment_matrix_formatter
 * \ingroup alignment_matrix
 *
 *
 * With seqan3::alignment_matrix_format you can style:
 *  * the epsilon symbol (alignment_matrix_format::epsilon)
 *  * the column symbol that separates each cell in a row (alignment_matrix_format::col_sep)
 *  * the row symbol that divides each row (alignment_matrix_format::row_sep)
 *  * the column symbol that comes after the row symbol (alignment_matrix_format::row_col_sep)
 *  * the traces symbols of a traceback matrix (alignment_matrix_format::trace_dir)
 *
 * ## Example
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

    /*!\brief eight symbols for each combination of directions a trace can have (each entry can have multiple symbols)
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

    //!\brief how many bytes needs the epsilon symbol (necessary when using utf-8)
    uint8_t epsilon_bytes{1};

    /*!\brief how many bytes needs every trace_dir symbol (necessary when using utf-8)
     * \attention A restriction is that every symbol needs the same amount of
     * bytes to represent one char.
     */
    uint8_t trace_dir_bytes{1};

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
    " ", ";", "", "",
    {"N","D","U","DU","L","DL","UL","DUL"}
};

constexpr alignment_matrix_format alignment_matrix_format::ascii
{
    " ", "|", "-", "/",
    {" ","D","U","DU","L","DL","UL","DUL"} // 1byte per char
};

constexpr alignment_matrix_format alignment_matrix_format::unicode_block
{
    u8"ε", u8"║", u8"═", u8"╬",
    {u8"█",u8"▘",u8"▝",u8"▀",u8"▖",u8"▌",u8"▞",u8"▛"}, // 3bytes per char
    2, 3
};

constexpr alignment_matrix_format alignment_matrix_format::unicode_braille
{
    u8"ε", u8"║", u8"═", u8"╬",
    {u8"⠀",u8"⠁",u8"⠈",u8"⠉",u8"⠄",u8"⠅",u8"⠌",u8"⠍"}, // 3bytes per char
    2, 3
};

constexpr alignment_matrix_format alignment_matrix_format::unicode_arrows
{
    u8"ε", u8"║", u8"═", u8"╬",
    {u8"↺",u8"↖",u8"↑",u8"↖↑",u8"←",u8"↖←",u8"↑←",u8"↖↑←"}, // 3bytes per char
    2, 3
};

/*!\brief Formats and prints trace and score matrices that satisfy the
 *        seqan3::alignment_matrix_concept.
 * \ingroup alignment_matrix
 * \tparam alignment_matrix_t The matrix class which satisfies the seqan3::alignment_matrix_concept
 *
 *
 * ## Example
 *
 * \snippet test/snippet/alignment/matrix/alignment_score_matrix.cpp code
 *
 * ### Output
 * \include test/snippet/alignment/matrix/alignment_score_matrix.out
 */
template <alignment_matrix_concept alignment_matrix_t>
struct alignment_matrix_formatter
{
    //!\brief The type of the #matrix
    using alignment_matrix_type = alignment_matrix_t;

    //!\brief The matrix to format.
    alignment_matrix_type const & matrix;

    //!\brief The actual format used by the formatter.
    alignment_matrix_format symbols{alignment_matrix_format::unicode_arrows};

    //!\brief Whether #alignment_matrix_type is a traceback matrix.
    //!\hideinitializer
    static constexpr bool is_traceback_matrix = []
    {
        using trace_matrix_t = detail::transfer_template_args_onto_t<alignment_matrix_type, alignment_trace_matrix>;
        return std::is_same_v<alignment_matrix_type, trace_matrix_t>;
    }();

    //!\brief Determines the largest width of all entries in the #matrix,
    //!       e.g. `-152` has width 4.
    std::size_t auto_width() const noexcept
    {
        std::size_t col_width = 1;
        for(unsigned row = 0; row < matrix.rows(); ++row)
            for(unsigned col = 0; col < matrix.cols(); ++col)
                col_width = std::max(col_width, char_length(char_at(row, col)));
        return col_width;
    }

    //!\brief Print the formatted #matrix to std::cout.
    //!\param[in]  column_width    width of each cell, std::nullopt defaults to auto_width()
    void format(std::optional<std::size_t> column_width = std::nullopt) const noexcept
    {
        format(std::cout, column_width);
    }

    /*!\brief Print the formatted #matrix to the \a cout stream.
     * \tparam         char_t          char type of std::basic_ostream
     * \tparam         traits_t        traits type of std::basic_ostream
     * \param[in,out]  cout            print formatted #matrix into this ostream
     * \param[in]      column_width    width of each cell, std::nullopt defaults to auto_width()
     */
    template <typename char_t, typename traits_t>
    void format(std::basic_ostream<char_t, traits_t> & cout, std::optional<std::size_t> column_width) const noexcept
    {
        std::size_t _column_width = column_width.has_value() ? column_width.value() : auto_width();

        auto cell = [&](auto && symbol, std::optional<std::size_t> symbol_bytes = std::nullopt)
        {
            // deal with unicode chars that mess up std::setw
            std::size_t bytes = symbol_bytes.value_or(
                is_traceback_matrix ? symbols.trace_dir_bytes : 1u);
            std::size_t length_bytes = char_length_bytes(symbol);
            std::size_t length = length_bytes / bytes;
            std::size_t offset = length_bytes - length;

            cout << std::left
                      << std::setw(_column_width + offset)
                      << symbol
                      << symbols.col_sep;
        };

        auto first_cell = [&](auto && symbol)
        {
            cout << symbol << symbols.col_sep;
        };

        // |_|d|a|t|a|b|a|s|e|
        auto first_row = [&]
        {
            first_cell(" ");
            cell(symbols.epsilon, symbols.epsilon_bytes);

            for(unsigned col = 0; col < matrix.cols()-1; ++col)
                cell(matrix.database()[col], 1);
            cout << "\n";
        };

        // |-|-|-|-|-|-|-|-|-|
        auto divider = [&]
        {
            cout << " " << symbols.row_col_sep;
            for(unsigned col = 0; col < matrix.cols(); ++col)
            {
                for(unsigned i = 0; i < _column_width; ++i)
                    cout << symbols.row_sep;
                cout << symbols.row_col_sep;
            }
            cout << "\n";
        };

        first_row();
        for(unsigned row = 0; row < matrix.rows(); ++row)
        {
            if (symbols.row_sep[0] != '\0')
                divider();

            // one query letter + one row of scores / traces
            if (row == 0)
                first_cell(symbols.epsilon);
            else
                first_cell(matrix.query()[row-1]);
            for(unsigned col = 0; col < matrix.cols(); ++col)
                cell(char_at(row, col));
            cout << "\n";
        }
    }

private:

    //!\brief Same as #matrix\.at(\a row, \a col), but converts the value to
    //!       a trace symbol (alignment_matrix_format::trace_dir) if the
    //!       #matrix is a traceback matrix.
    auto char_at(unsigned row, unsigned col) const noexcept
    {
        if constexpr(is_traceback_matrix)
        {
            trace_matrix_directions direction = matrix.at(row, col);
            return symbols.trace_dir[(unsigned)(direction) % 8u];
        } else
        {
            return matrix.at(row, col);
        }
    }

    //!\brief The length of the \a str (traceback symbols are unicode aware)
    std::size_t char_length(auto && str) const noexcept
    {
        if constexpr(is_traceback_matrix)
            return char_length_bytes(str) / symbols.trace_dir_bytes;
        return char_length_bytes(str);
    }

    //!\brief The number of bytes the \a str uses
    std::size_t char_length_bytes(auto && str) const noexcept
    {
        std::stringstream stream;
        stream << str;
        return stream.str().length();
    }
};

/*!\name Type deduction guides
 * \relates seqan3::alignment_matrix_formatter
 * \{
 */
template<typename alignment_matrix_t>
alignment_matrix_formatter(alignment_matrix_t const &) -> alignment_matrix_formatter<alignment_matrix_t>;

template<typename alignment_matrix_t>
alignment_matrix_formatter(alignment_matrix_t const &, alignment_matrix_format) -> alignment_matrix_formatter<alignment_matrix_t>;
//!\}
} // namespace seqan3
