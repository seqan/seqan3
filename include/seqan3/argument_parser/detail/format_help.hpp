// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Provides the format_help struct that print the help page to the command line
 *        and the two child formats (format_version, format_short_help) that print
 *        short help messages to the command line.
 */

#pragma once

#include <cassert>
#include <iostream>

#include <seqan3/argument_parser/detail/format_base.hpp>
#include <seqan3/argument_parser/detail/terminal.hpp>
#include <seqan3/core/detail/test_accessor.hpp>

namespace seqan3::detail
{

/*!\brief The format that prints the help page to std::cout.
 * \ingroup argument_parser
 *
 * \details
 * The help page printing is not done immediately, because the user might not
 * provide meta information, positional options, etc. in the correct order.
 * In addition the needed order would be different from the parse format.
 * Thus the calls are stored (parser_set_up_calls and positional_option_calls)
 * and only evaluated when calling format_help::parse().
 *
 * \remark For a complete overview, take a look at \ref argument_parser
 */
class format_help : public format_help_base<format_help>
{
    //!\brief The CRTP base class type.
    using base_type = format_help_base<format_help>;

    //!\brief Befriend the base class to give access to the private member functions.
    friend base_type;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    format_help() = default;                                //!< Defaulted.
    format_help(format_help const & pf) = default;          //!< Defaulted.
    format_help & operator=(format_help const &) = default; //!< Defaulted.
    format_help(format_help &&) = default;                  //!< Defaulted.
    format_help & operator=(format_help &&) = default;      //!< Defaulted.
    ~format_help() = default;                               //!< Defaulted.

    //!\copydoc format_help_base(std::vector<std::string> const &, bool const)
    format_help(std::vector<std::string> const & names, bool const advanced = false) : base_type{names, advanced} {};
    //!\}

protected:
    //!\privatesection
    //!\brief Stores the relevant parameters of the documentation on the screen.
    struct console_layout_struct
    {
        //!\brief The screen width.
        uint32_t screenWidth;
        //!\brief The default screen width.
        uint32_t defaultScreenWidth;
        //!\brief The maximal screen width.
        uint32_t maximalScreenWidth;
        //!\brief The minimal screen width.
        uint32_t minimalScreenWidth;
        //!\brief The left Padding.
        uint32_t leftPadding;
        //!\brief The center Padding.
        uint32_t centerPadding;
        //!\brief The right Padding.
        uint32_t rightPadding;
        //!\brief The left Column Width.
        uint32_t leftColumnWidth;
        //!\brief The right Column Width.
        uint32_t rightColumnWidth;
        //!\brief The right Column Tab.
        uint32_t rightColumnTab;

        //!\brief The constructor.
        //!\param[in] terminal_width The width of the terminal.
        console_layout_struct(uint32_t const terminal_width) :
            screenWidth{0},
            defaultScreenWidth{80},
            maximalScreenWidth{120},
            minimalScreenWidth{40},
            leftPadding{4},
            centerPadding{2},
            rightPadding{2},
            leftColumnWidth{4},
            rightColumnWidth{0}
        {
            // Guess terminal screen width and set into layout.
            screenWidth = (terminal_width > 0) ? terminal_width : defaultScreenWidth;
            screenWidth = std::max(screenWidth, minimalScreenWidth);
            screenWidth = std::min(screenWidth, maximalScreenWidth);
            screenWidth -= rightPadding;

            rightColumnWidth = screenWidth - leftPadding - leftColumnWidth - centerPadding - rightPadding;
            rightColumnTab = leftPadding + leftColumnWidth + centerPadding;
        }

        //!\brief The default constructor.
        console_layout_struct() : console_layout_struct{get_terminal_width()}
        {}
    };

    //!\brief Prints a help page header to std::cout.
    void print_header()
    {
        std::ostream_iterator<char> out(std::cout);

        std::cout << meta.app_name;
        if (!empty(meta.short_description))
            std::cout << " - " << meta.short_description;

        std::cout << "\n";
        unsigned len =
            text_width(meta.app_name) + (empty(meta.short_description) ? 0 : 3) + text_width(meta.short_description);
        std::fill_n(out, len, '=');
        std::cout << '\n';
    }

    /*!\brief Prints a help page section to std::cout.
     * \param[in] title The title of the subsection of the help page.
     */
    void print_section(std::string const & title)
    {
        std::ostream_iterator<char> out(std::cout);
        std::cout << '\n' << to_text("\\fB");
        std::transform(title.begin(),
                       title.end(),
                       out,
                       [](unsigned char c)
                       {
                           return std::toupper(c);
                       });
        std::cout << to_text("\\fP") << '\n';
        prev_was_paragraph = false;
    }

    /*!\brief Prints a help page subsection to std::cout.
     * \param[in] title The title of the subsection of the help page.
     */
    void print_subsection(std::string const & title)
    {
        std::ostream_iterator<char> out(std::cout);
        std::cout << '\n';
        std::fill_n(out, layout.leftPadding / 2, ' ');
        std::cout << in_bold(title) << '\n';
        prev_was_paragraph = false;
    }

    /*!\brief Prints a text to std::cout.
     * \param[in] text The text to print.
     * \param[in] line_is_paragraph Whether to insert as paragraph
     *            or just a line (only one line break if not a paragraph).
     */
    void print_line(std::string const & text, bool const line_is_paragraph)
    {
        if (prev_was_paragraph)
            std::cout << '\n';

        std::ostream_iterator<char> out(std::cout);
        std::fill_n(out, layout.leftPadding, ' ');
        print_text(text, layout.leftPadding);
        prev_was_paragraph = line_is_paragraph;
    }

    /*!\brief Prints a help page list_item to std::cout.
     * \param[in] term The key of the key-value pair of the list item.
     * \param[in] desc The value of the key-value pair of the list item.
     *
     * \details
     *
     * A list item is composed of a key (term) and value (desc)
     * and usually used for option identifier-description-pairs.
     * E.g.:
     *```console
     *     -a, --age LONG
     *            Super important integer for age.
     *```
     */
    void print_list_item(std::string const & term, std::string const & desc)
    {
        if (prev_was_paragraph)
            std::cout << '\n';

        std::ostream_iterator<char> out(std::cout);

        // Print term.
        std::fill_n(out, layout.leftPadding, ' ');
        std::cout << to_text(term);
        unsigned pos = layout.leftPadding + term.size();
        if (pos + layout.centerPadding > layout.rightColumnTab)
        {
            std::cout << '\n';
            pos = 0;
        }
        std::fill_n(out, layout.rightColumnTab - pos, ' ');
        print_text(desc, layout.rightColumnTab);

        prev_was_paragraph = false;
    }

    //!\brief Prints a help page footer to std::cout.
    void print_footer()
    {
        // no footer
    }

    /*!\brief Formats text for pretty command line printing.
     * \param[in] str The input string to format for correct command line printing.
     */
    std::string to_text(std::string const & str)
    {
        std::string result;

        for (auto it = str.begin(); it != str.end(); ++it)
        {
            if (*it == '\\')
            {
                // Handle escape sequence, we interpret only "\-", "\fI", and "\fB".
                ++it;
                assert(it != str.end());
                if (*it == '-')
                {
                    result.push_back(*it);
                }
                else if (*it == 'f')
                {
                    ++it;
                    assert(it != str.end());
                    if (*it == 'I')
                    {
                        if (is_terminal())
                            result.append("\033[4m");
                    }
                    else if (*it == 'B')
                    {
                        if (is_terminal())
                            result.append("\033[1m");
                    }
                    else if (*it == 'P')
                    {
                        if (is_terminal())
                            result.append("\033[0m");
                    }
                    else
                    {
                        result.append("\\f");
                        result.push_back(*it);
                    }
                }
                else
                {
                    result.push_back('\\');
                    result.push_back(*it);
                }
            }
            else
            {
                result.push_back(*it);
            }
        }

        return result;
    }

    /*!\brief Returns width of text if printed.
     * \param[in] text The string to compute the width for on the command line.
     * /detail Note: "\-" has length 1, "\fI", "\fB", "\fP" have length 0.
     */
    unsigned text_width(std::string const & text)
    {
        unsigned result = 0;

        for (unsigned i = 0; i < text.size(); ++i)
        {
            if (text[i] != '\\')
            {
                result += 1;
                continue;
            }

            if (i + 1 == text.size())
            {
                result += 1; // Will print "\\".
                continue;
            }

            if (text[i + 1] == '\\' || text[i + 1] == '-')
            {
                i += 1;
                result += 1;
                continue; // Will print '\\' or '-'.
            }

            if (i + 2 == text.size())
            {
                i += 1;
                result += 2; // Will print two chars.
                continue;
            }

            if (text[i + 1] == 'f')
            {
                if (text[i + 2] == 'B' || text[i + 2] == 'I' || text[i + 2] == 'P')
                    i += 2; // Skip f and {B, I, P}.
                else
                    result += 1;
            }
        }

        return result;
    }

    /*!\brief Prints text with correct line wrapping to the command line (std::cout).
     * \param[in] text   The string to print on the command line.
     * \param[in] tab    The position offset (indentation) to start printing at.
     */
    void print_text(std::string const & text, unsigned const tab)
    {
        unsigned pos = tab;
        std::ostream_iterator<char> out(std::cout);

        // Tokenize the text.
        std::istringstream iss(text.c_str());
        std::vector<std::string> tokens;
        std::ranges::copy(std::istream_iterator<std::string>(iss),
                          std::istream_iterator<std::string>(),
                          std::back_inserter(tokens));

        // Print the text.
        assert(pos <= tab);
        std::fill_n(out, tab - pos, ' '); // go to tab

        pos = tab;
        typedef std::vector<std::string>::const_iterator TConstIter;
        for (TConstIter it = tokens.begin(); it != tokens.end(); ++it)
        {
            if (it == tokens.begin())
            {
                std::cout << to_text(*it);
                pos += text_width(*it);
                if (pos > layout.screenWidth)
                {
                    std::cout << '\n';
                    std::fill_n(out, tab, ' ');
                    pos = tab;
                }
            }
            else
            {
                if (pos + 1 + text_width(*it) > layout.screenWidth)
                {
                    // Would go over screen with next, print current word on next line.
                    std::cout << '\n';
                    fill_n(out, tab, ' ');
                    std::cout << to_text(*it);
                    pos = tab + text_width(*it);
                }
                else
                {
                    std::cout << ' ';
                    std::cout << to_text(*it);
                    pos += text_width(*it) + 1;
                }
            }
        }
        if (!empty(tokens))
            std::cout << '\n';
    }

    /*!\brief Format string in bold.
     * \param[in] str The input string to format in bold.
     * \returns The string `str` wrapped in bold formatting.
     */
    std::string in_bold(std::string const & str)
    {
        return to_text("\\fB") + str + to_text("\\fP");
    }

    //!\brief Needed for correct formatting while calling different print functions.
    bool prev_was_paragraph{false};

    //!\brief Befriend seqan3::detail::test_accessor to grant access to layout.
    friend struct ::seqan3::detail::test_accessor;

    //!\brief Stores the relevant parameters of the documentation on the screen.
    console_layout_struct layout{};
};

/*!\brief The format that prints a short help message to std::cout.
 * \ingroup argument_parser
 *
 * \details
 *
 * The short help message printing is not done immediately, because the user cannot provide
 * meta information (e.g. app_name) on construction of the parser. Thus the meta information is collected
 * and only evaluated when calling seqan3::detail::format_version::parse.
 *
 * \remark For a complete overview, take a look at \ref argument_parser
 */
class format_short_help : public format_help
{
public:
    /*!\brief Initiates the printing of a short help message to std::cout.
     * \param[in] parser_meta The meta information that are needed for a detailed version information.
     */
    void parse(argument_parser_meta_data const & parser_meta)
    {
        meta = parser_meta;

        print_header();

        if (!parser_meta.synopsis.empty())
            print_synopsis();

        print_line("Try -h or --help for more information.\n", true);

        std::exit(EXIT_SUCCESS);
    }
};

/*!\brief The format that prints the version to std::cout.
 * \ingroup argument_parser
 *
 * \details
 *
 * The version printing is not done immediately, because the user cannot provide
 * meta information on construction of the parser. Thus the meta information is collected
 * and only evaluated when calling format_version::parse().
 *
 * \remark For a complete overview, take a look at \ref argument_parser
 */
class format_version : public format_help
{
public:
    /*!\brief Initiates the printing of the version information to std::cout.
     * \param[in] parser_meta The meta information that are needed for a detailed version information.
     */
    void parse(argument_parser_meta_data & parser_meta)
    {
        meta = parser_meta;

        print_header();
        print_version();

        std::exit(EXIT_SUCCESS); // program should not continue from here
    }
};

/*!\brief The format that prints the copyright information to std::cout.
 * \ingroup argument_parser
 *
 * \details
 *
 * The copyright message printing is not done immediately, because the user cannot provide
 * meta information (e.g. long_copyright) on construction of the parser. Thus the meta information is collected
 * and only evaluated when calling seqan3::detail::format_version::parse.
 *
 * \remark For a complete overview, take a look at \ref argument_parser
 */
class format_copyright : public format_help
{
public:
    /*!\brief Initiates the printing of the copyright message to std::cout.
     * \param[in] parser_meta The meta information that are needed for a detailed version information.
     */
    void parse(argument_parser_meta_data const & parser_meta)
    {
        meta = parser_meta;
        debug_stream_type stream{std::cout};
        std::string seqan_license{
            R"(Copyright (c) 2006-2025, Knut Reinert & Freie Universität Berlin
Copyright (c) 2016-2025, Knut Reinert & MPI für molekulare Genetik
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Knut Reinert or the FU Berlin nor the names of
      its contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.)"};

        stream << std::string(80, '=') << "\n"
               << in_bold("Copyright information for " + meta.app_name + ":\n") << std::string(80, '-') << '\n';

        if (!empty(meta.long_copyright))
        {
            stream << to_text("\\fP") << meta.long_copyright << "\n";
        }
        else if (!empty(meta.short_copyright))
        {
            stream << in_bold(meta.app_name + " full copyright information not available. "
                              + "Displaying short copyright information instead:\n")
                   << meta.short_copyright << "\n";
        }
        else
        {
            stream << to_text("\\fP") << meta.app_name << " copyright information not available.\n";
        }

        stream << std::string(80, '=') << '\n'
               << in_bold("This program contains SeqAn code licensed under the following terms:\n")
               << std::string(80, '-') << '\n'
               << seqan_license << '\n';

        std::exit(EXIT_SUCCESS);
    }
};

} // namespace seqan3::detail
