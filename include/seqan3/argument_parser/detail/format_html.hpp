// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Provides the format_html struct and its helper functions.
 */

#pragma once

#include <iostream>

#include <seqan3/argument_parser/detail/format_base.hpp>
#include <seqan3/argument_parser/detail/terminal.hpp>
#include <seqan3/version.hpp>

namespace seqan3::detail
{

/*!\brief The format that prints the help page as html to std::cout.
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
class format_html : public format_help_base<format_html>
{
    //!\brief The CRTP base class type.
    using base_type = format_help_base<format_html>;

    //!\brief Befriend the base class to give access to the private member functions.
    friend base_type;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    format_html() = default;                                   //!< Defaulted.
    format_html(format_html const & pf) = default;             //!< Defaulted.
    format_html & operator=(format_html const & pf) = default; //!< Defaulted.
    format_html(format_html &&) = default;                     //!< Defaulted.
    format_html & operator=(format_html &&) = default;         //!< Defaulted.
    ~format_html() = default;                                  //!< Defaulted.

    //!\copydoc format_help_base(std::vector<std::string> const &, bool const)
    format_html(std::vector<std::string> const & names, bool const advanced = false) : base_type{names, advanced} {};
    //!\}

private:
    //!\brief Closes HTML list tag (dl) if needed.
    void maybe_close_list()
    {
        if (is_dl)
        {
            std::cout << "</dl>\n";
            is_dl = false;
        }
    }

    //!\brief Closes HTML paragraph tag (p) if needed
    void maybe_close_paragraph()
    {
        if (is_p)
        {
            std::cout << "</p>\n";
            is_p = false;
        }
    }

    //!\brief Prints a help page header in HTML format to std::cout.
    void print_header()
    {
        // Print HTML boilerplate header.
        std::cout << "<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01//EN\" "
                  << "http://www.w3.org/TR/html4/strict.dtd\">\n"
                  << "<html lang=\"en\">\n"
                  << "<head>\n"
                  << "<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\">\n"
                  << "<title>" << escape_special_xml_chars(meta.app_name) << " &mdash; "
                  << escape_special_xml_chars(meta.short_description) << "</title>\n"
                  << "</head>\n"
                  << "<body>\n";

        std::cout << "<h1>" << to_html(meta.app_name) << "</h1>\n"
                  << "<div>" << to_html(meta.short_description) << "</div>\n";
    }

    /*!\brief Prints a section title in HTML format to std::cout.
     * \param[in] title The title of the section of the help page.
     */
    void print_section(std::string const & title)
    {
        // SEQAN_ASSERT_NOT_MSG(isDl && isP, "Current <dl> and <p> are mutually exclusive.");
        maybe_close_list();
        maybe_close_paragraph();
        std::cout << "<h2>" << to_html(title) << "</h2>\n";
    }

    /*!\brief Prints a subsection title in HTML format to std::cout.
     * \param[in] title The title of the subsection of the help page.
     */
    void print_subsection(std::string const & title)
    {
        // SEQAN_ASSERT_NOT_MSG(isDl && isP, "Current <dl> and <p> are mutually exclusive.");
        maybe_close_list();
        maybe_close_paragraph();
        std::cout << "<h3>" << to_html(title) << "</h3>\n";
    }

    /*!\brief Prints a text in HTML format to std::cout.
     * \param[in] text The text to print.
     * \param[in] line_is_paragraph Whether to insert as paragraph
     *            or just a line (only one line break if not a paragraph).
     */
    void print_line(std::string const & text, bool line_is_paragraph)
    {
        // SEQAN_ASSERT_NOT_MSG(isDl && isP, "Current <dl> and <p> are mutually exclusive.");
        maybe_close_list();
        if (!is_p) // open parapgraph
        {
            std::cout << "<p>\n";
            is_p = true;
        }
        std::cout << to_html(text) << "\n";
        if (line_is_paragraph)
            maybe_close_paragraph();
        else
            std::cout << "<br>\n";
    }

    /*!\brief Prints a help page list_item in HTML format to std::cout.
     * \param[in] term The key of the key-value pair of the list item.
     * \param[in] desc The value of the key-value pair of the list item.
     *
     * \details
     *
     * A list item is composed of a key (term) and value (desc)
     * and usually used for option identifier-description-pairs.
     */
    void print_list_item(std::string const & term, std::string const & desc)
    {
        // SEQAN_ASSERT_NOT_MSG(isDl && isP, "Current <dl> and <p> are mutually exclusive.");
        maybe_close_paragraph();

        if (!is_dl)
        {
            std::cout << "<dl>\n";
            is_dl = true;
        }
        std::cout << "<dt>" << to_html(term) << "</dt>\n"
                  << "<dd>" << to_html(desc) << "</dd>\n";
    }

    //!\brief Prints a help page footer in HTML format to std::cout.
    void print_footer()
    {
        maybe_close_paragraph();

        // Print HTML boilerplate footer.
        std::cout << "</body></html>";
    }

    /*!\brief Converts console output formatting to the HTML equivalent.
     *
     * \param[in] input The text to transform.
     */
    std::string to_html(std::string const & input)
    {
        std::string buffer = escape_special_xml_chars(input);
        std::string result;
        std::vector<std::string> open_tags; // acts as a stack of html tags

        for (auto it = input.begin(); it != input.end(); ++it)
        {
            if (*it == '\\')
            {
                // Handle escape sequence, we interpret only "\-", "\fI", and "\fB".
                ++it;
                assert(!(it == input.end()));
                if (*it == '-')
                {
                    result.push_back(*it);
                }
                else if (*it == 'f')
                {
                    ++it;
                    assert(!(it == input.end()));
                    if (*it == 'I')
                    {
                        open_tags.push_back("em");
                        result.append("<em>");
                    }
                    else if (*it == 'B')
                    {
                        open_tags.push_back("strong");
                        result.append("<strong>");
                    }
                    else if (*it == 'P')
                    {
                        assert(!open_tags.empty());
                        result.append("</");
                        result.append(open_tags.back());
                        result.append(">");
                        open_tags.pop_back();
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

    /*!\brief Format string as in_bold.
     * \param[in] str The input string to format in bold.
     * \returns The string `str` wrapped in bold formatting.
     */
    std::string in_bold(std::string const & str)
    {
        return "<strong>" + str + "</strong>";
    }

    //!\brief Current state is either inside a html \<dl\> tag (true) or not (false).
    bool is_dl{false};
    //!\brief Current state is either inside a html \<p\> tag (true) or not (false).
    bool is_p{false};
};

} // namespace seqan3::detail
