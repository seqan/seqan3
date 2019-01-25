// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Contains the format_html struct and its helper functions.
 */

#pragma once

#include <iostream>

#include <seqan3/argument_parser/detail/format_base.hpp>
#include <seqan3/core/detail/terminal.hpp>
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
 */
class format_html : public format_base
{
public:
    /*!\brief Adds a print_list_item call to be evaluated later on.
     *
     * \tparam option_type    The type of variable in which to store the given command line argument.
     * \tparam validator_type The type of validator applied to the value after parsing.
     *
     * \param[out] value     The variable in which to store the given command line argument.
     * \param[in]  short_id  The short identifier for the option (e.g. 'i').
     * \param[in]  long_id   The long identifier for the option (e.g. "integer").
     * \param[in]  desc      The description of the option.
     * \param[in]  spec      Advanced option specification. see seqan3::option_spec.
     * \param[in]  validator The validator applied to the value after parsing (callable).
     */
    template <typename option_type, typename validator_type>
    void add_option(option_type & value,
                    char const short_id,
                    std::string const & long_id,
                    std::string const & desc,
                    option_spec const & spec,
                    validator_type && validator)
    {
        parser_set_up_calls.push_back([=, &value]()
        {
            if (!(spec & option_spec::HIDDEN) && (!(spec & option_spec::ADVANCED)))
                print_list_item(prep_id_for_help(short_id, long_id) + " " + option_type_and_list_info(value),
                                 (desc + " " + validator.get_help_page_message()));
        });
    }

    /*!\brief Adds a print_list_item call to be evaluated later on.
     *
     * \param[in]  short_id The short identifier for the flag (e.g. 'i').
     * \param[in]  long_id  The long identifier for the flag (e.g. "integer").
     * \param[in]  desc     The description of the flag.
     * \param[in]  spec     Advanced flag specification. see seqan3::option_spec.
     */
    void add_flag(bool & /*value*/,
                  char const short_id,
                  std::string const & long_id,
                  std::string const & desc,
                  option_spec const & spec)
    {
        parser_set_up_calls.push_back([=] ()
        {
            if (!(spec & option_spec::HIDDEN) && (!(spec & option_spec::ADVANCED)))
                print_list_item(prep_id_for_help(short_id, long_id), desc);
        });
    }

    /*!\brief Adds a print_list_item call to be evaluated later on.
     *
     * \tparam option_type    The type of variable in which to store the given command line argument.
     * \tparam validator_type The type of validator applied to the value after parsing.
     *
     * \param[out] value     The variable in which to store the given command line argument.
     * \param[in]  desc      The description of the positional option.
     * \param[in]  validator The validator applied to the value after parsing (callable).
     */
    template <typename option_type, typename validator_type>
    void add_positional_option(option_type & value,
                               std::string const & desc,
                               validator_type && validator)
    {
        ++positional_option_count;

        positional_option_calls.push_back([=, &value]()
        {
            std::string key{"\\fBARGUMENT " + std::to_string(positional_option_count) + "\\fP " + option_type_and_list_info(value)};
            print_list_item(key, (desc + " " + validator.get_help_page_message()));
        });
    }

    /*!\brief Initiates the printing of the help page to std::cout.
     * \param[in] parser_meta The meta information that are needed for a detailed help page.
     */
    void parse(argument_parser_meta_data const & parser_meta)
    {
        meta = parser_meta;

        print_header();

        print_section("Synopsis");
        _print_synopsis();

        if (!meta.description.empty())
        {
            print_section("Description");
            for (auto desc : meta.description)
                print_line(desc);
        }

        // add positional options if specified
        if (!positional_option_calls.empty())
            print_section("Positional Arguments");

        // each call will evaluate the function print_list_item()
        for (auto f : positional_option_calls)
            f();

        // add options and flags if specified
        if (!parser_set_up_calls.empty())
            print_section("Options");

        // each call will evaluate the function print_list_item()
        for (auto f : parser_set_up_calls)
            f();

        print_footer();

        throw parser_interruption();
    }

    /*!\brief Adds a print_section call to parser_set_up_calls.
     * \param[in] title The section title to print.
     */
    void add_section(std::string const & title)
    {
        parser_set_up_calls.push_back([=] ()
        {
            print_section(title);
        });
    }
    /*!\brief Adds a print_subsection call to parser_set_up_calls.
     * \param[in] title The subsection title to print.
     */
    void add_subsection(std::string const & title)
    {
        parser_set_up_calls.push_back([=] ()
        {
            print_subsection(title);
        });
    }

    /*!\brief Adds a print_line call to parser_set_up_calls.
     * \param[in] text              The text to print.
     * \param[in] line_is_paragraph Specify whether the text should be ended with
     *                              a new line or not.
     */
    void add_line(std::string const & text, bool const line_is_paragraph)
    {
        parser_set_up_calls.push_back([=] ()
        {
            print_line(text, line_is_paragraph);
        });
    }

    /*!\brief Adds a print_list_item call to parser_set_up_calls.
     * \param[in] key         The name that will be listed on the left hand side.
     * \param[in] description The description that will be listed on the right
     *                        hand side.
     */
    void add_list_item(std::string const & key, std::string const & description)
    {
        parser_set_up_calls.push_back([=] ()
        {
            print_list_item(key, description);
        });
    }

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

    //!\brief Prints a help page synopsis in HTML format section to std::cout.
    void _print_synopsis()
    {
        for (unsigned i = 0; i < meta.synopsis.size(); ++i)
        {
            std::string text = "\\fB";
            text.append(meta.app_name);
            text.append("\\fP ");
            text.append(meta.synopsis[i]);

            print_line(text, false);
        }
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
            std::cout << "<br />\n";
    }

    /*!\brief Prints a text in HTML format to std::cout.
     * \param[in] text The text to print.
     */
    void print_line(std::string const & text)
    {
        print_line(text, true);
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
        maybe_close_list();

        // Print version, date and url.
        std::cout << "<h2>Version</h2>\n"
                  << "<strong>Last update:</strong> " << to_html(meta.date) << "<br>\n<strong>"
                  << meta.app_name << " version:</strong> " << meta.version << "<br>\n"
                  << "<strong>SeqAn version:</strong> " << SEQAN3_VERSION_MAJOR << '.' <<  SEQAN3_VERSION_MINOR << '.'
                  << SEQAN3_VERSION_PATCH << "<br>\n";

        if (!meta.url.empty())
        {
            std::cout << "<h2>Url</h2>\n"
                    << meta.url << "<br>\n";
        }
        std::cout << "<br>\n";

        // Print legal stuff
        if ((!meta.short_copyright.empty()) || (!meta.long_copyright.empty()) || (!meta.citation.empty()))
        {
            std::cout << "<h2>Legal</h2>\n<strong>";

            if (!meta.short_copyright.empty())
                std::cout << meta.app_name << " Copyright: </strong>"
                       << meta.short_copyright << "<br>\n<strong>";

            std::cout << "SeqAn Copyright:</strong> 2006-2015 Knut Reinert, FU-Berlin; released under the 3-clause BSDL.<br>\n<strong>";

            if (!meta.citation.empty())
                std::cout << "In your academic works please cite:</strong> " << meta.citation << "<br>\n";
            else
                std::cout << "</strong>";

            if (!meta.long_copyright.empty())
                std::cout << "For full copyright and/or warranty information see <tt>--copyright</tt>.\n";
        }

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

    /*!\brief Stores all meta information about the application.
     *
     * \details
     *
     * This needs to be a member of format_parse, because it needs to present
     * (not filled) when the parser_set_up_calls vector is filled, since all
     * printing functions need some meta information.
     * The member variable itself is filled when copied over from the argument_parser
     * when calling format_parse::parse. That way all the information needed are
     * there, when the actual printing starts.
     *
     * This is not private because it is needed for short but nicely
     * formatted (error) output to the command line.
     */
    argument_parser_meta_data meta;
    //!\brief Vector of functions that stores all calls except add_positional_option.
    std::vector<std::function<void()>> parser_set_up_calls;
    //!\brief Vector of functions that stores add_positional_option calls.
    std::vector<std::function<void()>> positional_option_calls; // singled out to be printed on top
    //!\brief Current state is either inside a html \<dl\> tag (true) or not (false).
    bool is_dl{false};
    //!\brief Current state is either inside a html \<p\> tag (true) or not (false).
    bool is_p{false};
    //!\brief Keeps track of the number of positional options.
    unsigned positional_option_count{0};
};

} // namespace seqan3
