// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Contains the format_help struct that print the help page to the command line
 *        and the two child formats (format_version, format_short_help) that print
 *        short help messages to the command line.
 */

#pragma once

#include <iostream>

#include <seqan3/argument_parser/detail/format_base.hpp>
#include <seqan3/core/detail/terminal.hpp>
#include <seqan3/version.hpp>

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
 */
class format_help : public format_base
{
public:
    /*!\name Constructors, destructor and assignment
     * \brief The (default) constructors.
     * \{
     */
    format_help() = default;
    format_help(format_help const & pf) = default;
    format_help & operator=(format_help const & pf) = default;
    format_help(format_help &&) = default;
    format_help & operator=(format_help &&) = default;

    /*!\brief Initializes a format_help object.
     *
     * \param advanced Set to `true` to show advanced options.
     */
    format_help(bool advanced) :
        show_advanced_options(advanced)
    {}

    //!\brief The destructor.
    ~format_help() = default;
    //!\}

    /*!\brief Adds a seqan3::print_list_item call to be evaluated later on.
     *
     * \tparam option_type The type of variable in which to store the given command line argument.
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
            if (!(spec & option_spec::HIDDEN) && (!(spec & option_spec::ADVANCED) || show_advanced_options))
                print_list_item(prep_id_for_help(short_id, long_id) + " " + option_type_and_list_info(value),
                                 (desc + " " + validator.get_help_page_message()));
        });
    }

    /*!\brief Adds a seqan3::print_list_item call to be evaluated later on.
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
            if (!(spec & option_spec::HIDDEN) && (!(spec & option_spec::ADVANCED) || show_advanced_options))
                print_list_item(prep_id_for_help(short_id, long_id), desc);
        });
    }

    /*!\brief Adds a seqan3::print_list_item call to be evaluated later on.
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
                               validator_type & validator)
    {
        positional_option_calls.push_back([=, &value]()
        {
            ++positional_option_count;
            std::string key{"\\fBARGUMENT-" + std::to_string(positional_option_count) + "\\fP " + option_type_and_list_info(value)};
            print_list_item(key, (desc + " " + validator.get_help_page_message()));
        });
    }

    /*!\brief Initiates the printing of the help page to std::cout.
     * \param[in] parser_meta The meta information that are needed for a detailed help page.
     */
    void parse(argument_parser_meta_data & parser_meta)
    {
        meta = parser_meta;

        print_header();

        if (!meta.synopsis.empty())
        {
            print_section("Synopsis");
            print_synopsis();
        }

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

        if (!meta.examples.empty())
        {
            print_section("Examples");
            for (auto example : meta.examples)
                print_line(example);
        }

        print_footer();

        throw parser_interruption(); // program should not continue from here
    }

    /*!\brief Adds a print_section call to parser_set_up_calls.
     * \param[in] title The title of the section of the help page.
     */
    void add_section(std::string const & title)
    {
        parser_set_up_calls.push_back([=] ()
        {
            print_section(title);
        });
    }

    /*!\brief Adds a print_subsection call to parser_set_up_calls.
     * \param[in] title The title of the subsection of the help page.
     */
    void add_subsection(std::string const & title)
    {
        parser_set_up_calls.push_back([=] ()
        {
            print_subsection(title);
        });
    }

    /*!\brief Adds a print_line call to parser_set_up_calls.
     * \param[in] text The line text to be printed to the help page.
     * \param[in] line_is_paragraph True if you want a new line at the end.
     */
    void add_line(std::string const & text, bool line_is_paragraph)
    {
        parser_set_up_calls.push_back([=] ()
        {
            print_line(text, line_is_paragraph);
        });
    }

    /*!\brief Adds a seqan3::print_list_item call to parser_set_up_calls.
     * \param[in] key The key of the key-value pair list item.
     * \param[in] desc The key of the key-value pair list item.
     */
    void add_list_item(std::string const & key, std::string const & desc)
    {
        parser_set_up_calls.push_back([=] ()
        {
            print_list_item(key, desc);
        });
    }

    /*!\brief Stores all meta information about the application
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
     * This function is not private because it is needed for short but nicely
     * formatted (error) output to the command line.
     */
    argument_parser_meta_data meta;

protected:
    //!\privatesection
    //!\brief Stores the relevant parameters of the documentation on the screen.
    struct console_layout_struct
    {
        //!\brief The screen width.
        unsigned screenWidth;
        //!\brief The default screen width.
        unsigned defaultScreenWidth;
        //!\brief The maximal screen width.
        unsigned maximalScreenWidth;
        //!\brief The minimal screen width.
        unsigned minimalScreenWidth;
        //!\brief The left Padding.
        unsigned leftPadding;
        //!\brief The center Padding.
        unsigned centerPadding;
        //!\brief The right Padding.
        unsigned rightPadding;
        //!\brief The left Column Width.
        unsigned leftColumnWidth;
        //!\brief The right Column Width.
        unsigned rightColumnWidth;
        //!\brief The right Column Tab.
        unsigned rightColumnTab;

        //!\brief The constructor.
        console_layout_struct() :
            screenWidth{0}, defaultScreenWidth{80}, maximalScreenWidth{120}, minimalScreenWidth{40},
            leftPadding{4}, centerPadding{2}, rightPadding{2}, leftColumnWidth{4}, rightColumnWidth{0}
        {
            // Guess terminal screen width and set into layout.
            unsigned cols = get_terminal_width();
            screenWidth = (cols > 0) ? cols : defaultScreenWidth;
            screenWidth = std::max(screenWidth, minimalScreenWidth);
            screenWidth = std::min(screenWidth, maximalScreenWidth);
            screenWidth -= rightPadding;

            rightColumnWidth = screenWidth - leftPadding - leftColumnWidth - centerPadding - rightPadding;
            rightColumnTab = leftPadding + leftColumnWidth + centerPadding;
        }
    };

    //!\brief Prints a help page header to std::cout.
    void print_header()
    {
        std::ostream_iterator<char> out(std::cout);

        std::cout << meta.app_name;
        if (!empty(meta.short_description))
            std::cout << " - " << meta.short_description;

        std::cout << "\n";
        unsigned len = text_width(meta.app_name) + (empty(meta.short_description) ? 0 : 3) +
                       text_width(meta.short_description);
        std::fill_n(out, len, '=');
        std::cout << '\n';
    }

    //!\brief Prints a help page synopsis section to std::cout.
    void print_synopsis()
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

    /*!\brief Prints a help page section to std::cout.
     * \param[in] title The title of the subsection of the help page.
     */
    void print_section(std::string const & title)
    {
        std::ostream_iterator<char> out(std::cout);
        std::cout << '\n' << to_text("\\fB");
        std::transform(title.begin(), title.end(), out, [] (unsigned char c) { return std::toupper(c); });
        std::cout << to_text("\\fP") << '\n';
        prev_was_paragraph = false;
    }

    /*!\brief Prints a help page subsection to std::cout.
     * \param[in] title The title of the subsection of the help page.
     */
    void print_subsection(std::string const & title)
    {
        std::ostream_iterator<char> out(std::cout);
        std::cout << '\n' << to_text("\\fB");
        std::fill_n(out, layout.leftPadding / 2, ' ');
        std::cout << title << to_text("\\fP") << '\n';
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
        print_text(text, layout, layout.leftPadding);
        prev_was_paragraph = line_is_paragraph;
    }

    /*!\brief Prints a text to std::cout.
     * \param[in] text The line of text to print to the help page.
     */
    void print_line(std::string const & text)
    {
        print_line(text, true);
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
     *-------------------------------------------
     *     -a, --age LONG
     *            Super important integer for age.
     *-------------------------------------------
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
        print_text(desc, layout, layout.rightColumnTab);

        prev_was_paragraph = false;
    }

    //!\brief Prints the version information to std::cout.
    void print_version()
    {
        std::ostream_iterator<char> out(std::cout);

        // Print version, date and url.
        std::cout << "\n" << to_text("\\fB") << "VERSION" << to_text("\\fP") << "\n";
        std::fill_n(out, layout.leftPadding, ' ');
        std::cout << to_text("\\fB") << "Last update: " << to_text("\\fP") << meta.date << "\n";
        std::fill_n(out, layout.leftPadding, ' ');
        std::cout << to_text("\\fB") << meta.app_name << " version: " << to_text("\\fP") << meta.version << "\n";
        std::fill_n(out, layout.leftPadding, ' ');
        std::cout << to_text("\\fB") << "SeqAn version: " << to_text("\\fP") << SEQAN3_VERSION_MAJOR << '.'
                  <<  SEQAN3_VERSION_MINOR << '.' << SEQAN3_VERSION_PATCH;

        if (!empty(meta.url))
        {
            std::cout <<  "\n" << to_text("\\fB") << "URL" << to_text("\\fP") << "\n";
            std::fill_n(out, layout.leftPadding, ' ');
            std::cout << meta.url << "\n";
        }
        std::cout << "\n";
    }

    //!\brief Prints a help page footer to std::cout.
    void print_footer()
    {
        print_version();

        std::ostream_iterator<char> out(std::cout);

        // Print legal stuff
        if ((!empty(meta.short_copyright)) || (!empty(meta.long_copyright)) || (!empty(meta.citation)))
        {
            std::cout << "\n" << to_text("\\fB") << "LEGAL" << to_text("\\fP") << "\n";

            if (!empty(meta.short_copyright))
            {
                std::fill_n(out, layout.leftPadding, ' ');
                std::cout << to_text("\\fB") << meta.app_name << " Copyright: "
                          << to_text("\\fP") << meta.short_copyright << "\n";
            }
            std::fill_n(out, layout.leftPadding, ' ');
            std::cout << to_text("\\fB") << "SeqAn Copyright: " << to_text("\\fP")
                      << "2006-2015 Knut Reinert, FU-Berlin; released under the 3-clause BSDL.\n";
            if (!empty(meta.citation))
            {
                std::fill_n(out, layout.leftPadding, ' ');
                std::cout << to_text("\\fB") << "In your academic works please cite: " << to_text("\\fP")
                          << meta.citation << "\n";
            }
            if (!empty(meta.long_copyright))
            {
                std::fill_n(out, layout.leftPadding, ' ');
                std::cout << "For full copyright and/or warranty information see " << to_text("\\fB")
                          << "--copyright" << to_text("\\fP") << ".\n";
            }
        }
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
                result += 1;  // Will print "\\".
                continue;
            }

            if (text[i + 1] == '\\' || text[i + 1] == '-')
            {
                i += 1;
                result += 1;
                continue;  // Will print '\\' or '-'.
            }

            if (i + 2 == text.size())
            {
                i += 1;
                result += 2;  // Will print two chars.
                continue;
            }

            if (text[i + 1] == 'f')
            {
                if (text[i + 2] == 'B' || text[i + 2] == 'I' || text[i + 2] == 'P')
                    i += 2;  // Skip f and {B, I, P}.
                else
                    result += 1;
            }
        }

        return result;
    }

    /*!\brief Prints text with correct line wrapping to the command line (std::cout).
     * \param[in] text   The string to print on the command line.
     * \param[in] layout The command line parameters (e.g. padding) for correct printing.
     * \param[in] tab    The position offset (indentation) to start printing at.
     */
    void print_text(std::string const & text,
                    console_layout_struct const & layout,
                    unsigned const tab)
    {
        unsigned pos = tab;
        std::ostream_iterator<char> out(std::cout);

        // Tokenize the text.
        std::istringstream iss(text.c_str());
        std::vector<std::string> tokens;
        std::copy(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>(),
                  std::back_inserter<std::vector<std::string> >(tokens));

        // Print the text.
        assert(pos <= tab);
        std::fill_n(out, tab - pos, ' ');  // go to tab

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

    //!\brief Vector of functions that stores all calls except add_positional_option.
    std::vector<std::function<void()>> parser_set_up_calls;
    //!\brief Vector of functions that stores add_positional_option calls.
    std::vector<std::function<void()>> positional_option_calls; // singled out to be printed on top
    //!\brief Keeps track of the number of positional options
    unsigned positional_option_count{0};
    //!\brief Needed for correct formatting while calling different print functions.
    bool prev_was_paragraph{false};
    //!\brief Whether to show advanced options or not.
    bool show_advanced_options{false};
    //!\brief Stores the relevant parameters of the documentation on the screen.
    console_layout_struct layout;
};

/*!\brief The format that prints a short help message to std::cout.
 * \ingroup argument_parser
 *
 * \details
 *
 * The short help message printing is not done immediately, because the user cannot provide
 * meta information (e.g. app_name) on construction of the parser. Thus the meta information is collected
 * and only evaluated when calling seqan3::format_version::parse.
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

        print_line("Try -h or --help for more information.\n");

        throw parser_interruption();
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
 */
class format_version : public format_help
{
public:
     /*!\brief Adds a seqan3::print_list_item call to be evaluated later on.
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
            if (!(spec & option_spec::HIDDEN) && (!(spec & option_spec::ADVANCED) || show_advanced_options))
                print_list_item(prep_id_for_help(short_id, long_id) + " " + option_type_and_list_info(value),
                                 (desc + " " + validator.get_help_page_message()));
        });
    }

    /*!\brief Adds a seqan3::print_list_item call to be evaluated later on.
     *
     * \param[out] value    The variable in which to store the given command line argument
     *                      (which is only needed for the interface here).
     * \param[in]  short_id The short identifier for the flag (e.g. 'i').
     * \param[in]  long_id  The long identifier for the flag (e.g. "integer").
     * \param[in]  desc     The description of the flag.
     * \param[in]  spec     Advanced flag specification. see seqan3::option_spec.
     */
    void add_flag(bool & value,
                  char const short_id,
                  std::string const & long_id,
                  std::string const & desc,
                  option_spec const & spec)
    {
        parser_set_up_calls.push_back([=, &value] ()
        {
            if (!(spec & option_spec::HIDDEN) && (!(spec & option_spec::ADVANCED) || show_advanced_options))
                print_list_item(prep_id_for_help(short_id, long_id), desc);
        });
    }

    /*!\brief Adds a seqan3::print_list_item call to be evaluated later on.
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
        positional_option_calls.push_back([=, &value]()
        {
            ++positional_option_count;
            std::string key{"\\fBARGUMENT " + std::to_string(positional_option_count) + "\\fP " + option_type_and_list_info(value)};
            print_list_item(key, (desc + " " + validator.get_help_page_message()));
        });
    }

    /*!\brief Initiates the printing of the version information to std::cout.
     * \param[in] parser_meta The meta information that are needed for a detailed version information.
     */
    void parse(argument_parser_meta_data & parser_meta)
    {
        meta = parser_meta;

        print_header();
        print_version();

        throw parser_interruption(); // program should not continue from here
    }
};

} // namespace seqan3::detail
