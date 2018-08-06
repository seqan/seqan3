// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
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
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Contains the format_man struct and its helper functions.
 */

#pragma once

#include <iostream>

#include <seqan3/argument_parser/detail/format_base.hpp>
#include <seqan3/core/detail/terminal.hpp>
#include <seqan3/version.hpp>

namespace seqan3::detail
{

/*!\brief The format that prints the help page information formatted for a man page to std::cout.
 * \ingroup argument_parser
 *
 * \details
 *
 * The help page printing is not done immediately, because the user might not
 * provide meta information, positional options, etc. in the correct order.
 * In addition the needed order would be different from the parse format.
 * Thus the calls are stored (parser_set_up_calls and positional_option_calls)
 * and only evaluated when calling seqan3::format_help::parse.
 */
class format_man : public format_base
{
public:
    /*!\brief Adds a seqan3::detail::print_list_item call to be evaluated later on.
     *
     * \param[out] value     The variable in which to store the given command line argument.
     * \param[in]  short_id  The short identifier for the option (e.g. `'i'`).
     * \param[in]  long_id   The long identifier for the option (e.g. `"integer"`).
     * \param[in]  desc      The description of the option.
     * \param[in]  spec      Advanced option specification. see seqan3::option_spec.
     * \param[in]  validator The validator applied to the value after parsing (callable).
     */
    template <typename option_t, typename validator_type>
    void add_option(option_t & value,
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

    /*!\brief Adds a seqan3::detail::print_list_item call to be evaluated later on.
     *
     * \param[in]  short_id The short identifier for the flag (e.g. `'i'`).
     * \param[in]  long_id  The long identifier for the flag (e.g. `"integer"`).
     * \param[in]  desc     The description of the flag (which is only needed for the interface here).
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

    /*!\brief Adds a seqan3::detail::print_list_item call to be evaluated later on.
     *
     * \param[out] value     The variable in which to store the given command line argument.
     * \param[in]  desc      The description of the positional option (which is only needed
     *                       for the interface here).
     * \param[in]  validator The validator applied to the value after parsing (callable).
     */
    template <typename option_t, typename validator_type>
    void add_positional_option(option_t & value,
                               std::string const & desc,
                               validator_type && validator)
    {
        ++positional_option_count;

        positional_option_calls.push_back([=, &value]()
        {
            std::string key{"\\fBARGUMENT-" + std::to_string(positional_option_count) +
                            "\\fP " + option_type_and_list_info(value)};
            print_list_item(key, (desc + " " + validator.get_help_page_message()));
        });
    }

    /*!\brief Initiates the printing of the help page to std::cout.
     * \param[in] parser_meta The meta information that are needed for a detailed help page.
     */
    void parse(argument_parser_meta_data const & parser_meta)
    {
        print_header();

        print_section("Synopsis");
        print_synopsis();

        if (!parser_meta.description.empty())
        {
            print_section("Description");
            for (auto desc : parser_meta.description)
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

    /*!\brief Adds a print_section call to seqan3::format_man::parser_set_up_calls.
     * \param[in] title The title of the section to print.
     */
    void add_section(std::string const & title)
    {
        parser_set_up_calls.push_back([&] ()
        {
            print_section(title);
        });
    }

    /*!\brief Adds a print_subsection call to seqan3::format_man::parser_set_up_calls.
     * \param[in] title The title of the subsection to print.
     */
    void add_subsection(std::string const & title)
    {
        parser_set_up_calls.push_back([&] ()
        {
            print_subsection(title);
        });
    }

    /*!\brief Adds a print_line call to seqan3::format_man::parser_set_up_calls.
     * \param[in] text              The text to print.
     * \param[in] line_is_paragraph Whether to insert a new line at the end of the text or not.
     */
    void add_line(std::string const & text, bool line_is_paragraph)
    {
        parser_set_up_calls.push_back([&] ()
        {
            print_line(text, line_is_paragraph);
        });
    }

    /*!\brief Adds a seqan3::detail::print_list_item call to seqan3::format_man::parser_set_up_calls.
     * \param[in] key         The name that will be listed on the left hand side.
     * \param[in] description The description that will be listed on the right
     *                        hand side.
     */
    void add_list_item(std::string const & key, std::string const & description)
    {
        parser_set_up_calls.push_back([&] ()
        {
            print_list_item(key, description);
        });
    }

private:
    //!\brief Prints a help page header in man page format to std::cout.
    void print_header()
    {
        std::ostream_iterator<char> out(std::cout);

        // Print .TH line.
        std::cout << ".TH ";
        std::transform(meta.app_name.begin(), meta.app_name.end(), out, [] (unsigned char c) { return std::toupper(c); });
        std::cout << " " << std::to_string(meta.man_page_section) << " \"" << meta.date << "\" \"";
        std::transform(meta.app_name.begin(), meta.app_name.end(), out, [] (unsigned char c) { return std::tolower(c); });
        std::cout << " " << meta.version << "\" \"" << meta.man_page_title << "\"\n";

        // Print NAME section.
        std::cout << ".SH NAME\n"
                  << meta.app_name << " \\- " << meta.short_description << std::endl;
    }

    //!\brief Prints a help page synopsis in man page format section to std::cout.
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

    /*!\brief Prints a section title in man page format to std::cout.
     * \param[in] title The title of the section to print.
     */
    void print_section(std::string const & title)
    {
        std::ostream_iterator<char> out(std::cout);
        std::cout << ".SH ";
        std::transform(title.begin(), title.end(), out, [] (unsigned char c) { return std::toupper(c); });
        std::cout << "\n";
        is_first_in_section = true;
    }

    /*!\brief Prints a subsection title in man page format to std::cout.
     * \param[in] title The title of the subsection to print.
     */
    void print_subsection(std::string const & title)
    {
            std::cout << ".SS " << title << "\n";
            is_first_in_section = true;
    }

    /*!\brief Prints a help page section in man page format to std::cout.
     *
     * \param[in] text The text to print.
     * \param[in] line_is_paragraph Whether to insert as paragraph
     *            or just a line (only one line break if not a paragraph).
     */
    void print_line(std::string const & text, bool const line_is_paragraph)
    {
        if (!is_first_in_section && line_is_paragraph)
            std::cout << ".sp\n";
        else if (!is_first_in_section && !line_is_paragraph)
            std::cout << ".br\n";

        std::cout << text << "\n";
        is_first_in_section = false;
    }

    /*!\brief Prints a text in man page format to std::cout.
     * \param[in] text The text to print.
     */
    void print_line(std::string const & text)
    {
        print_line(text, true);
    }

    /*!\brief Prints a help page list_item in man page format to std::cout.
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
            std::cout << ".TP\n"
                   << term << "\n"
                   << desc << "\n";
            is_first_in_section = false;
    }

    //!\brief Prints a help page footer in man page format.
    void print_footer()
    {
        // Print legal stuff
        if ((!empty(meta.short_copyright)) || (!empty(meta.long_copyright)) || (!empty(meta.citation)))
        {
            std::cout << ".SH LEGAL\n";

            if (!empty(meta.short_copyright))
                std::cout << "\\fB" << meta.app_name << " Copyright:\\fR " << meta.short_copyright << "\n.br\n";

            std::cout << "\\fBSeqAn Copyright:\\fR 2006-2015 Knut Reinert, FU-Berlin; released under the 3-clause BSDL.\n.br\n";

            if (!empty(meta.citation))
                std::cout << "\\fBIn your academic works please cite:\\fR " << meta.citation << "\n.br\n";

            if (!empty(meta.long_copyright))
                std::cout << "For full copyright and/or warranty information see \\fB--copyright\\fR.\n";
        }
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
    //!\brief Needed for correct indentation and line breaks.
    bool is_first_in_section{true};
    //!\brief Keeps track of the number of positional options.
    unsigned positional_option_count{0};
};

} // namespace seqan3
