// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Provides the format_man struct and its helper functions.
 */

#pragma once

#include <iostream>

#include <seqan3/argument_parser/detail/format_base.hpp>
#include <seqan3/argument_parser/detail/terminal.hpp>
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
 * and only evaluated when calling seqan3::detail::format_help_base::parse.
 */
class format_man : public format_help_base<format_man>
{
    //!\brief The CRTP base class type.
    using base_type = format_help_base<format_man>;

    //!\brief Befriend the base class to give access to the private member functions.
    friend base_type;

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    format_man() = default;                                  //!< Defaulted.
    format_man(format_man const & pf) = default;             //!< Defaulted.
    format_man & operator=(format_man const & pf) = default; //!< Defaulted.
    format_man(format_man &&) = default;                     //!< Defaulted.
    format_man & operator=(format_man &&) = default;         //!< Defaulted.
    ~format_man() = default;                                 //!< Defaulted.

    //!\copydoc format_help_base(std::vector<std::string> const &, bool const)
    format_man(std::vector<std::string> const & names, bool const advanced = false) : base_type{names, advanced}
    {};
    //!\}

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

    //!\brief Needed for correct indentation and line breaks.
    bool is_first_in_section{true};
};

} // namespace seqan3
