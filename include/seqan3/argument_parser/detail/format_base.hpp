// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Provides the format_base struct containing all helper functions
 *        that are needed in all formats.
 */

#pragma once

#include <filesystem>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include <seqan3/argument_parser/auxiliary.hpp>
#include <seqan3/argument_parser/detail/concept.hpp>
#include <seqan3/argument_parser/exceptions.hpp>
#include <seqan3/argument_parser/validators.hpp>
#include <seqan3/utility/detail/type_name_as_string.hpp>
#include <seqan3/version.hpp>

namespace seqan3::detail
{

/*!\brief The format that contains all helper functions needed in all formats.
 * \ingroup argument_parser
 * \remark For a complete overview, take a look at \ref argument_parser
 */
class format_base
{
protected:
    /*!\brief Returns the input type as a string (reflection).
     * \tparam value_type The type whose name is converted std::string.
     * \returns The type of the value as a string.
     */
    template <typename value_type>
    static std::string get_type_name_as_string(value_type const & /**/)
    {
        using type = std::decay_t<value_type>;
        using types = type_list<int8_t,
                                uint8_t,
                                int16_t,
                                uint16_t,
                                int32_t,
                                uint32_t,
                                int64_t,
                                uint64_t,
                                double,
                                float,
                                bool,
                                char,
                                std::string,
                                std::filesystem::path>;
        std::vector<std::string> names{"signed 8 bit integer",
                                       "unsigned 8 bit integer",
                                       "signed 16 bit integer",
                                       "unsigned 16 bit integer",
                                       "signed 32 bit integer",
                                       "unsigned 32 bit integer",
                                       "signed 64 bit integer",
                                       "unsigned 64 bit integer",
                                       "double",
                                       "float",
                                       "bool",
                                       "char",
                                       "std::string",
                                       "std::filesystem::path"};

        if constexpr (list_traits::contains<type, types>)
            return names[list_traits::find<type, types>];
        else
            return detail::type_name_as_string<value_type>;
    }

    /*!\brief Returns the `value_type` of the input container as a string (reflection).
     * \tparam container_type The container type for which to query it's value_type.
     * \returns The type of the container value_type as a string.
     */
    template <detail::is_container_option container_type>
    static std::string get_type_name_as_string(container_type const & /**/)
    {
        typename container_type::value_type tmp{};
        return get_type_name_as_string(tmp);
    }

    /*!\brief Formats the type of a value for the help page printing.
     * \tparam option_value_type The type of the option value to get the info for.
     * \param[in] value The value to deduct the type from.
     * \returns The type of the value as string.
     */
    template <typename option_value_type>
    static std::string option_type_and_list_info(option_value_type const & value)
    {
        return ("(\\fI" + get_type_name_as_string(value) + "\\fP)");
    }

    /*!\brief Formats the container and its value_type for the help page printing.
     * \tparam container_type A type that must satisfy the seqan3::detail::is_container_option.
     * \param[in] container The container to deduct the type from.
     *
     * \returns The type of the container value type as a string, encapsulated in "List of".
     */
    template <detail::is_container_option container_type>
    static std::string option_type_and_list_info(container_type const & container)
    {
        return ("(\\fIList\\fP of \\fI" + get_type_name_as_string(container) + "\\fP)");
    }

    /*!\brief Formats the option/flag identifier pair for the help page printing.
     * \param[in] short_id The short identifier of the option/flag.
     * \param[in] long_id  The long identifier of the option/flag.
     * \returns The name of the short and long id, prepended with (double)dash.
     *
     * \details  e.g. "-i,--integer", "-i", or "--integer".
     */
    static std::string prep_id_for_help(char const short_id, std::string const & long_id)
    {
        // Build list item term.
        std::string term;
        if (short_id != '\0')
            term = "\\fB-" + std::string(1, short_id) + "\\fP";

        if (short_id != '\0' && !long_id.empty())
            term.append(", ");

        if (!long_id.empty())
            term.append("\\fB--" + long_id + "\\fP");

        return term;
    }

    /*!\brief Escapes certain characters for correct output.
     * \param[in] original The string containing characters to be escaped.
     * \returns The original string as their corresponding html/xml representation.
     *
     * \details Special characters considered are `"`, `\`, `&`, `<` and `>`.
     */
    std::string escape_special_xml_chars(std::string const & original)
    {
        std::string escaped;
        escaped.reserve(original.size()); // will be at least as long

        for (auto c : original)
        {
            if (c == '"')
                escaped.append("&quot;");
            else if (c == '\'')
                escaped.append("&apos;");
            else if (c == '&')
                escaped.append("&amp;");
            else if (c == '<')
                escaped.append("&lt;");
            else if (c == '>')
                escaped.append("&gt;");
            else
                escaped.push_back(c);
        }

        return escaped;
    }

    /*!\brief Expands multiple one character flag identifiers for pretty help output.
     * \param[in] flag_cluster The string of one character flags.
     * \returns A string that lists all flags as a comma separated list.
     *
     * \details e.g. "-agdg" becomes "-a, -g, -d and -g".
     */
    static std::string expand_multiple_flags(std::string const & flag_cluster)
    {
        std::string tmp;
        auto it{flag_cluster.begin()};

        if (flag_cluster[0] == '-')
            ++it;

        for (; it != flag_cluster.end() - 1; ++it)
            tmp.append({'-', *it, ',', ' '});

        tmp.erase(tmp.find_last_of(',')); // remove last ', '
        tmp.append({'a', 'n', 'd', ' ', '-', flag_cluster[flag_cluster.size() - 1]});

        return tmp;
    }
};

/*!\brief The format that contains all helper functions needed in all formats for
 *        printing the interface description of the application (to std::cout).
 * \ingroup argument_parser
 * \remark For a complete overview, take a look at \ref argument_parser
 */
template <typename derived_type>
class format_help_base : public format_base
{
private:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    format_help_base() = default;                                        //!< Defaulted.
    format_help_base(format_help_base const & pf) = default;             //!< Defaulted.
    format_help_base & operator=(format_help_base const & pf) = default; //!< Defaulted.
    format_help_base(format_help_base &&) = default;                     //!< Defaulted.
    format_help_base & operator=(format_help_base &&) = default;         //!< Defaulted.
    ~format_help_base() = default;                                       //!< Defaulted.

    /*!\brief Initializes a format_help_base object.
     * \param[in] names    A list of subcommands (see \link subcommand_arg_parse subcommand parsing \endlink).
     * \param[in] advanced Set to `true` to show advanced options.
     */
    format_help_base(std::vector<std::string> const & names, bool const advanced) :
        command_names{names},
        show_advanced_options{advanced}
    {}
    //!\}

public:
    /*!\brief Adds a seqan3::print_list_item call to be evaluated later on.
     * \copydetails seqan3::argument_parser::add_option
     */
    template <typename option_type, typename validator_type>
    void add_option(option_type & value,
                    char const short_id,
                    std::string const & long_id,
                    std::string const & desc,
                    option_spec const spec,
                    validator_type && option_validator)
    {
        std::string id = prep_id_for_help(short_id, long_id) + " " + option_type_and_list_info(value);
        std::string info{desc};
        info += ((spec & option_spec::required) ? std::string{" "} : detail::to_string(" Default: ", value, ". "));
        info += option_validator.get_help_page_message();
        store_help_page_element(
            [this, id, info]()
            {
                derived_t().print_list_item(id, info);
            },
            spec);
    }

    /*!\brief Adds a seqan3::print_list_item call to be evaluated later on.
     * \copydetails seqan3::argument_parser::add_flag
     */
    void add_flag(bool & SEQAN3_DOXYGEN_ONLY(value),
                  char const short_id,
                  std::string const & long_id,
                  std::string const & desc,
                  option_spec const spec)
    {
        std::string id = prep_id_for_help(short_id, long_id);
        store_help_page_element(
            [this, id, desc]()
            {
                derived_t().print_list_item(id, desc);
            },
            spec);
    }

    /*!\brief Adds a seqan3::print_list_item call to be evaluated later on.
     * \copydetails seqan3::argument_parser::add_positional_option
     */
    template <typename option_type, typename validator_type>
    void add_positional_option(option_type & value, std::string const & desc, validator_type & option_validator)
    {
        std::string msg = option_validator.get_help_page_message();

        positional_option_calls.push_back(
            [this, &value, desc, msg]()
            {
                ++positional_option_count;
                derived_t().print_list_item(detail::to_string("\\fBARGUMENT-",
                                                              positional_option_count,
                                                              "\\fP ",
                                                              option_type_and_list_info(value)),
                                            desc +
                                                // a list at the end may be empty and thus have a default value
                                                ((detail::is_container_option<option_type>)
                                                     ? detail::to_string(" Default: ", value, ". ")
                                                     : std::string{" "})
                                                + msg);
            });
    }

    /*!\brief Initiates the printing of the help page to std::cout.
     * \param[in] parser_meta The meta information that are needed for a detailed help page.
     */
    void parse(argument_parser_meta_data & parser_meta)
    {
        meta = parser_meta;

        derived_t().print_header();

        if (!meta.synopsis.empty())
        {
            derived_t().print_section("Synopsis");
            derived_t().print_synopsis();
        }

        if (!meta.description.empty())
        {
            derived_t().print_section("Description");
            for (auto desc : meta.description)
                print_line(desc);
        }

        if (!command_names.empty())
        {
            derived_t().print_section("Subcommands");
            derived_t().print_line("This program must be invoked with one of the following subcommands:", false);
            for (std::string const & name : command_names)
                derived_t().print_line("- \\fB" + name + "\\fP", false);
            derived_t().print_line("See the respective help page for further details (e.g. by calling " + meta.app_name
                                       + " " + command_names[0] + " -h).",
                                   true);
            derived_t().print_line("The following options below belong to the top-level parser and need to be "
                                   "specified \\fBbefore\\fP the subcommand key word. Every argument after the "
                                   "subcommand key word is passed on to the corresponding sub-parser.",
                                   true);
        }

        // add positional options if specified
        if (!positional_option_calls.empty())
            derived_t().print_section("Positional Arguments");

        // each call will evaluate the function derived_t().print_list_item()
        for (auto f : positional_option_calls)
            f();

        // add options and flags if specified
        if (!parser_set_up_calls.empty())
            derived_t().print_section("Options");

        // each call will evaluate the function derived_t().print_list_item()
        for (auto f : parser_set_up_calls)
            f();

        if (!meta.examples.empty())
        {
            derived_t().print_section("Examples");
            for (auto example : meta.examples)
                print_line(example);
        }

        print_version();

        print_legal();

        derived_t().print_footer();

        std::exit(EXIT_SUCCESS); // program should not continue from here
    }

    /*!\brief Adds a print_section call to parser_set_up_calls.
     * \copydetails seqan3::argument_parser::add_section
     */
    void add_section(std::string const & title, option_spec const spec)
    {
        store_help_page_element(
            [this, title]()
            {
                derived_t().print_section(title);
            },
            spec);
    }

    /*!\brief Adds a print_subsection call to parser_set_up_calls.
     * \copydetails seqan3::argument_parser::add_subsection
     */
    void add_subsection(std::string const & title, option_spec const spec)
    {
        store_help_page_element(
            [this, title]()
            {
                derived_t().print_subsection(title);
            },
            spec);
    }

    /*!\brief Adds a print_line call to parser_set_up_calls.
     * \copydetails seqan3::argument_parser::add_line
     */
    void add_line(std::string const & text, bool is_paragraph, option_spec const spec)
    {
        store_help_page_element(
            [this, text, is_paragraph]()
            {
                derived_t().print_line(text, is_paragraph);
            },
            spec);
    }

    /*!\brief Adds a seqan3::print_list_item call to parser_set_up_calls.
     * \copydetails seqan3::argument_parser::add_list_item
     */
    void add_list_item(std::string const & key, std::string const & desc, option_spec const spec)
    {
        store_help_page_element(
            [this, key, desc]()
            {
                derived_t().print_list_item(key, desc);
            },
            spec);
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

    //!\brief Befriend the derived type so it can access private functions.
    friend derived_type;

protected:
    //!\brief Returns the derived type.
    derived_type & derived_t()
    {
        return static_cast<derived_type &>(*this);
    }

    //!\brief Prints a synopsis in any format.
    void print_synopsis()
    {
        for (unsigned i = 0; i < meta.synopsis.size(); ++i)
        {
            std::string text = "\\fB";
            text.append(meta.synopsis[i]);
            text.insert(text.find_first_of(" \t"), "\\fP");

            derived_t().print_line(text, false);
        }
    }

    /*!\brief Delegates to seqan3::print_line(std::string const & text, true) of each format.
     * \param[in] text The text to print.
     */
    void print_line(std::string const & text)
    {
        derived_t().print_line(text, true);
    }

    //!\brief Prints the version information.
    void print_version()
    {
        std::string const version_str{seqan3_version_cstring};

        // Print version, date and url.
        derived_t().print_section("Version");
        derived_t().print_line(derived_t().in_bold("Last update: ") + meta.date, false);
        derived_t().print_line(derived_t().in_bold(meta.app_name + " version: ") + meta.version, false);
        derived_t().print_line(derived_t().in_bold("SeqAn version: ") + version_str, false);

        if (!empty(meta.url))
        {
            derived_t().print_section("Url");
            derived_t().print_line(meta.url, false);
        }
    }

    //!\brief Prints the legal information.
    void print_legal()
    {
        // Print legal stuff
        if ((!empty(meta.short_copyright)) || (!empty(meta.long_copyright)) || (!empty(meta.citation))
            || (!empty(meta.author)) || (!empty(meta.email)))
        {
            derived_t().print_section("Legal");

            if (!empty(meta.short_copyright))
            {
                derived_t().print_line(derived_t().in_bold(meta.app_name + " Copyright: ") + meta.short_copyright,
                                       false);
            }

            if (!empty(meta.author))
            {
                derived_t().print_line(derived_t().in_bold("Author: ") + meta.author, false);
            }

            if (!empty(meta.email))
            {
                derived_t().print_line(derived_t().in_bold("Contact: ") + meta.email, false);
            }

            derived_t().print_line(derived_t().in_bold("SeqAn Copyright: ")
                                       + "2006-2025 Knut Reinert, FU-Berlin; released under the 3-clause BSDL.",
                                   false);

            if (!empty(meta.citation))
            {
                derived_t().print_line(derived_t().in_bold("In your academic works please cite: ") + meta.citation,
                                       false);
            }

            if (!empty(meta.long_copyright))
            {
                derived_t().print_line("For full copyright and/or warranty information see "
                                           + derived_t().in_bold("--copyright") + ".",
                                       false);
            }
        }
    }

    //!\brief Vector of functions that stores all calls except add_positional_option.
    std::vector<std::function<void()>> parser_set_up_calls;
    //!\brief Vector of functions that stores add_positional_option calls.
    std::vector<std::function<void()>> positional_option_calls; // singled out to be printed on top
    //!\brief Keeps track of the number of positional options
    unsigned positional_option_count{0};
    //!\brief The names of subcommand programs.
    std::vector<std::string> command_names{};
    //!\brief Whether to show advanced options or not.
    bool show_advanced_options{true};

private:
    /*!\brief Adds a function object to parser_set_up_calls **if** the annotation in `spec` does not prevent it.
     * \param[in] printer The invokable that, if added to `parser_set_up_calls`, prints information to the help page.
     * \param[in] spec The option specification deciding whether to add the information to the help page.
     *
     * \details
     *
     * If `spec` equals `seqan3::option_spec::hidden`, the information is never added to the help page.
     * If `spec` equals `seqan3::option_spec::advanced`, the information is only added to the help page if
     * the advanced help page has been queried on the command line (`show_advanced_options == true`).
     */
    void store_help_page_element(std::function<void()> printer, option_spec const spec)
    {
        if (!(spec & option_spec::hidden) && (!(spec & option_spec::advanced) || show_advanced_options))
            parser_set_up_calls.push_back(std::move(printer));
    }
};

} // namespace seqan3::detail
