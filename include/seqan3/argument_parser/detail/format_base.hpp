// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Contains the format_base struct containing all helper functions
 *        that are needed in all formats.
 */

 #pragma once

#include <iostream>
#include <limits>
#include <sstream>
#include <string>

#include <meta/meta.hpp>

#include <seqan3/argument_parser/auxiliary.hpp>
#include <seqan3/argument_parser/exceptions.hpp>
#include <seqan3/argument_parser/validators.hpp>

namespace seqan3::detail
{

/*!\brief The format that contains all helper functions needed in all formats for
 *        printing the interface description of the application (to std::cout).
 * \ingroup argument_parser
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
        using types = meta::list<int8_t,
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
                                 std::string>;
        std::vector<std::string> names{"INT (8 bit)",
                                       "UNSIGNED (8 bit)",
                                       "INT (16 bit)",
                                       "UNSIGNED (16 bit)",
                                       "INT (32 bit)",
                                       "UNSIGNED (32 bit)",
                                       "INT (64 bit)",
                                       "UNSIGNED (64 bit)",
                                       "DOUBLE",
                                       "FLOAT",
                                       "BOOL",
                                       "CHAR",
                                       "STRING"};

        if constexpr (meta::in<types, type>::value)
            return names[meta::find_index<types, type>::value];
        else
            return "UNKNOWN_TYPE";
    }

    /*!\brief Returns the `value_type` of the input container as a string (reflection).
     * \tparam container_type The container type for which to query it's value_type.
     * \returns The type of the container value_type as a string.
     */
    template <sequence_container_concept container_type>
    //!\cond
        requires !std::is_same_v<container_type, std::string>
    //!\endcond
    static std::string get_type_name_as_string(container_type const & /**/)
    {
        typename container_type::value_type tmp;
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
        return ("\\fI" + get_type_name_as_string(value) + "\\fP");
    }

    /*!\brief Formats the container and its value_type for the help page printing.
     * \tparam container_type A type that must satisfy the seqan3::sequence_container_concept.
     * \param[in] container The container to deduct the type from.
     *
     * \returns The type of the container value type as a string, encapsulated in "List of".
     */
    template <typename container_type>
    //!\cond
        requires sequence_container_concept<container_type> && !std::is_same_v<container_type, std::string>
    //!\endcond
    static std::string option_type_and_list_info(container_type const & container)
    {
        return ("List of \\fI" + get_type_name_as_string(container) + "\\fP's");
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
            tmp.append("-" + std::string(1, *it) + ", ");

        tmp.erase(tmp.find_last_of(',')); // remove last ', '
        tmp.append(" and -" + std::string(1, flag_cluster[flag_cluster.size() - 1]));

        return tmp;
    }
};

} // namespace seqan3::detail
