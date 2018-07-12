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
 * \brief Contains the format_parse class.
 */

#pragma once

#include <sstream>
#include <string>
#include <vector>

#include <seqan3/argument_parser/detail/format_base.hpp>
#include <seqan3/std/concept/core_language.hpp>
#include <seqan3/std/concept/container.hpp>

namespace seqan3::detail
{

/*!\brief The format that organizes the actual parsing of command line arguments.
 * \ingroup argument_parser
 *
 * \details
 *
 * In order to be independent of the options value type, we do not want to store
 * parameters/options/flags/.. directly (though a variant might work, it is hacky).
 * Directly parsing is also difficult, since the order of parsing options/flags
 * is non trivial (e.g. ambiguousness of '-g 4' => option+value or flag+positional).
 * Therefore, we store the parsing calls of the developer in a function object,
 * (format_parse::option_and_flag_calls, seqan3::detail::format_parse::positional_option_calls)
 * executing them in a new order when calling format_parse::parse().
 * This enables us to parse any option type and resolve any ambiguousness, so no
 * additional restrictions apply to the developer when setting up the parser.
 *
 * Order of parsing:
 * -#. Options            (order within as specified by the developer)
 * -#. Flags              (order within as specified by the developer)
 * -#. Positional Options (order within as specified by the developer)
 *
 * When parsing flags and options, the identifiers (and values) are removed from
 * the vector format_parse::argv. That way, options that are specified multiple times,
 * but are no container type, can be identified and an error is reported.
 */
class format_parse : public format_base
{
public:
    /*!\name Constructors, destructor and assignment
     * \brief The (default) constructors.
     * \{
     */
    format_parse() = delete;
    format_parse(format_parse const & pf) = default;
    format_parse & operator=(format_parse const & pf) = default;
    format_parse(format_parse &&) = default;
    format_parse & operator=(format_parse &&) = default;

    /*!\brief The constructor of the parse format.
     * \param[in] argc_ The number of command line arguments.
     * \param[in] argv_ The command line arguments to parse.
     */
    format_parse(int const argc_, const char ** argv_) :
        argc(argc_ - 1)
    {
        init(argc_, argv_);
    }

    ~format_parse() = default;
    //!\}

    /*!\brief Adds an seqan3::detail::get_option call to be evaluated later on.
     *
     * \tparam option_type    The type of variable in which to store the given command line argument.
     * \tparam validator_type The type of validator applied to the value after parsing.
     *
     * \param[out] value     The variable in which to store the given command line argument.
     * \param[in]  short_id  The short identifier for the option (e.g. `'i'`).
     * \param[in]  long_id   The long identifier for the option (e.g. `"integer"`).
     * \param[in]  spec      Advanced option specification. see seqan3::option_spec.
     * \param[in]  validator The validator applied to the value after parsing (callable).
     *
     * \throws seqan3::parser_design_error
     */
    template <typename option_type, typename validator_type>
    void add_option(option_type & value,
                    char const short_id,
                    std::string const & long_id,
                    std::string const & /*desc*/,
                    option_spec const & spec,
                    validator_type && validator)
    {
        option_calls.push_back([=, &value]()
        {
            get_option(value, short_id, long_id, spec, validator);
        });
    }

    /*!\brief Adds a get_flag call to be evaluated later on.
     *
     * \param[out] value    The variable in which to store the given command line argument.
     * \param[in]  short_id The short identifier for the flag (e.g. `'i'`).
     * \param[in]  long_id  The long identifier for the flag (e.g. `"integer"`).
     *
     * \throws seqan3::parser_design_error
     */
    void add_flag(bool & value,
                  char const short_id,
                  std::string const & long_id,
                  std::string const & /*desc*/,
                  option_spec const & /*spec*/)
    {
        flag_calls.push_back([=, &value]()
        {
            get_flag(value, short_id, long_id);
        });
    }

    /*!\brief Adds a get_positional_option call to be evaluated later on.
     *
     * \tparam option_type    The type of variable in which to store the given command line argument.
     * \tparam validator_type The type of validator applied to the value after parsing.
     *
     * \param[out] value     The variable in which to store the given command line argument.
     * \param[in]  validator The validator applied to the value after parsing (callable).
     *
     * \throws seqan3::parser_design_error
     */
    template <typename option_type, typename validator_type>
    void add_positional_option(option_type & value,
                               std::string const & /*desc*/,
                               validator_type && validator)
    {
        positional_option_calls.push_back([=, &value]()
        {
            get_positional_option(value, validator);
        });
    }

    //!\brief Initiates the actual command line parsing.
    void parse(argument_parser_meta_data const & /*meta*/)
    {
        end_of_options_it = std::find(argv.begin(), argv.end(), "--");

        // parse options first, because we need to rule out -keyValue pairs
        // (e.g. -AnoSpaceAfterIdentifierA) before parsing flags
        for (auto && f : option_calls)
            f();

        for (auto && f : flag_calls)
            f();

        check_for_unkown_ids();

        if (end_of_options_it != argv.end())
            *end_of_options_it = ""; // remove -- before parsing positional arguments

        for (auto && f : positional_option_calls)
            f();

        check_for_left_over_args();
    }

    // functions are not needed for command line parsing but are part of the format interface.
    //!\cond
    void add_section(std::string const &) {}
    void add_subsection(std::string const &) {}
    void add_line(std::string const &, bool) {}
    void add_list_item(std::string const &, std::string const &) {}
    //!\endcond

private:
    /*!\brief Initializes the format_parse on construction.
     *
     * \param[in] argc_ Number of command line arguments.
     * \param[in] argv_ Vector of command line arguments.
     *
     * \details
     * Adds all command line arguments to the member format_parse::argv,
     * but splits all values (options) containing an equality sign.
     */
    void init(int argc_, const char ** argv_)
    {
        argv.resize(argc_ - 1); // -1 because of the binary name

        for(int i = 1; i < argc_; ++i) // start at 1 to skip binary name
            argv.push_back(argv_[i]);
    }

    /*!\brief Appends a double dash to a long identifier and returns it.
    * \param[in] long_id The name of the long identifier.
    * \returns The input long name prepended with a double dash.
    */
    std::string prepend_dash(std::string const & long_id)
    {
        return ("--" + long_id);
    }

    /*!\brief Appends a double dash to a short identifier and returns it.
    * \param[in] short_id The name of the short identifier.
    * \returns The input short name prepended with a single dash.
    */
    std::string prepend_dash(char const short_id)
    {
        return ("-" + std::string(1, short_id));
    }

    /*!\brief Finds the position of a short/long identifier in format_parse::argv.
     *
     * \param[in] begin_it The iterator where to start the search of the identifier.
     *                     Note that the end iterator is kept as a member variable.
     * \param[in] id       The identifier to search for (must not contain dashes).
     * \returns An iterator pointing to the first occurrence of the identifier in
     *          the list pointed to by begin_t. If the list does not contain the
     *          identifier `id`, the member variable `end_of_options_it` is returned.
     *
     * Note: The `id` is compared to the prefix of each value in the list, such
     *       that "-idValue" arguments are correctly identified.
     */
    template <typename id_type>
    std::vector<std::string>::iterator find_option_id(std::vector<std::string>::iterator const begin_it, id_type const & id)
    {
        return(std::find_if(begin_it, end_of_options_it,
            [&] (const std::string & v)
            {
                size_t id_size{(prepend_dash(id)).size()};
                if (v.size() < id_size)
                    return false; // cannot be the correct identifier

                return v.substr(0, id_size) == prepend_dash(id); // check if prefix of v is the same
            }));
    }

    /*!\brief Returns true and removes the long identifier if it is in format_parse::argv.
     * \param[in] long_id The long identifier of the flag to check.
     */
    bool flag_is_set(std::string const & long_id)
    {
        auto it = std::find(argv.begin(), end_of_options_it, prepend_dash(long_id));

        if (it != end_of_options_it)
            *it = ""; // remove seen flag

        return(it != end_of_options_it);
    }

    /*!\brief Returns true and removes the short identifier if it is in format_parse::argv.
     * \param[in] short_id The short identifier of the flag to check.
     */
    bool flag_is_set(char const short_id)
    {
        // short flags need special attention, since they could be grouped (-rGv <=> -r -G -v)
        for (std::string & arg : argv)
        {
            if (arg[0] == '-' && arg.size() > 1 && arg[1] != '-') // is option && not dash && no long option
            {
                auto pos = arg.find(short_id);

                if (pos != std::string::npos)
                {
                    arg.erase(pos, 1); // remove seen bool

                    if (arg == "-") // if flag is empty now
                        arg = "";

                    return true;
                }
            }
        }
        return false;
    }

    /*!\brief Tries to cast an input string into a value.
     *
     * \tparam option_t Must satisfy the seqan3::istream_concept.
     *
     * \param[out] value Stores the casted value.
     * \param[in]  in    The input argument to be casted.
     *
     * \throws seqan3::parser_invalid_argument
     */
    template <typename option_t>
    //!\cond
        requires istream_concept<std::istringstream, option_t>
    //!\endcond
    void retrieve_value(option_t & value, std::string const & in)
    {
        std::istringstream stream{in};
        stream >> value;

        if (stream.fail() || !stream.eof())
            throw type_conversion_failed("Argument " + in + " could not be casted to type " +
                                         get_type_name_as_string(value) + ".");
    }

    //!\cond
    void retrieve_value(std::string & value, std::string const & in)
    {
        value = in;
    }
    //!\endcond

    /*!\brief Appends a casted value to its container.
     *
     * \tparam container_option_t Must satisfy the seqan3::sequence_container_concept and
     *                            its value_type must satisfy the seqan3::istream_concept
     *
     * \param[out] value Container that stores the casted value.
     * \param[in]  in    The input argument to be casted.
     */
    template <sequence_container_concept container_option_t>
    //!\cond
        requires istream_concept<std::istringstream, typename container_option_t::value_type>
    //!\cond
    void retrieve_value(container_option_t & value, std::string const & in)
    {
        typename container_option_t::value_type tmp;

        retrieve_value(tmp, in); // throws on failure
        value.push_back(tmp);
    }

    /*!\brief Tries to cast an input string into a (integral) value and tests for overflow errors.
     *
     * \tparam option_t Must be stream-convertible.
     *
     * \param[out] value Stores the casted value.
     * \param[in]  in    The input argument to be casted.
     *
     * \throws seqan3::type_conversion_failed
     * \throws seqan3::overflow_error_on_conversion
     *
     * \details
     *
     * This function tests for possible overflow errors:
     * It first casts into a signed long long int (to also catch negative values for unsigned ints)
     * and then checks the min-max range given by numeric limits.
     * Exception: Due to using signed long long int, this does not catch overflow errors
     * for unsigned long int's but that kind of large a number will hopefully never be supplied
     * via the command line..
     */
    template <integral_concept option_t>
    //!\cond
        requires istream_concept<std::istringstream, option_t>
    //!\endcond
    void retrieve_value(option_t & value, std::string const & in)
    {
        int64_t tmp;
        std::istringstream stream{in};
        stream >> tmp;

        if (stream.fail() || !stream.eof()) // !stream.eof() catches 13a as wrong
            throw type_conversion_failed("Argument " + in + " could not be casted to type " +
                                         get_type_name_as_string(value) + ".");

        if (tmp > static_cast<int64_t>(std::numeric_limits<option_t>::max()) ||
            tmp < static_cast<int64_t>(std::numeric_limits<option_t>::min()))
            throw overflow_error_on_conversion("Argument " + in + " is not in integer range [" +
                                               std::to_string(std::numeric_limits<option_t>::min()) + "," +
                                               std::to_string(std::numeric_limits<option_t>::max()) + "].");

        value = static_cast<option_t>(tmp);
    }

    /*!\brief Handles value retrieval for options based on different kev value pairs.
     *
     * \param[out] value     Stores the value found in argv, casted by retrieve_value.
     * \param[in]  option_it The iterator where the option identifier was found.
     * \param[in]  id        The option identifier supplied on the command line.
     *
     * \throws seqan3::parser_invalid_argument
     *
     * \details
     *
     * The value at option_it is inspected whether it is an '-key value', '-key=value'
     * or '-keyValue' pair and the input is extracted accordingly. The input
     * will then be tried to be casted into the `value` parameter.
     *
     * Returns true on success and false otherwise.
     */
    template <typename option_type, typename id_type>
    bool identify_and_retrieve_option_value(option_type & value,
                                            std::vector<std::string>::iterator & option_it,
                                            id_type const & id)
    {
        if (option_it != end_of_options_it)
        {
            std::string input_value;
            size_t id_size = (prepend_dash(id)).size();

            if ((*option_it).size() > id_size) // identifier includes value (-keyValue or -key=value)
            {
                if ((*option_it)[id_size] == '=') // -key=value
                {
                    if ((*option_it).size() == id_size + 1) // malformed because no value follows '-i='
                        throw parser_invalid_argument("Value cast failed for option " +
                                                      prepend_dash(id) +
                                                      ": No value was provided.");
                    input_value = (*option_it).substr(id_size + 1);
                }
                else // -kevValue
                {
                    input_value = (*option_it).substr(id_size);
                }

                *option_it = ""; // remove used identifier-value pair
            }
            else // -key value
            {
                *option_it = ""; // remove used identifier
                ++option_it;
                if (option_it == end_of_options_it) // should not happen
                    throw parser_invalid_argument("Value cast failed for option " +
                                                  prepend_dash(id) +
                                                  ": No value was provided.");
                input_value = *option_it;
                *option_it = ""; // remove value
            }

            try
            {
                retrieve_value(value, input_value);
            }
            catch (parser_invalid_argument const & ex)
            {
                throw parser_invalid_argument("Value cast failed for option " + prepend_dash(id) + ": " + ex.what());
            }

            return true;
        }
        return false;
    }

    /*!\brief Handles value retrieval (non container type) options.
     *
     * \param[out] value Stores the value found in argv, casted by retrieve_value.
     * \param[in]  id    The option identifier supplied on the command line.
     *
     * \throws seqan3::option_declared_multiple_times
     *
     * \details
     *
     * If the option identifier is found in format_parse::argv, the value of
     * the following position in argv is tried to be casted into value
     * and the identifier and value argument are removed from argv.
     *
     * Returns true on success and false otherwise. This is needed to catch
     * the user error of supplying multiple arguments for the same
     * (non container!) option by specifying the short AND long identifier.
     */
    template <typename option_type, typename id_type>
    bool get_option_by_id(option_type & value, id_type const & id)
    {
        auto it = find_option_id(argv.begin(), id);

        if (it != end_of_options_it)
            identify_and_retrieve_option_value(value, it, id);

        if (find_option_id(it, id) != end_of_options_it) // should not be found again
           throw option_declared_multiple_times("Option " + prepend_dash(id) +
                                                "is no list/container but declared multiple times.");

       return (it != end_of_options_it); // first search was successful or not
    }

    /*!\brief Handles value retrieval (container type) options.
     *
     * \param[out] value Stores all values found in argv, casted by retrieve_value.
     * \param[in]  id    The option identifier supplied on the command line.
     *
     * \details
     *
     * Since option_type is a container, the option is a list and can be parsed
     * multiple times.
     *
     */
    template <sequence_container_concept option_type, typename id_type>
    //!cond
        requires !std::is_same_v<option_type, std::string>
    //!\endcond
    bool get_option_by_id(option_type & value, id_type const & id)
    {
        auto it = find_option_id(argv.begin(), id);
        bool seen_at_least_once{it != end_of_options_it};

        while (it != end_of_options_it)
        {
            identify_and_retrieve_option_value(value, it, id);
            it = find_option_id(it, id);
        }

        return seen_at_least_once;
    }

    /*!\brief Checks format_parse::argv for unknown options/flags.
     *
     * \throws seqan3::unknown_option
     *
     * \details
     *
     * This function is used by format_parse::parse() AFTER all flags and options
     * specified by the developer were parsed and therefore removed from argv.
     * Thus, all remaining flags/options are unknown.
     *
     * In addition this function removes "--" (if specified) from argv to
     * clean argv for positional option retrieval.
     */
    void check_for_unkown_ids()
    {
        for (auto it = argv.begin(); it != end_of_options_it; ++it)
        {
            std::string arg{*it};
            if (!arg.empty() && arg[0] == '-') // may be an identifier
            {
                if (arg == "-")
                {
                    continue; // positional option
                }
                else if (arg[1] != '-' && arg.size() > 2) // one dash, but more than one character (-> multiple flags)
                {
                    throw unknown_option("Unknown flags " + expand_multiple_flags(arg) +
                                         ". In case this is meant to be a non-option/argument/parameter, " +
                                         "please specify the start of arguments with '--'. " +
                                         "See -h/--help for program information.");
                }
                else // unknown short or long option
                {
                    throw unknown_option("Unknown option " + arg +
                                         ". In case this is meant to be a non-option/argument/parameter, " +
                                         "please specify the start of non-options with '--'. " +
                                         "See -h/--help for program information.");
                }
            }
        }
    }

    /*!\brief Checks format_parse::argv for unknown options/flags.
     *
     * \throws seqan3::too_many_arguments
     *
     * \details
     *
     * This function is used by format_parse::parse() AFTER all flags, options
     * and positional options specified by the developer were parsed and
     * therefore removed from argv.
     * Thus, all remaining non-empty arguments are too much.
     */
    void check_for_left_over_args()
    {
        if (std::find_if(argv.begin(), argv.end(), [](std::string const & s){return (s != "");}) != argv.end())
            throw too_many_arguments("Too many arguments provided. Please see -h/--help for more information.");
    }

    /*!\brief Handles command line option retrieval.
     *
     * \param[out] value     The variable in which to store the given command line argument.
     * \param[in]  short_id  The short identifier for the option (e.g. `'i'`).
     * \param[in]  long_id   The long identifier for the option (e.g. `"integer"`).
     * \param[in]  spec      Advanced option specification. see seqan3::option_spec.
     * \param[in]  validator The validator applied to the value after parsing (callable).
     *
     * \throws seqan3::option_declared_multiple_times
     * \throws seqan3::validation_failed
     * \throws seqan3::required_option_missing
     *
     * \details
     *
     * This function
     * - checks if the option is required but not set,
     * - retrieves any value found by the short or long identifier,
     * - throws on (mis)use of both identifiers for non-container type values,
     * - re-throws the validation exception with appended option information
     */
    template <typename option_type, typename validator_type>
    void get_option(option_type & value,
                     char const short_id,
                     std::string const & long_id,
                     option_spec const & spec,
                     validator_type && validator)
    {
        bool short_id_is_set{get_option_by_id(value, short_id)};
        bool long_id_is_set{get_option_by_id(value, long_id)};

        // if value is no container we need to check for multiple declarations
        if (short_id_is_set && long_id_is_set &&
            !(sequence_container_concept<option_type> && !std::is_same_v<option_type, std::string>))
            throw option_declared_multiple_times("Option " + prepend_dash(short_id) + "/" + prepend_dash(long_id) +
                                                 " is no list/container but specified multiple times");

        if (short_id_is_set || long_id_is_set)
        {
            try
            {
                validator(value);
            }
            catch (std::exception & ex)
            {
                throw validation_failed(std::string("Validation failed for option ") + prepend_dash(short_id) + "/" +
                                        prepend_dash(long_id) + ": " + ex.what());
            }
        }
        else // option is not set
        {
            // check if option is required
            if (spec & option_spec::REQUIRED)
                throw required_option_missing("Option " + prepend_dash(short_id) + "/" + prepend_dash(long_id) +
                                              " is required but not set.");
        }
    }

    /*!\brief Handles command line flags, whether they are set or not.
     *
     * \param[out] value    The variable in which to store the given command line argument.
     * \param[in]  short_id The short identifier for the flag (e.g. `'i'`).
     * \param[in]  long_id  The long identifier for the flag (e.g. `"integer"`).
     *
     */
    void get_flag(bool & value,
                  char const short_id,
                  std::string const & long_id)
    {
        value = flag_is_set(short_id) || flag_is_set(long_id);
    }

    /*!\brief Handles command line positional option retrieval.
     *
     * \param[out] value     The variable in which to store the given command line argument.
     * \param[in]  validator The validator applied to the value after parsing (callable).
     *
     * \throws seqan3::parser_invalid_argument
     * \throws seqan3::too_few_arguments
     * \throws seqan3::validation_failed
     * \throws seqan3::parser_design_error
     *
     * \details
     *
     * This function assumes that
     * -#) argv has been stripped from all known options and flags
     * -#) argv has been checked for unknown options
     * -#) argv does not contain "--" anymore
     *  Thus we can simply iterate over non empty entries of argv.
     *
     * This function
     * - checks if the user did not provide enough arguments,
     * - retrieves the next(no container type) or all (container type),
     *   remaining non empty value/s in argv,
     * - re-throws the value cast exception with appended positional option information,
     * - and re-throws the validation exception with appended positional option information
     */
    template <typename option_type, typename validator_type>
    void get_positional_option(option_type & value,
                               validator_type && validator)
    {
        ++positional_option_count;
        auto it = std::find_if(argv.begin(), argv.end(), [](std::string const & s){return (s != "");});

        if (it == argv.end())
            throw too_few_arguments("Not enough positional arguments provided (Need at least " +
                                    std::to_string(positional_option_calls.size()) + "). See -h/--help for more information.");

        if (sequence_container_concept<option_type> && !std::is_same_v<option_type, std::string>) // vector/list will be filled with all remaining arguments
        {
            if (positional_option_count != (positional_option_calls.size()))
                throw parser_design_error("Lists are only allowed as the last positional option!");

            while (it != argv.end())
            {
                try
                {
                    retrieve_value(value, *it);
                }
                catch (parser_invalid_argument const & ex)
                {
                    throw parser_invalid_argument("Value cast failed for positional option " +
                                                  std::to_string(positional_option_count) + ": " + ex.what());
                }

                *it = ""; // remove arg from argv
                it = std::find_if(it, argv.end(), [](std::string const & s){return (s != "");});
                ++positional_option_count;
            }
        }
        else
        {
            try
            {
                retrieve_value(value, *it);
            }
            catch (parser_invalid_argument const & ex)
            {
                throw parser_invalid_argument("Value cast failed for positional option " +
                                              std::to_string(positional_option_count) + ": " + ex.what());
            }

            *it = ""; // remove arg from argv
        }

        try
        {
            validator(value);
        }
        catch (std::exception & ex)
        {
            throw validation_failed("Validation failed for positional option " +
                                    std::to_string(positional_option_count) + ": " + ex.what());
        }
    }

    //!\brief Stores get_option calls to be evaluated when calling format_parse::parse().
    std::vector<std::function<void()>> option_calls;
    //!\brief Stores get_flag calls to be evaluated when calling format_parse::parse().
    std::vector<std::function<void()>> flag_calls;
    //!\brief Stores get_positional_option calls to be evaluated when calling format_parse::parse().
    std::vector<std::function<void()>> positional_option_calls;
    //!\brief Keeps track of the number of specified positional options.
    unsigned positional_option_count{0};
    //!\brief Vector of command line arguments.
    std::vector<std::string> argv;
    //!\brief Number of command line arguments.
    int argc;
    //!\brief Artificial end of argv if -- was seen.
    std::vector<std::string>::iterator end_of_options_it;
};

} // namespace seqan3
