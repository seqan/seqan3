// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Provides the format_parse class.
 */

#pragma once

#include <sstream>
#include <string>
#include <vector>

#include <seqan3/argument_parser/detail/format_base.hpp>
#include <seqan3/core/char_operations/predicate.hpp>
#include <seqan3/core/detail/type_inspection.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/charconv>

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
     * \{
     */
    format_parse() = delete;                                     //!< Deleted.
    format_parse(format_parse const & pf) = default;             //!< Defaulted.
    format_parse & operator=(format_parse const & pf) = default; //!< Defaulted.
    format_parse(format_parse &&) = default;                     //!< Defaulted.
    format_parse & operator=(format_parse &&) = default;         //!< Defaulted.
    ~format_parse() = default;                                   //!< Defaulted.

    /*!\brief The constructor of the parse format.
     * \param[in] argc_ The number of command line arguments.
     * \param[in] argv_ The command line arguments to parse.
     */
    format_parse(int const argc_, std::vector<std::string> && argv_) :
        argc{argc_ - 1}, argv{std::move(argv_)}
    {}
    //!\}

    /*!\brief Adds an seqan3::detail::get_option call to be evaluated later on.
     * \copydetails seqan3::argument_parser::add_option
     */
    template <typename option_type, typename validator_type>
    void add_option(option_type & value,
                    char const short_id,
                    std::string const & long_id,
                    std::string const & SEQAN3_DOXYGEN_ONLY(desc),
                    option_spec const spec,
                    validator_type && validator)
    {
        option_calls.push_back([this, &value, short_id, long_id, spec, validator]()
        {
            get_option(value, short_id, long_id, spec, validator);
        });
    }

    /*!\brief Adds a get_flag call to be evaluated later on.
     * \copydetails seqan3::argument_parser::add_flag
     */
    void add_flag(bool & value,
                  char const short_id,
                  std::string const & long_id,
                  std::string const & SEQAN3_DOXYGEN_ONLY(desc),
                  option_spec const & SEQAN3_DOXYGEN_ONLY(spec))
    {
        flag_calls.push_back([this, &value, short_id, long_id]()
        {
            get_flag(value, short_id, long_id);
        });
    }

    /*!\brief Adds a get_positional_option call to be evaluated later on.
     * \copydetails seqan3::argument_parser::add_positional_option
     */
    template <typename option_type, typename validator_type>
    void add_positional_option(option_type & value,
                               std::string const & SEQAN3_DOXYGEN_ONLY(desc),
                               validator_type && validator)
    {
        positional_option_calls.push_back([this, &value, validator]()
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

        check_for_unknown_ids();

        if (end_of_options_it != argv.end())
            *end_of_options_it = ""; // remove -- before parsing positional arguments

        for (auto && f : positional_option_calls)
            f();

        check_for_left_over_args();
    }

    // functions are not needed for command line parsing but are part of the format interface.
    //!\cond
    void add_section(std::string const &, option_spec const) {}
    void add_subsection(std::string const &, option_spec const) {}
    void add_line(std::string const &, bool, option_spec const) {}
    void add_list_item(std::string const &, std::string const &, option_spec const) {}
    //!\endcond

    //!\brief Checks whether `id` is empty.
    template <typename id_type>
    static bool is_empty_id(id_type const & id)
    {
        if constexpr (std::same_as<remove_cvref_t<id_type>, std::string>)
            return id.empty();
        else // char
            return is_char<'\0'>(id);
    }

private:
    //!\brief Describes the result of parsing the user input string given the respective option value type.
    enum class option_parse_result
    {
        success, //!< Parsing of user input was successful.
        error, //!< There was some error while trying to parse the user input.
        overflow_error //!< Parsing was successful but the arithmetic value would cause an overflow.
    };

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

    /*!\brief Returns "-[short_id]/--[long_id]" if both are non-empty or just one of them if the other is empty.
    * \param[in] short_id The name of the short identifier.
    * \param[in] long_id  The name of the long identifier.
    * \returns The short_id prepended with a single dash and the long_id prepended with a double dash, separated by '/'.
    */
    std::string combine_option_names(char const short_id, std::string const & long_id)
    {
        if (short_id == '\0')
            return prepend_dash(long_id);
        else if (long_id.empty())
            return prepend_dash(short_id);
        else // both are set (note: both cannot be empty, this is caught before)
            return prepend_dash(short_id) + "/" + prepend_dash(long_id);
    }

    /*!\brief Finds the position of a short/long identifier in format_parse::argv.
     * \tparam id_type The identifier type; must be either of type `char` if it denotes a short identifier or
     *                 std::string if it denotes a long identifier.
     * \param[in] begin_it The iterator where to start the search of the identifier.
     *                     Note that the end iterator is kept as a member variable.
     * \param[in] id       The identifier to search for (must not contain dashes).
     * \returns An iterator pointing to the first occurrence of the identifier in
     *          the list pointed to by begin_t. If the list does not contain the
     *          identifier `id`, the member variable `end_of_options_it` is returned.
     *
     * \details
     *
     * **Valid short-id value pairs are: `-iValue`, `-i=Value`, or `-i Value`**
     * If the `id` passed to this function is of type `char`, it is assumed to be a short identifier.
     * The `id` is found by comparing the prefix of every argument in argv to the `id` prepended with a single `-`.
     *
     * **Valid long id value pairs are: `--id=Value`, `--id Value`**.
     * If the `id` passed to this function is of type `std::string`, it is assumed to be a long identifier.
     * The `id` is found by comparing every argument in argv to `id` prepended with two dashes (`--`)
     * or a prefix of such followed by the equal sign `=`.
     */
    template <typename id_type>
    std::vector<std::string>::iterator find_option_id(std::vector<std::string>::iterator const begin_it, id_type const & id)
    {
        if (is_empty_id(id))
            return end_of_options_it;

        return (std::find_if(begin_it, end_of_options_it,
            [&] (std::string const & current_arg)
            {
                std::string full_id = prepend_dash(id);

                if constexpr (std::same_as<id_type, char>) // short id
                {
                    // check if current_arg starts with "-o", i.e. it correctly identifies all short notations:
                    // "-ovalue", "-o=value", and "-o value".
                    return current_arg.substr(0, full_id.size()) == full_id;
                }
                else
                {
                    // only "--opt Value" or "--opt=Value" are valid
                    return current_arg.substr(0, full_id.size()) == full_id && // prefix is the same
                           (current_arg.size() == full_id.size() || current_arg[full_id.size()] == '='); // space or `=`
                }
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

    /*!\brief Tries to parse an input string into a value using the stream `operator>>`.
     * \tparam option_t Must model seqan3::input_stream_over.
     * \param[out] value Stores the parsed value.
     * \param[in] in The input argument to be parsed.
     * \returns seqan3::option_parse_result::error if `in` could not be parsed via the stream
     *          operator and otherwise seqan3::option_parse_result::success.
     */
    template <typename option_t>
    //!\cond
        requires input_stream_over<std::istringstream, option_t>
    //!\endcond
    option_parse_result parse_option_value(option_t & value, std::string const & in)
    {
        std::istringstream stream{in};
        stream >> value;

        if (stream.fail() || !stream.eof())
            return option_parse_result::error;

        return option_parse_result::success;
    }

    /*!\brief Sets an option value depending on the keys found in seqan3::enumeration_names<option_t>.
     * \tparam option_t Must model seqan3::named_enumeration.
     * \param[out] value Stores the parsed value.
     * \param[in] in The input argument to be parsed.
     * \returns seqan3::option_parse_result::error if `in` could not be found in the
     *          seqan3::enumeration_names<option_t> map and otherwise seqan3::option_parse_result::success.
     */
    template <named_enumeration option_t>
    option_parse_result parse_option_value(option_t & value, std::string_view const in)
    {
        auto map = seqan3::enumeration_names<option_t>;

        if (auto it = map.find(in); it == map.end())
            return option_parse_result::error;
        else
            value = it->second;

        return option_parse_result::success;
    }

    //!\cond
    option_parse_result parse_option_value(std::string & value, std::string const & in)
    {
        value = in;
        return option_parse_result::success;
    }
    //!\endcond

    /*!\brief Parses the given option value and appends it to the target container.
     * \tparam container_option_t Must model the seqan3::sequence_container and
     *                            its value_type must model the seqan3::input_stream_over
     *
     * \param[out] value The container that stores the parsed value.
     * \param[in] in The input argument to be parsed.
     * \returns A seqan3::option_parse_result whether parsing was successful or not.
     */
    template <sequence_container container_option_t>
    //!\cond
        requires input_stream_over<std::istringstream, typename container_option_t::value_type>
    //!\endcond
    option_parse_result parse_option_value(container_option_t & value, std::string const & in)
    {
        typename container_option_t::value_type tmp{};

        auto res = parse_option_value(tmp, in);

        if (res == option_parse_result::success)
            value.push_back(tmp);

        return res;
    }

    /*!\brief Tries to parse an input string into an arithmetic value.
     * \tparam option_t The option value type; must model seqan3::arithmetic.
     * \param[out] value Stores the parsed value.
     * \param[in] in The input argument to be parsed.
     * \returns seqan3::option_parse_result::error if `in` could not be parsed to an arithmetic type
     *          via std::from_chars, seqan3::option_parse_result::overflow_error if `in` could be parsed but the
     *          value is too large for the respective type, and otherwise seqan3::option_parse_result::success.
     *
     * \details
     *
     * This function delegates to std::from_chars.
     */
    template <arithmetic option_t>
    //!\cond
        requires input_stream_over<std::istringstream, option_t>
    //!\endcond
    option_parse_result parse_option_value(option_t & value, std::string const & in)
    {
        auto res = std::from_chars(&in[0], &in[in.size()], value);

        if (res.ec == std::errc::result_out_of_range)
            return option_parse_result::overflow_error;
        else if (res.ec == std::errc::invalid_argument || res.ptr != &in[in.size()])
            return option_parse_result::error;

        return option_parse_result::success;
    }

    /*!\brief Tries to parse an input string into a boolean value.
     * \param[out] value Stores the parsed value.
     * \param[in] in The input argument to be parsed.
     * \returns A seqan3::option_parse_result whether parsing was successful or not.
     *
     * \details
     *
     * This function accepts the strings "0" or "false" which sets sets `value` to `false` or "1" or "true" which
     * sets `value` to `true`.
     */
    option_parse_result parse_option_value(bool & value, std::string const & in)
    {
        if (in == "0")
            value = false;
        else if (in == "1")
            value = true;
        else if (in == "true")
            value = true;
        else if (in == "false")
            value = false;
        else
            return option_parse_result::error;

        return option_parse_result::success;
    }

    /*!\brief Tries to parse an input string into boolean value.
     * \param[in] res A result value of parsing an input string to the respective option value type.
     * \param[in] option_name The name of the option whose input was parsed.
     * \param[in] input_value The original user input in question.
     *
     * \throws seqan3::user_input_error if `res` was not seqan3::option_parse_result::success.
     */
    template <typename option_type>
    void throw_on_input_error(option_parse_result const res,
                              std::string const & option_name,
                              std::string const & input_value)
    {
        std::string msg{"Value parse failed for " + option_name + ": "};

        if (res == option_parse_result::error)
        {
            throw user_input_error{msg + "Argument " + input_value + " could not be parsed as type " +
                                   get_type_name_as_string(input_value) + "."};
        }

        if constexpr (arithmetic<option_type>)
        {
            if (res == option_parse_result::overflow_error)
            {
                throw user_input_error{msg + "Numeric argument " + input_value + " is not in the valid range [" +
                                       std::to_string(std::numeric_limits<option_type>::min()) + "," +
                                       std::to_string(std::numeric_limits<option_type>::max()) + "]."};
            }
        }

        assert(res == option_parse_result::success); // if nothing was thrown, the result must have been a success
    }

    /*!\brief Handles value retrieval for options based on different kev value pairs.
     *
     * \param[out] value     Stores the value found in argv, parsed by parse_option_value.
     * \param[in]  option_it The iterator where the option identifier was found.
     * \param[in]  id        The option identifier supplied on the command line.
     *
     * \throws seqan3::too_few_arguments if the option was not followed by a value.
     * \throws seqan3::user_input_error if the given option value was invalid.
     *
     * \details
     *
     * The value at option_it is inspected whether it is an '-key value', '-key=value'
     * or '-keyValue' pair and the input is extracted accordingly. The input
     * will then be tried to be parsed into the `value` parameter.
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
                        throw too_few_arguments("Missing value for option " + prepend_dash(id));
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
                    throw too_few_arguments("Missing value for option " + prepend_dash(id));
                input_value = *option_it;
                *option_it = ""; // remove value
            }

            auto res = parse_option_value(value, input_value);
            throw_on_input_error<option_type>(res, prepend_dash(id), input_value);

            return true;
        }
        return false;
    }

    /*!\brief Handles value retrieval (non container type) options.
     *
     * \param[out] value Stores the value found in argv, parsed by parse_option_value.
     * \param[in] id The option identifier supplied on the command line.
     *
     * \throws seqan3::option_declared_multiple_times
     *
     * \details
     *
     * If the option identifier is found in format_parse::argv, the value of
     * the following position in argv is tried to be parsed given the respective option value type
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
                                                " is no list/container but declared multiple times.");

       return (it != end_of_options_it); // first search was successful or not
    }

    /*!\brief Handles value retrieval (container type) options.
     *
     * \param[out] value Stores all values found in argv, parsed by parse_option_value.
     * \param[in]  id    The option identifier supplied on the command line.
     *
     * \details
     *
     * Since option_type is a container, the option is a list and can be parsed
     * multiple times.
     *
     */
    template <sequence_container option_type, typename id_type>
    //!\cond
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
    void check_for_unknown_ids()
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
     * \param[in]  short_id  The short identifier for the option (e.g. 'i').
     * \param[in]  long_id   The long identifier for the option (e.g. "integer").
     * \param[in]  spec      Advanced option specification, see seqan3::option_spec.
     * \param[in]  validator The validator applied to the value after parsing (callable).
     *
     * \throws seqan3::option_declared_multiple_times
     * \throws seqan3::validation_error
     * \throws seqan3::required_option_missing
     *
     * \details
     *
     * This function
     * - checks if the option is required but not set,
     * - retrieves any value found by the short or long identifier,
     * - throws on (mis)use of both identifiers for non-container type values,
     * - re-throws the validation exception with appended option information.
     */
    template <typename option_type, typename validator_type>
    void get_option(option_type & value,
                     char const short_id,
                     std::string const & long_id,
                     option_spec const spec,
                     validator_type && validator)
    {
        bool short_id_is_set{get_option_by_id(value, short_id)};
        bool long_id_is_set{get_option_by_id(value, long_id)};

        // if value is no container we need to check for multiple declarations
        if (short_id_is_set && long_id_is_set &&
            !(sequence_container<option_type> && !std::is_same_v<option_type, std::string>))
            throw option_declared_multiple_times("Option " + combine_option_names(short_id, long_id) +
                                                 " is no list/container but specified multiple times");

        if (short_id_is_set || long_id_is_set)
        {
            try
            {
                validator(value);
            }
            catch (std::exception & ex)
            {
                throw validation_error(std::string("Validation failed for option ") +
                                        combine_option_names(short_id, long_id) + ": " + ex.what());
            }
        }
        else // option is not set
        {
            // check if option is required
            if (spec & option_spec::REQUIRED)
                throw required_option_missing("Option " + combine_option_names(short_id, long_id) +
                                              " is required but not set.");
        }
    }

    /*!\brief Handles command line flags, whether they are set or not.
     *
     * \param[out] value    The variable in which to store the given command line argument.
     * \param[in]  short_id The short identifier for the flag (e.g. 'i').
     * \param[in]  long_id  The long identifier for the flag (e.g. "integer").
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
     * \throws seqan3::argument_parser_error
     * \throws seqan3::too_few_arguments
     * \throws seqan3::validation_error
     * \throws seqan3::design_error
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
     * - retrieves the next (no container type) or all (container type) remaining non empty value/s in argv
     */
    template <typename option_type, typename validator_type>
    void get_positional_option(option_type & value,
                               validator_type && validator)
    {
        ++positional_option_count;
        auto it = std::find_if(argv.begin(), argv.end(), [](std::string const & s){return (s != "");});

        if (it == argv.end())
            throw too_few_arguments("Not enough positional arguments provided (Need at least " +
                                    std::to_string(positional_option_calls.size()) +
                                    "). See -h/--help for more information.");

        if (sequence_container<option_type> && !std::is_same_v<option_type, std::string>) // vector/list will be filled with all remaining arguments
        {
            assert(positional_option_count == positional_option_calls.size()); // checked on set up.

            while (it != argv.end())
            {
                auto res = parse_option_value(value, *it);
                std::string id = "positional option" + std::to_string(positional_option_count);
                throw_on_input_error<option_type>(res, id, *it);

                *it = ""; // remove arg from argv
                it = std::find_if(it, argv.end(), [](std::string const & s){return (s != "");});
                ++positional_option_count;
            }
        }
        else
        {
            auto res = parse_option_value(value, *it);
            std::string id = "positional option" + std::to_string(positional_option_count);
            throw_on_input_error<option_type>(res, id, *it);

            *it = ""; // remove arg from argv
        }

        try
        {
            validator(value);
        }
        catch (std::exception & ex)
        {
            throw validation_error("Validation failed for positional option " +
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
    //!\brief Number of command line arguments.
    int argc;
    //!\brief Vector of command line arguments.
    std::vector<std::string> argv;
    //!\brief Artificial end of argv if \-- was seen.
    std::vector<std::string>::iterator end_of_options_it;
};

} // namespace seqan3
