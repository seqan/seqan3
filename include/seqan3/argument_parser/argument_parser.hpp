// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Provides seqan3::argument_parser class.
 */

#pragma once

#include <future>
#include <iostream>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <variant>
#include <vector>

// #include <seqan3/argument_parser/detail/format_ctd.hpp>
#include <seqan3/argument_parser/detail/format_help.hpp>
#include <seqan3/argument_parser/detail/format_html.hpp>
#include <seqan3/argument_parser/detail/format_man.hpp>
#include <seqan3/argument_parser/detail/format_parse.hpp>
#include <seqan3/argument_parser/detail/terminal.hpp>
#include <seqan3/argument_parser/detail/version_check.hpp>
#include <seqan3/core/debug_stream/detail/to_string.hpp>
#include <seqan3/core/detail/test_accessor.hpp>

namespace seqan3
{

/*!\brief The SeqAn command line parser.
 * \ingroup argument_parser
 *
 * \details
 *
 * The seqan3::argument_parser is a general purpose argument parser that provides
 * convenient access to the command line arguments passed to the program.
 * It automatically generates a help page and can export manual-pages as well
 * as HTML documentation.
 *
 * Furthermore common tool descriptor (CTD) files can be
 * exported (<a href="https://github.com/seqan/seqan3/issues/89">soon</a>)
 * and a server be queried for available application updates.
 *
 * ### Terminology
 *
 * Since the terms option and arguments are frequently used in different contexts
 * we want to first clarify our usage:
 *
 * - **options** [e.g. `-i myfile`, `--infile myfile`] refer to key-value pairs.
 *               The key is either a short indentifier, restricted to a single
 *               character `-i`, or a long identifier `--infile`.
 *
 * - **positional options** [e.g. `arg1`] refers to command line arguments that
 *                          are specified without an identifier/key, are always
 *                          required and are identified by their position.
 *
 * - **flags** [e.g. `-b`] refers to identifiers that are not followed by a
 *                         value (booleans) and therefore only indicate whether
 *                         they are set or not.
 *
 * ### Add/get options, flags or positional Options
 *
 * Adding an option is done in a single call. You simply
 * need to provide a predeclared variable and some additional information like
 * the identifier, description or advanced restrictions. To actually retrieve
 * the value from the command line and enable every other mechanism you need
 * to call the seqan3::argument_parser::parse function in the end.
 *
 * \include test/snippet/argument_parser/argument_parser_1.cpp
 *
 * Now you can call your application via the command line:
 *
 * ```console
 * MaxMuster% ./grade_avg_app -n Peter --bonus 1.0 2.0 1.7
 * Peter has an average grade of 1.425
 * ```
 * You can also display the help page automatically:
 *
 * ```console
 * MaxMuster% ./grade_avg_app --help
 * Grade-Average
 * =============
 *
 * POSITIONAL ARGUMENTS
 *     ARGUMENT-1 List of DOUBLE's
 *           Please specify your grades.
 *
 * OPTIONS
 *     -n, --name STRING
 *           Please specify your name.
 *     -b, --bonus
 *           Please specify if you got the bonus.
 *
 * VERSION
 *     Last update:
 *     Grade-Average version:
 *     SeqAn version: 3.0.0
 *
 * ```
 *
 * ### Errors that are caught by the argument_parser
 *
 * There are two different kinds of errors: Developer errors and user errors.
 *
 * Developer errors are those that violate the seqan3::argument_parser design
 * (e.g. calling the seqan3::argument_parser::parse function twice or specifying
 * two different options with the same identifier.)
 * In this case, a seqan3::design_error is thrown.
 *
 * The second kind are user errors, due to invalid command line calls. In this
 * case a seqan3::argument_parser_error is thrown.
 *
 * For example:
 *
 * ```console
 * MaxMuster% ./grade_avg_app -n Peter 2.0 abc 1.7
 * [PARSER ERROR] Value cast failed for positional option 2: Argument abc could not be casted to type DOUBLE.
 * ```
 *
 * See the seqan3::argument_parser::parse documentation for a detailed list of
 * which exceptions are caught.
 *
 * ### Update Notifications
 *
 * SeqAn applications that are using the seqan3::argument_parser can check SeqAn servers for version updates.
 * The functionality helps getting new versions out to users faster.
 * It is also used to inform application developers of new versions of the SeqAn library
 * which means that applications ship with less bugs.
 * For privacy implications, please see: https://github.com/seqan/seqan3/wiki/Update-Notifications.
 *
 * Users of applications that have this feature activated can opt-out, by either:
 *
 *  * disabling it for a specific application simply by setting the option `--version-check false/0` or
 *  * disabling it for all applications by setting the `SEQAN3_NO_VERSION_CHECK` environment variable.
 *
 * Note that in case there is no `--version-check` option (display available options with `-h/--help)`,
 * then the developer already disabled the version check functionality.
 *
 */
class argument_parser
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    argument_parser() = delete;                                    //!< Deleted.
    argument_parser(argument_parser const &) = delete;             //!< Deleted. Holds std::future.
    argument_parser & operator=(argument_parser const &) = delete; //!< Deleted. Holds std::future.
    argument_parser(argument_parser &&) = default;                 //!< Defaulted.
    argument_parser & operator=(argument_parser &&) = default;     //!< Defaulted.

    /*!\brief Initializes an seqan3::argument_parser object from the command line arguments.
     *
     * \param[in] app_name The name of the app that is displayed on the help page.
     * \param[in] argc The number of command line arguments.
     * \param[in] argv The command line arguments to parse.
     * \param[in] version_updates Notify users about version updates (default seqan3::update_notifications::on).
     * \param[in] subcommands A list of subcommands (see \link subcommand_arg_parse subcommand parsing \endlink).
     *
     * \throws seqan3::design_error if the application name contains illegal characters.
     *
     * The application name must only contain alpha-numeric characters, `_` or `-`,
     * i.e. the following regex must evaluate to true: `"^[a-zA-Z0-9_-]+$"`.
     *
     * See the [argument parser tutorial](https://docs.seqan.de/seqan3/main_dev/tutorial_argument_parser.html)
     * for more information about the version check functionality.
     */
    argument_parser(std::string const app_name,
                    int const argc,
                    char const * const * const argv,
                    update_notifications version_updates = update_notifications::on,
                    std::vector<std::string> subcommands = {}) :
        version_check_dev_decision{version_updates},
        subcommands{std::move(subcommands)}
    {
        if (!std::regex_match(app_name, app_name_regex))
        {
            throw design_error{("The application name must only contain alpha-numeric characters or '_' and '-' "
                                "(regex: \"^[a-zA-Z0-9_-]+$\").")};
        }

        for (auto & sub : this->subcommands)
        {
            if (!std::regex_match(sub, app_name_regex))
            {
                throw design_error{"The subcommand name must only contain alpha-numeric characters or '_' and '-' "
                                   "(regex: \"^[a-zA-Z0-9_-]+$\")."};
            }
        }

        info.app_name = std::move(app_name);

        init(argc, argv);
    }

    //!\brief The destructor.
    ~argument_parser()
    {
        // wait for another 3 seconds
        if (version_check_future.valid())
            version_check_future.wait_for(std::chrono::seconds(3));
    }
    //!\}

    /*!\name Adding options
     * \brief Add (positional) options and flags to the parser.
     * \{
     */
    /*!\brief Adds an option to the seqan3::argument_parser.
     *
     * \tparam option_type Must have a formatted input function (stream >> value).
     *                     If option_type is a container, its value type must have the
     *                     formatted input function (exception: std::string is not
     *                     regarded as a container).
     *                     See <a href="https://en.cppreference.com/w/cpp/named_req/FormattedInputFunction"> FormattedInputFunction </a>.
     * \tparam validator_type The type of validator to be applied to the option
     *                        value. Must model seqan3::validator.
     *
     * \param[in, out] value The variable in which to store the given command line argument.
     * \param[in] short_id The short identifier for the option (e.g. 'a').
     * \param[in] long_id The long identifier for the option (e.g. "age").
     * \param[in] desc The description of the option to be shown in the help page.
     * \param[in] spec Advanced option specification, see seqan3::option_spec.
     * \param[in] option_validator A seqan3::validator that verifies the value after parsing (callable).
     *
     * \throws seqan3::design_error
     */
    template <typename option_type, validator validator_type = detail::default_validator<option_type>>
        requires (argument_parser_compatible_option<option_type>
                  || argument_parser_compatible_option<std::ranges::range_value_t<option_type>>)
              && std::invocable<validator_type, option_type>
    void add_option(option_type & value,
                    char const short_id,
                    std::string const & long_id,
                    std::string const & desc,
                    option_spec const spec = option_spec::standard,
                    validator_type option_validator = validator_type{}) // copy to bind rvalues
    {
        if (sub_parser != nullptr)
            throw design_error{"You may only specify flags for the top-level parser."};

        verify_identifiers(short_id, long_id);
        // copy variables into the lambda because the calls are pushed to a stack
        // and the references would go out of scope.
        std::visit(
            [=, &value](auto & f)
            {
                f.add_option(value, short_id, long_id, desc, spec, option_validator);
            },
            format);
    }

    /*!\brief Adds a flag to the seqan3::argument_parser.
     *
     * \param[in, out] value     The variable which shows whether the flag is turned off (default) or on.
     * \param[in]      short_id  The short identifier for the flag (e.g. 'i').
     * \param[in]      long_id   The long identifier for the flag (e.g. "integer").
     * \param[in]      desc      The description of the flag to be shown in the help page.
     * \param[in]      spec      Advanced flag specification, see seqan3::option_spec.
     *
     * \throws seqan3::design_error if `value` is true.
     *
     */
    void add_flag(bool & value,
                  char const short_id,
                  std::string const & long_id,
                  std::string const & desc,
                  option_spec const spec = option_spec::standard)
    {
        if (value)
            throw design_error("A flag's default value must be false.");

        verify_identifiers(short_id, long_id);
        // copy variables into the lambda because the calls are pushed to a stack
        // and the references would go out of scope.
        std::visit(
            [=, &value](auto & f)
            {
                f.add_flag(value, short_id, long_id, desc, spec);
            },
            format);
    }

    /*!\brief Adds a positional option to the seqan3::argument_parser.
     *
     * \tparam option_type Must have a formatted input function (stream >> value).
     *                     If option_type is a container, its value type must have the
     *                     formatted input function (exception: std::string is not
     *                     regarded as a container).
     *                     See <a href="https://en.cppreference.com/w/cpp/named_req/FormattedInputFunction"> FormattedInputFunction </a>.
     * \tparam validator_type The type of validator to be applied to the option
     *                        value. Must model seqan3::validator.
     *
     * \param[in, out] value The variable in which to store the given command line argument.
     * \param[in] desc The description of the positional option to be shown in the help page.
     * \param[in] option_validator A seqan3::validator that verifies the value after parsing (callable).
     *
     * \throws seqan3::design_error
     *
     * \details
     *
     * The validator must be applicable to the given output variable (\p value).
     */
    template <typename option_type, validator validator_type = detail::default_validator<option_type>>
        requires (argument_parser_compatible_option<option_type>
                  || argument_parser_compatible_option<std::ranges::range_value_t<option_type>>)
              && std::invocable<validator_type, option_type>
    void add_positional_option(option_type & value,
                               std::string const & desc,
                               validator_type option_validator = validator_type{}) // copy to bind rvalues
    {
        if (sub_parser != nullptr)
            throw design_error{"You may only specify flags for the top-level parser."};

        if (has_positional_list_option)
            throw design_error{"You added a positional option with a list value before so you cannot add "
                               "any other positional options."};

        if constexpr (detail::is_container_option<option_type>)
            has_positional_list_option = true; // keep track of a list option because there must be only one!

        // copy variables into the lambda because the calls are pushed to a stack
        // and the references would go out of scope.
        std::visit(
            [=, &value](auto & f)
            {
                f.add_positional_option(value, desc, option_validator);
            },
            format);
    }
    //!\}

    /*!\brief Initiates the actual command line parsing.
     *
     * \attention The function must be called at the very end of all parser
     * related code and should be enclosed in a try catch block as the argument parser may throw.
     *
     * \throws seqan3::design_error if this function was already called before.
     *
     * \throws seqan3::option_declared_multiple_times if an option that is not a list was declared multiple times.
     * \throws seqan3::user_input_error if an incorrect argument is given as (positional) option value.
     * \throws seqan3::required_option_missing if the user did not provide a required option.
     * \throws seqan3::too_many_arguments if the command line call contained more arguments than expected.
     * \throws seqan3::too_few_arguments if the command line call contained less arguments than expected.
     * \throws seqan3::validation_error if the argument was not excepted by the provided validator.
     *
     * \details
     *
     * When no specific key words are supplied, the seqan3::argument_parser
     * starts to process the command line for specified options, flags and
     * positional options.
     *
     * If the given command line input (`argv`) contains the following keywords (in order of checking), the parser
     * will exit (std::exit) with error code 0 after doing the following:
     *
     * - **-h/\--help** Prints the help page.
     * - **-hh/\--advanced-help** Prints the help page including advanced options.
     * - <b>\--version</b> Prints the version information.
     * - <b>\--export-help [format]</b> Prints the application description in the given format (html/man/ctd).
     * - <b>\--version-check false/0/true/1</b> Disable/enable update notifications.
     *
     * Example:
     *
     * \include test/snippet/argument_parser/argument_parser_2.cpp
     *
     * The code above gives the following output when calling `--help`:
     *
     * ```console
     * MaxMuster$ ./age_app --help
     * The-Age-App
     * ===========
     *
     * OPTIONS
     *     -a, --user-age (signed 32 bit integer)
     *           Please specify your age.
     *
     * VERSION
     *     Last update:
     *     The-Age-App version:
     *     SeqAn version: 3.0.0
     *
     * Thanks for using The-Age-App!
     *
     * ```
     *
     * If the user does something wrong, it looks like this:
     * ```console
     * MaxMuster$ ./age_app --foo
     * The Age App - [PARSER ERROR] Unknown option --foo. Please see -h/--help for more information.
     * ```
     *
     * ```console
     * MaxMuster$ ./age_app -a abc
     * The Age App - [PARSER ERROR] Value cast failed for option -a: Argument abc
     *                              could not be casted to type (signed 32 bit integer).
     * ```
     */
    void parse()
    {
        if (parse_was_called)
            throw design_error("The function parse() must only be called once!");

        detail::version_checker app_version{info.app_name, info.version, info.url};

        if (std::holds_alternative<detail::format_parse>(format) && !subcommands.empty() && sub_parser == nullptr)
        {
            throw too_few_arguments{detail::to_string("You either forgot or misspelled the subcommand! Please specify"
                                                      " which sub-program you want to use: one of ",
                                                      subcommands,
                                                      ". Use -h/--help for more information.")};
        }

        if (app_version.decide_if_check_is_performed(version_check_dev_decision, version_check_user_decision))
        {
            // must be done before calling parse on the format because this might std::exit
            std::promise<bool> app_version_prom;
            version_check_future = app_version_prom.get_future();
            app_version(std::move(app_version_prom));
        }

        std::visit(
            [this](auto & f)
            {
                f.parse(info);
            },
            format);
        parse_was_called = true;
    }

    //!\brief Returns a reference to the sub-parser instance if
    //!       \link subcommand_arg_parse subcommand parsing \endlink was enabled.
    argument_parser & get_sub_parser()
    {
        if (sub_parser == nullptr)
        {
            throw design_error("No subcommand was provided at the construction of the argument parser!");
        }

        return *sub_parser;
    }

    /*!\brief Checks whether the option identifier (`id`) was set on the command line by the user.
     * \tparam id_type Either type `char` or a type that a `std::string` is constructible from.
     * \param[in] id The short (`char`) or long (`std::string`) option identifier to search for.
     * \returns `true` if option identifier `id` was set on the command line by the user.
     * \throws seqan3::design_error if the function is used incorrectly (see details below).
     *
     * \details
     *
     * You can only ask for option identifiers that were added to the parser beforehand via
     * `seqan3::argument_parser::add_option`.
     * As in the `seqan3::argument_parser::add_option` call, pass short identifiers as a `char` and long identifiers
     * as a `std::string` or a type that a `std::string` is constructible from (e.g. a `const char *`).
     *
     * ### Example
     *
     * \include test/snippet/argument_parser/is_option_set.cpp
     *
     * ### Exceptions
     *
     * This function throws a seqan3::design_error if
     * * `seqan3::argument_parser::parse()` was not called before.
     * * a long identifier was passed (e.g. a `std::string`) that only consists of a single character. If you mean to
     *   pass a short identifier, please pass it as a `char` not a `std::string`.
     * * the option identifier cannot be found in the list of valid option identifiers that were added to the parser
     *   via `seqan3::argument_parser::add_option()` calls beforehand.
     */
    template <typename id_type>
        requires std::same_as<id_type, char> || std::constructible_from<std::string, id_type>
    bool is_option_set(id_type const & id) const
    {
        if (!parse_was_called)
            throw design_error{"You can only ask which options have been set after calling the function `parse()`."};

        // the detail::format_parse::find_option_id call in the end expects either a char or std::string
        using char_or_string_t = std::conditional_t<std::same_as<id_type, char>, char, std::string>;
        char_or_string_t short_or_long_id = {id}; // e.g. convert char * to string here if necessary

        if constexpr (!std::same_as<id_type, char>) // long id was given
        {
            if (short_or_long_id.size() == 1)
            {
                throw design_error{"Long option identifiers must be longer than one character! If " + short_or_long_id
                                   + "' was meant to be a short identifier, please pass it as a char ('') not a string"
                                     " (\"\")!"};
            }
        }

        if (std::find(used_option_ids.begin(), used_option_ids.end(), std::string{id}) == used_option_ids.end())
            throw design_error{"You can only ask for option identifiers that you added with add_option() before."};

        // we only need to search for an option before the `end_of_options_indentifier` (`--`)
        auto end_of_options = std::find(cmd_arguments.begin(), cmd_arguments.end(), end_of_options_indentifier);
        auto option_it = detail::format_parse::find_option_id(cmd_arguments.begin(), end_of_options, short_or_long_id);
        return option_it != end_of_options;
    }

    //!\name Structuring the Help Page
    //!\{

    /*!\brief Adds an help page section to the seqan3::argument_parser.
     * \param[in] title The title of the section.
     * \param[in] spec Whether to always display this section title (seqan3::option_spec::standard), only when showing
     *                 the advanced help page (seqan3::option_spec::advanced) or never (seqan3::option_spec::hidden).
     * \details This only affects the help page and other output formats.
     */
    void add_section(std::string const & title, option_spec const spec = option_spec::standard)
    {
        std::visit(
            [&](auto & f)
            {
                f.add_section(title, spec);
            },
            format);
    }

    /*!\brief Adds an help page subsection to the seqan3::argument_parser.
     * \param[in] title The title of the subsection.
     * \param[in] spec Whether to always display this subsection title (seqan3::option_spec::standard), only when showing
     *                 the advanced help page (seqan3::option_spec::advanced) or never (seqan3::option_spec::hidden).
     * \details This only affects the help page and other output formats.
     */
    void add_subsection(std::string const & title, option_spec const spec = option_spec::standard)
    {
        std::visit(
            [&](auto & f)
            {
                f.add_subsection(title, spec);
            },
            format);
    }

    /*!\brief Adds an help page text line to the seqan3::argument_parser.
     * \param[in] text The text to print.
     * \param[in] is_paragraph Whether to insert as paragraph or just a line (Default: false).
     * \param[in] spec Whether to always display this line (seqan3::option_spec::standard), only when showing
     *                 the advanced help page (seqan3::option_spec::advanced) or never (seqan3::option_spec::hidden).
     * \details
     * If the line is not a paragraph (false), only one line break is appended, otherwise two line breaks are appended.
     * This only affects the help page and other output formats.
     */
    void add_line(std::string const & text, bool is_paragraph = false, option_spec const spec = option_spec::standard)
    {
        std::visit(
            [&](auto & f)
            {
                f.add_line(text, is_paragraph, spec);
            },
            format);
    }

    /*!\brief Adds an help page list item (key-value) to the seqan3::argument_parser.
     * \param[in] key  The key of the key-value pair of the list item.
     * \param[in] desc The value of the key-value pair of the list item.
     * \param[in] spec Whether to always display this list item (seqan3::option_spec::standard), only when showing
     *                 the advanced help page (seqan3::option_spec::advanced) or never (seqan3::option_spec::hidden).
     *
     * \details
     *
     * Note: This only affects the help page and other output formats.
     *
     * A list item is composed of a key (`key`) and value (`desc`)
     * and usually used for option identifier-description-pairs.
     * E.g.:
     *```console
     *     -a, --age LONG
     *            Super important integer for age.
     *```
     */
    void
    add_list_item(std::string const & key, std::string const & desc, option_spec const spec = option_spec::standard)
    {
        std::visit(
            [&](auto & f)
            {
                f.add_list_item(key, desc, spec);
            },
            format);
    }
    //!\}

    /*!\brief Aggregates all parser related meta data (see seqan3::argument_parser_meta_data struct).
     *
     * \attention You should supply as much information as possible to help users
     * of the application.
     *
     * \details
     *
     * The meta information is assembled in a struct to provide a central access
     * point that can be easily extended.
     *
     * You can access the members directly:
     * (see seqan3::argument_parser_meta_data for a list of the info members)
     *
     * \include test/snippet/argument_parser/argument_parser_3.cpp
     *
     * This will produce a nice help page when the user calls `-h` or `--help`:
     *
     * ```console
     * MaxMuster% ./penguin_app --help
     * Penguin_Parade - Organize your penguin parade
     * =============================================
     *
     * DESCRIPTION
     *     First Paragraph.
     *
     *     Second Paragraph.
     *
     * POSITIONAL ARGUMENTS
     *     ARGUMENT-1 List of STRING's
     *           Specify the names of the penguins.
     *
     * OPTIONS
     *     -d, --day (signed 32 bit integer)
     *           Please specify your preferred day.
     *     -m, --month (signed 32 bit integer)
     *           Please specify your preferred month.
     *     -y, --year (signed 32 bit integer)
     *           Please specify your preferred year.
     *
     * EXAMPLES
     *     ./penguin_parade Skipper Kowalski Rico Private -d 10 -m 02 -y 2017
     *
     * VERSION
     *     Last update: 12.01.2017
     *     Penguin_Parade version: 2.0.0
     *     SeqAn version: 3.0.0
     * ```
     */
    argument_parser_meta_data info;

private:
    //!\brief Keeps track of whether the parse function has been called already.
    bool parse_was_called{false};

    //!\brief Keeps track of whether the user has added a positional list option to check if this was the very last.
    bool has_positional_list_option{false};

    //!\brief Set on construction and indicates whether the developer deactivates the version check calls completely.
    update_notifications version_check_dev_decision{};

    //!\brief Whether the **user** specified to perform the version check (true) or not (false), default unset.
    std::optional<bool> version_check_user_decision;

    //!\brief Befriend seqan3::detail::test_accessor to grant access to version_check_future and format.
    friend struct ::seqan3::detail::test_accessor;

    //!\brief The future object that keeps track of the detached version check call thread.
    std::future<bool> version_check_future;

    //!\brief Validates the application name to ensure an escaped server call.
    std::regex app_name_regex{"^[a-zA-Z0-9_-]+$"};

    //!\brief Signals the argument parser that no options follow this string but only positional arguments.
    static constexpr std::string_view const end_of_options_indentifier{"--"};

    //!\brief Stores the sub-parser in case \link subcommand_arg_parse subcommand parsing \endlink is enabled.
    std::unique_ptr<argument_parser> sub_parser{nullptr};

    //!\brief Stores the sub-parser names in case \link subcommand_arg_parse subcommand parsing \endlink is enabled.
    std::vector<std::string> subcommands{};

    /*!\brief The format of the argument parser that decides the behavior when
     *        calling the seqan3::argument_parser::parse function.
     *
     * \details
     *
     * The format is set in the function argument_parser::init.
     */
    std::variant<detail::format_parse,
                 detail::format_help,
                 detail::format_short_help,
                 detail::format_version,
                 detail::format_html,
                 detail::format_man,
                 detail::format_copyright/*,
                 detail::format_ctd*/> format{detail::format_help{{}, false}}; // Will be overwritten in any case.

    //!\brief List of option/flag identifiers that are already used.
    std::set<std::string> used_option_ids{"h", "hh", "help", "advanced-help", "export-help", "version", "copyright"};

    //!\brief The command line arguments.
    std::vector<std::string> cmd_arguments{};

    /*!\brief Initializes the seqan3::argument_parser class on construction.
     *
     * \param[in] argc        The number of command line arguments.
     * \param[in] argv        The command line arguments.
     *
     * \throws seqan3::too_few_arguments if option --export-help was specified without a value
     * \throws seqan3::too_few_arguments if option --version-check was specified without a value
     * \throws seqan3::validation_error if the value passed to option --export-help was invalid.
     * \throws seqan3::validation_error if the value passed to option --version-check was invalid.
     * \throws seqan3::too_few_arguments if a sub parser was configured at construction but a subcommand is missing.
     *
     * \details
     *
     * This function adds all command line parameters to the cmd_arguments member variable
     * to take advantage of the vector functionality later on. Additionally,
     * the format member variable is set, depending on which parameters are given
     * by the user:
     *
     * - **no arguments** If no arguments are provided on the commandline, the
     *                    seqan3::detail::format_short_help is set.
     * - **-h/\--help** sets the format to seqan3::detail::format_help
     * - **-hh/\--advanced-help** sets the format to seqan3::detail::format_help
     *                           and show_advanced_options to `true`.
     * - <b>\--version</b> sets the format to seqan3::detail::format_version.
     * - <b>\--export-help html</b> sets the format to seqan3::detail::format_html.
     * - <b>\--export-help man</b> sets the format to seqan3::detail::format_man.
     * - <b>\--export-help ctd</b> sets the format to seqan3::detail::format_ctd.
     * - else the format is that to seqan3::detail::format_parse
     *
     * If `--export-help` is specified with a value other than html/man or ctd
     * an seqan3::argument_parser_error is thrown.
     */
    void init(int argc, char const * const * const argv)
    {
        if (argc <= 1) // no arguments provided
        {
            format = detail::format_short_help{};
            return;
        }

        bool special_format_was_set{false};

        for (int i = 1, argv_len = argc; i < argv_len; ++i) // start at 1 to skip binary name
        {
            std::string arg{argv[i]};

            if (std::ranges::find(subcommands, arg) != subcommands.end())
            {
                sub_parser = std::make_unique<argument_parser>(info.app_name + "-" + arg,
                                                               argc - i,
                                                               argv + i,
                                                               update_notifications::off);
                break;
            }

            if (arg == "-h" || arg == "--help")
            {
                format = detail::format_help{subcommands, false};
                init_standard_options();
                special_format_was_set = true;
            }
            else if (arg == "-hh" || arg == "--advanced-help")
            {
                format = detail::format_help{subcommands, true};
                init_standard_options();
                special_format_was_set = true;
            }
            else if (arg == "--version")
            {
                format = detail::format_version{};
                special_format_was_set = true;
            }
            else if (arg.substr(0, 13) == "--export-help") // --export-help=man is also allowed
            {
                std::string export_format;

                if (arg.size() > 13)
                {
                    export_format = arg.substr(14);
                }
                else
                {
                    if (argv_len <= i + 1)
                        throw too_few_arguments{"Option --export-help must be followed by a value."};
                    export_format = argv[i + 1];
                }

                if (export_format == "html")
                    format = detail::format_html{subcommands};
                else if (export_format == "man")
                    format = detail::format_man{subcommands};
                // TODO (smehringer) use when CTD support is available
                // else if (export_format == "ctd")
                //     format = detail::format_ctd{};
                else
                    throw validation_error{"Validation failed for option --export-help: "
                                           "Value must be one of [html, man]"};
                init_standard_options();
                special_format_was_set = true;
            }
            else if (arg == "--copyright")
            {
                format = detail::format_copyright{};
                special_format_was_set = true;
            }
            else if (arg == "--version-check")
            {
                if (++i >= argv_len)
                    throw too_few_arguments{"Option --version-check must be followed by a value."};

                arg = argv[i];

                if (arg == "1" || arg == "true")
                    version_check_user_decision = true;
                else if (arg == "0" || arg == "false")
                    version_check_user_decision = false;
                else
                    throw validation_error{"Value for option --version-check must be true (1) or false (0)."};

                // in case --version-check is specified it shall not be passed to format_parse()
                argc -= 2;
            }
            else
            {
                cmd_arguments.push_back(std::move(arg));
            }
        }

        if (!special_format_was_set)
            format = detail::format_parse(argc, cmd_arguments);
    }

    //!\brief Adds standard options to the help page.
    void init_standard_options()
    {
        add_subsection("Basic options:");
        add_list_item("\\fB-h\\fP, \\fB--help\\fP", "Prints the help page.");
        add_list_item("\\fB-hh\\fP, \\fB--advanced-help\\fP", "Prints the help page including advanced options.");
        add_list_item("\\fB--version\\fP", "Prints the version information.");
        add_list_item("\\fB--copyright\\fP", "Prints the copyright/license information.");
        add_list_item("\\fB--export-help\\fP (std::string)",
                      "Export the help page information. Value must be one of [html, man].");
        if (version_check_dev_decision == update_notifications::on)
            add_list_item("\\fB--version-check\\fP (bool)",
                          "Whether to check for the newest app version. Default: true.");
    }

    /*!\brief Checks whether the long identifier has already been used before.
    * \param[in] id The long identifier of the command line option/flag.
    * \returns `true` if an option or flag with the long identifier exists or `false`
    *          otherwise.
    */
    template <typename id_type>
    bool id_exists(id_type const & id)
    {
        if (detail::format_parse::is_empty_id(id))
            return false;
        return (!(used_option_ids.insert(std::string({id}))).second);
    }

    /*!\brief Verifies that the short and the long identifiers are correctly formatted.
     * \param[in] short_id The short identifier of the command line option/flag.
     * \param[in] long_id  The long identifier of the command line option/flag.
     * \throws seqan3::design_error
     * \details Specifically, checks that identifiers haven't been used before,
     *          the length of long IDs is either empty or longer than one char,
     *          the characters used in the identifiers are all valid,
     *          and at least one of short_id or long_id is given.
     */
    void verify_identifiers(char const short_id, std::string const & long_id)
    {
        constexpr auto allowed = is_alnum || is_char<'_'> || is_char<'@'>;

        if (id_exists(short_id))
            throw design_error("Option Identifier '" + std::string(1, short_id) + "' was already used before.");
        if (id_exists(long_id))
            throw design_error("Option Identifier '" + long_id + "' was already used before.");
        if (long_id.length() == 1)
            throw design_error("Long IDs must be either empty, or longer than one character.");
        if (!allowed(short_id) && !is_char<'\0'>(short_id))
            throw design_error("Option identifiers may only contain alphanumeric characters, '_', or '@'.");
        if (long_id.size() > 0 && is_char<'-'>(long_id[0]))
            throw design_error("First character of long ID cannot be '-'.");

        std::for_each(long_id.begin(),
                      long_id.end(),
                      [&allowed](char c)
                      {
                          if (!(allowed(c) || is_char<'-'>(c)))
                              throw design_error(
                                  "Long identifiers may only contain alphanumeric characters, '_', '-', or '@'.");
                      });
        if (detail::format_parse::is_empty_id(short_id) && detail::format_parse::is_empty_id(long_id))
            throw design_error("Option Identifiers cannot both be empty.");
    }
};

} // namespace seqan3
