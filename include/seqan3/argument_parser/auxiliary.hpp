// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Provides auxiliary information.
 */

#pragma once

#include <sstream>
#include <vector>

#include <seqan3/core/detail/customisation_point.hpp>
#include <seqan3/core/detail/debug_stream_type.hpp>
#include <seqan3/core/type_traits/basic.hpp>
#include <seqan3/io/stream/concept.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/type_traits>

namespace seqan3::custom
{

/*!\brief A type that can be specialised to provide customisation point implementations for the seqan3::argument_parser
 *        such that third party types may be adapted.
 * \tparam t The type you wish to specialise for.
 * \ingroup argument_parser
 *
 * \details
 *
 * ### Named Enumerations
 *
 * In order to use a third party type within the seqan3::argument_parser::add_option or
 * seqan3::argument_parser::add_positional_option call, you can specialise this struct in the following way:
 *
 * \include test/snippet/argument_parser/custom_argument_parsing_enumeration.cpp
 *
 * Please note that by default the `t const`, `t &` and `t const &` specialisations of this class inherit the
 * specialisation for `t` so you usually only need to provide a specialisation for `t`.
 *
 * \note Only use this if you cannot provide respective functions in your namespace. See the tutorial
 * \ref tutorial_argument_parser for an example of customising a type within your own namespace.
 */
template <typename t>
struct argument_parsing
{}; // forward

//!\cond
template <typename t>
struct argument_parsing<t const> : argument_parsing<t>
{};

template <typename t>
struct argument_parsing<t &> : argument_parsing<t>
{};

template <typename t>
struct argument_parsing<t const &> : argument_parsing<t>
{};
//!\endcond

} // seqan3::custom

namespace seqan3::detail::adl_only
{

//!\brief Poison-pill overload to prevent non-ADL forms of unqualified lookup.
template <typename t>
std::unordered_map<std::string_view, t> enumeration_names(t) = delete;

/*!\brief Functor definition for seqan3::enumeration_names.
 *
 * We need a class template here because we need the original option_t next to s_option_t that may be wrapped
 * std::type_identity to be default constructible. We need option_t to be default constructible because the
 * respective CPO should be callable without a function parameter.
 */
template <typename option_t>
struct enumeration_names_fn
{
    //!\brief `option_t` with cvref removed and possibly wrapped in std::type_identity.
    using s_option_t = std::conditional_t<std::is_nothrow_default_constructible_v<remove_cvref_t<option_t>> &&
                                          seqan3::is_constexpr_default_constructible_v<remove_cvref_t<option_t>>,
                                          remove_cvref_t<option_t>,
                                          std::type_identity<option_t>>;

    SEQAN3_CPO_IMPL(1, (deferred_type_t<seqan3::custom::argument_parsing<option_t>, decltype(v)>::enumeration_names))
    SEQAN3_CPO_IMPL(0, (enumeration_names(v))) // ADL

    //!\brief Operator definition.
    template <typename dummy = int> // need to make this a template to enforce deferred initialisation
    //!\cond
        requires requires
        {
            { impl(priority_tag<1>{}, s_option_t{}, dummy{}) };
            std::same_as<decltype(impl(priority_tag<1>{}, s_option_t{}, dummy{})),
                         std::unordered_map<std::string_view, remove_cvref_t<option_t>>>;
        }
    //!\endcond
    auto operator()() const
    {
        return impl(priority_tag<1>{}, s_option_t{});
    }
};

} // namespace seqan3::detail::adl_only

namespace seqan3
{

/*!\name Customisation Points
 * \{
 */

/*!\brief Return a conversion map from std::string_view to option_type.
 * \tparam your_type Type of the value to retrieve the conversion map for.
 * \param value The value is not used, just its type.
 * \returns A std::unordered_map<std::string_view, your_type> that maps a string identifier to a value of your_type.
 * \ingroup argument_parser
 * \details
 *
 * This is a function object. Invoke it with the parameter(s) specified above.
 *
 * It acts as a wrapper and looks for two possible implementations (in this order):
 *
 *   1. A static member `enumeration_names` in `seqan3::custom::argument_parsing<your_type>` that is of type
 *      `std::unordered_map<std::string_view, your_type>>`.
 *   2. A free function `enumeration_names(your_type const a)` in the namespace of your type (or as `friend`) which
 *      returns a `std::unordered_map<std::string_view, your_type>>`.
 *
 * ### Example
 *
 * If you are working on a type in your own namespace, you should implement a free function like this:
 *
 * \include test/snippet/argument_parser/custom_enumeration.cpp
 *
 * **Only if you cannot access the namespace of your type to customize** you may specialize
 * the seqan3::custom::argument_parsing struct like this:
 *
 * \include test/snippet/argument_parser/custom_argument_parsing_enumeration.cpp
 *
 * ### Customisation point
 *
 * This is a customisation point (see \ref about_customisation). To specify the behaviour for your own type,
 * simply provide one of the two functions specified above.
 */
template <typename option_type>
//!\cond
    requires requires { { detail::adl_only::enumeration_names_fn<option_type>{}() }; }
//!\endcond
inline auto const enumeration_names = detail::adl_only::enumeration_names_fn<option_type>{}();
//!\}

/*!\interface seqan3::named_enumeration <>
 * \brief Checks whether the free function seqan3::enumeration_names can be called on the type.
 * \ingroup argument_parser
 * \tparam option_type The type to check.
 *
 * ### Requirements
 *
 * * A instance of seqan3::enumeration_names<option_type> must exist and be of type
 *   `std::unordered_map<std::string, option_type>`.
 */
//!\cond
template <typename option_type>
SEQAN3_CONCEPT named_enumeration = requires
{
    { seqan3::enumeration_names<option_type> };
};
//!\endcond

/*!\interface seqan3::argument_parser_compatible_option <>
 * \brief Checks whether the the type can be used in an add_(positional_)option call on the argument parser.
 * \ingroup argument_parser
 * \tparam option_type The type to check.
 *
 * ### Requirements
 *
 * In order to model this concept, the type must either be streamable to std::istringstream or
 * model seqan3::named_enumeration<option_type>.
 */
//!\cond
template <typename option_type>
SEQAN3_CONCEPT argument_parser_compatible_option = input_stream_over<std::istringstream, option_type> ||
                                                   named_enumeration<option_type>;
//!\endcond

/*!\name Formatted output overloads
 * \{
 */
/*!\brief A type (e.g. an enum) can be made debug streamable by customizing the seqan3::enumeration_names.
 * \tparam option_type Type of the enum to be printed.
 * \param s  The seqan3::debug_stream.
 * \param op The value to print.
 * \relates seqan3::debug_stream_type
 *
 * \details
 *
 * This searches the seqan3::enumeration_names of the respective type for the value \p op and prints the
 * respective string if found or '\<UNKNOWN_VALUE\>' if the value cannot be found in the map.
 */
template <typename char_t, typename option_type>
//!\cond
    requires named_enumeration<remove_cvref_t<option_type>>
//!\endcond
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & s, option_type && op)
{
    for (auto & [key, value] : enumeration_names<option_type>)
    {
        if (op == value)
            return s << key;
    }

    return s << "<UNKNOWN_VALUE>";
}
//!\}

/*!\brief Used to further specify argument_parser options/flags.
 * \ingroup argument_parser
 *
 * \details
 *
 * All options and flags are set to option_spec::DEFAULT unless specified
 * otherwise by the developer, e.g. when calling argument_parser::add_option().
 *
 * \include test/snippet/argument_parser/auxiliary.cpp
 */
enum option_spec
{
    DEFAULT  = 0, //!< The default were no checking or special displaying is happening.
    REQUIRED = 1, /*!< Set an option as required if you want to enforce that the user
                   * supplies this option when calling the program via the command line.
                   * If the option is missing, the argument_parser will automatically
                   * detect this and throw a invalid_argument exception.
                   */
    ADVANCED = 2, /*!< Set an option/flag to advanced if you do not want the option to
                   * be displayed in the normal help page (`-h/--help`). Instead, the
                   * advanced options are only displayed when calling `-hh/--advanced-help`
                   */
    HIDDEN   = 4, /*!< Set an option/flag to hidden, if you want to completely hide it from
                   * the user. It will never appear on the help page nor any export format.
                   * For example, this can be useful for debugging reasons.
                   * (e.g. "A tool for mapping reads to the genome").
                   */
};

/*!\brief Stores all parser related meta information of the seqan3::argument_parser.
 * \ingroup argument_parser
 *
 * \attention You should supply as much information as possible to help the users
 * of your application.
 *
 * \details
 *
 * The meta information is assembled in a struct to provide a central access
 * point that can be easily extended.
 */
struct argument_parser_meta_data // holds all meta information
{
    /*!\brief The application name that will be displayed on the help page.
     *
     * The application name must only contain alpha-numeric characters, '_' or '-',
     * i.e. the following regex must evaluate to true: `\"^[a-zA-Z0-9_-]+$\"`.
     */
    std::string app_name;
    //!\brief The version information `MAJOR.MINOR.PATH` (e.g. 3.1.3)
    std::string version;
    //!\brief A short description of the application (e.g. "A tool for mapping reads to the genome").
    std::string short_description;
    //!\brief Your name ;-)
    std::string author;
    //!\brief The author's e-mail address for correspondence.
    std::string email;
    /*!\brief The date that the application was last updated. Keep this updated,
     *!          since it will tell your users that the application is maintained.
     */
    std::string date;
    //!\brief A link to  your github/gitlab project with the newest release.
    std::string url;
    //!\brief Brief copyright (and/or license) information.
    std::string short_copyright;
    /*!\brief Detailed copyright information that will be displayed
     *        when the user specifies "--copyright" on the command line.
     */
    std::string long_copyright;
    //!\brief How  users shall cite your application.
    std::string citation;
    /*!\brief The title of your man page when exported by specifying
     *        "--export-help man" on the common line.
     */
    std::string man_page_title;
    //!\brief The man page section info (type `man man` on the command line for more information).
    unsigned man_page_section{1};
    /*!\brief A more detailed description that is displayed on the help
     *        page in the section "DESCRIPTION". Each `std::string` appended
     *        to the description vector will be treated as a paragraph and
     *        is separated by a new line.
     */
    std::vector<std::string> description;
    /*!\brief Add lines of usage to the synopsis section of the help page (e.g.
     *        "./my_read_mapper [OPTIONS] FILE1 FILE1").
     */
    std::vector<std::string> synopsis;
    /*!\brief Provide some examples on how to use your tool and what standard
     *        parameters might be appropriate in different cases (e.g.
     *        "./my_read_mapper -s 3 --my_flag path/infile1").
     */
    std::vector<std::string> examples;
};

} // namespace seqan3
