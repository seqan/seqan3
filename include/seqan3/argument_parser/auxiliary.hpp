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
 * \brief Contains auxiliary information.
 */

#pragma once

#include <seqan3/io/concept.hpp>

namespace seqan3
{

/*!\brief Used to further specify argument_parser options/flags.
 * \ingroup argument_parser
 *
 * \details
 *
 * All options and flags are set to option_spec::DEFAULT unless specified
 * otherwise by the developer, e.g. when calling argument_parser::add_option().
 *
 * ```cpp
 *     myparser.add_option(myvar, 's', "special-op", "You know what you doin'?",
 *                         option_spec::ADVANCED);
 * ```
 */
enum option_spec
{
    DEFAULT  = 0, //!< The default were no checking or special displaying is happening.
    REQUIRED = 1, /*!< Set an option as required if you want to enforce that the user
                   * supplies this option when calling the program via the command line.
                   * If the option is missing, the argument_parser will automatically
                   * detect this and throw a invalid_argument exception.
                   */
    ADVANCED = 2, /* Set an option/flag to advanced if you do not want the option to
                   * be displayed in the normal help page (`-h/--help`). Instead, the
                   * advanced options are only displayed when calling `-hh/--advanced-help`
                   */
    HIDDEN   = 4, /* Set an option/flag to hidden, if you want to completely hide it from
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
    //!\brief The application name that will be displayed on the help page.
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
    //!\brief A link to github, your website or a wiki page.
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
    //!\brief Add lines of usage to the synopsis section of the help page (e.g. "[OPTIONS] FILE1 FILE1").
    std::vector<std::string> synopsis;
    /*!\brief Provide some examples on how to use your tool and what standard
     *        parameters might be appropriate in different cases (e.g.
     *        "./my_read_mapper -s 3 --my_flag path/infile1").
     */
    std::vector<std::string> examples;
};

namespace detail
{

/*!\brief Streams all parameters and returns a concatenated string.
 *
 * \tparam value_type Must be streamable (stream << value).
 * \param  values Variable number of parameters of any type that
 *                implement the stream operator.
 */
template <typename ...value_type>
//!\cond
requires (ostream_concept<std::iostream, std::remove_reference_t<value_type>> && ...)
//!\endcond
std::string to_string(value_type && ...values)
{
    std::stringstream stream;
    (stream << ... << std::forward<value_type>(values));
    return stream.str();
}

} // namespace detail

} // namespace seqan3
