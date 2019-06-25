// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
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

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/stream/concept.hpp>

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
