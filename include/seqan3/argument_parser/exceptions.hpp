// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Provides parser related exceptions.
 */

#pragma once

#include <stdexcept>

#include <seqan3/core/platform.hpp>

namespace seqan3
{
/*!\brief Argument parser exception that is thrown whenever there is an error
 * while parsing the command line arguments.
 * \ingroup argument_parser
 *
 * \details
 *
 * Errors caught by the argument_parser:
 *
 * - Unknown option/flag (not specified by developer but set by user)
 * - Too many positional options
 * - Too few positional options
 * - Option that was declared as required (option_spec::required) was not set
 * - Option is not a list but specified multiple times
 * - Type conversion failed
 * - Validation failed (as defined by the developer)
 *
 * \remark For a complete overview, take a look at \ref argument_parser
 */
class argument_parser_error : public std::runtime_error
{
public:
    /*!\brief The constructor.
     * \param[in] s The error message.
     */
    argument_parser_error(std::string const & s) : std::runtime_error(s)
    {}
};

//!\brief Argument parser exception thrown when encountering unknown option.
//!\ingroup argument_parser
//!\remark For a complete overview, take a look at \ref argument_parser
class unknown_option : public argument_parser_error
{
public:
    /*!\brief The constructor.
     * \param[in] s The error message.
     */
    unknown_option(std::string const & s) : argument_parser_error(s)
    {}
};

//!\brief Argument parser exception thrown when too many arguments are provided.
//!\ingroup argument_parser
//!\remark For a complete overview, take a look at \ref argument_parser
class too_many_arguments : public argument_parser_error
{
public:
    /*!\brief The constructor.
     * \param[in] s The error message.
     */
    too_many_arguments(std::string const & s) : argument_parser_error(s)
    {}
};

//!\brief Argument parser exception thrown when too few arguments are provided.
//!\ingroup argument_parser
//!\remark For a complete overview, take a look at \ref argument_parser
class too_few_arguments : public argument_parser_error
{
public:
    /*!\brief The constructor.
     * \param[in] s The error message.
     */
    too_few_arguments(std::string const & s) : argument_parser_error(s)
    {}
};

//!\brief Argument parser exception thrown when a required option is missing.
//!\ingroup argument_parser
//!\remark For a complete overview, take a look at \ref argument_parser
class required_option_missing : public argument_parser_error
{
public:
    /*!\brief The constructor.
     * \param[in] s The error message.
     */
    required_option_missing(std::string const & s) : argument_parser_error(s)
    {}
};

//!\brief Argument parser exception thrown when a non-list option is declared multiple times.
//!\ingroup argument_parser
//!\remark For a complete overview, take a look at \ref argument_parser
class option_declared_multiple_times : public argument_parser_error
{
public:
    /*!\brief The constructor.
     * \param[in] s The error message.
     */
    option_declared_multiple_times(std::string const & s) : argument_parser_error(s)
    {}
};

//!\brief Argument parser exception thrown when an incorrect argument is given as (positional) option value.
//!\ingroup argument_parser
//!\remark For a complete overview, take a look at \ref argument_parser
class user_input_error : public argument_parser_error
{
public:
    /*!\brief The constructor.
     * \param[in] s The error message.
     */
    user_input_error(std::string const & s) : argument_parser_error(s)
    {}
};

//!\brief Argument parser exception thrown when an argument could not be casted to the according type.
//!\ingroup argument_parser
//!\remark For a complete overview, take a look at \ref argument_parser
class validation_error : public argument_parser_error
{
public:
    /*!\brief The constructor.
     * \param[in] s The error message.
     */
    validation_error(std::string const & s) : argument_parser_error(s)
    {}
};

/*!\brief Argument parser exception that is thrown whenever there is an design
 * error directed at the developer of the application (e.g. Reuse of option).
 *
 * \details
 *
 * Errors caught by the argument_parser:
 *
 * - Reuse of a short or long identifier (must be unique)
 * - Both identifiers must not be empty (one is ok)
 */
class design_error : public argument_parser_error
{
public:
    /*!\brief The constructor.
     * \param[in] s The error message.
     */
    design_error(std::string const & s) : argument_parser_error(s)
    {}
};

} // namespace seqan3
