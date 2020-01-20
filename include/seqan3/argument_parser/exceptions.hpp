// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

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
 *
 * \details
 *
 * Errors caught by the argument_parser:
 *
 * - Unknown option/flag (not specified by developer but set by user)
 * - Too many positional options
 * - Too few positional options
 * - Option that was declared as required (option_spec::REQUIRED) was not set
 * - Option is not a list but specified multiple times
 * - Type conversion failed
 * - Validation failed (as defined by the developer)
 */
class parser_invalid_argument : public std::invalid_argument
{
public:
    /*!\brief The constructor.
     * \param[in] s The error message.
     */
    parser_invalid_argument(std::string const & s) : std::invalid_argument(s) {}
};

//!\brief Argument parser exception thrown when encountering unknown option.
class unknown_option : public parser_invalid_argument
{
public:
    /*!\brief The constructor.
     * \param[in] s The error message.
     */
    unknown_option(std::string const & s) : parser_invalid_argument(s) {}
};

//!\brief Argument parser exception thrown when too many arguments are provided.
class too_many_arguments : public parser_invalid_argument
{
public:
    /*!\brief The constructor.
     * \param[in] s The error message.
     */
    too_many_arguments(std::string const & s) : parser_invalid_argument(s) {}
};

//!\brief Argument parser exception thrown when too few arguments are provided.
class too_few_arguments : public parser_invalid_argument
{
public:
    /*!\brief The constructor.
     * \param[in] s The error message.
     */
    too_few_arguments(std::string const & s) : parser_invalid_argument(s) {}
};

//!\brief Argument parser exception thrown when a required option is missing.
class required_option_missing : public parser_invalid_argument
{
public:
    /*!\brief The constructor.
     * \param[in] s The error message.
     */
    required_option_missing(std::string const & s) : parser_invalid_argument(s) {}
};

//!\brief Argument parser exception thrown when a non-list option is declared multiple times.
class option_declared_multiple_times : public parser_invalid_argument
{
public:
    /*!\brief The constructor.
     * \param[in] s The error message.
     */
    option_declared_multiple_times(std::string const & s) : parser_invalid_argument(s) {}
};

//!\brief Argument parser exception thrown when an argument could not be casted to the according type.
class type_conversion_failed : public parser_invalid_argument
{
public:
    /*!\brief The constructor.
     * \param[in] s The error message.
     */
    type_conversion_failed(std::string const & s) : parser_invalid_argument(s) {}
};

//!\brief Argument parser exception thrown when an argument could not be casted to the according type.
class overflow_error_on_conversion : public parser_invalid_argument
{
public:
    /*!\brief The constructor.
     * \param[in] s The error message.
     */
    overflow_error_on_conversion(std::string const & s) : parser_invalid_argument(s) {}
};

//!\brief Argument parser exception thrown when an argument could not be casted to the according type.
class validation_failed : public parser_invalid_argument
{
public:
    /*!\brief The constructor.
     * \param[in] s The error message.
     */
    validation_failed(std::string const & s) : parser_invalid_argument(s) {}
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
class parser_design_error : public std::logic_error
{
public:
    /*!\brief The constructor.
     * \param[in] s The error message.
     */
    parser_design_error(std::string const & s) : std::logic_error(s) {}
};

} // namespace seqan3
