// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
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
 * \brief Contains parser related exceptions.
 */

#pragma once

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

//!\brief Argument parser exception thrown when encountering unkown option.
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

/*!\brief This exception is not an error but expected behavior that shall
 * terminate the program (e.g. when printing the help page).
 *
 * \details
 *
 * Behavior that triggers a parser interruption:
 *
 * - **--version** Prints the version information.
 * - **--copyright** Prints the copyright information.
 * - **-h/--help** Prints the help page (excluding advanced options).
 * - **-hh/--advanced-help** Prints the help page including advanced options.
 * - **--export-help [format]** Prints the help page information in the
 *                              given format (html/man/ctd).
 */
class parser_interruption : public std::exception
{
public:
    //!\brief Returns the error message.
    char const * what()
    {
        return error.c_str();
    }

private:
    //!\brief Message to the developer when exception is not caught.
    std::string error{std::string("ATTENTION: The parser printed or exported the help page/interface information.") +
                      "This behaviour is expected but the exception should be caught by the developer through " +
                      "a try-catch block (see documentation) and the program correctly terminated."};
};

} // namespace seqan3
