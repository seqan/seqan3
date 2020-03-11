// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides exceptions used in the I/O module.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <ios>
#include <stdexcept>

#include <seqan3/core/platform.hpp>

namespace seqan3
{

/*!\addtogroup io
 * \{
 */
// ----------------------------------------------------------------------------
// file open exceptions
// ----------------------------------------------------------------------------

//!\brief Thrown if there is no format that accepts a given file extension.
struct unhandled_extension_error : std::invalid_argument
{
    //!\brief Constructor that forwards the exception string.
    unhandled_extension_error(std::string const & s) : std::invalid_argument{s}
    {}
};

//!\brief Thrown if there is an unspecified filesystem or stream error while opening, e.g. permission problem.
struct file_open_error : std::runtime_error
{
    //!\brief Constructor that forwards the exception string.
    file_open_error(std::string const & s) : std::runtime_error{s}
    {}
};

//!\brief Thrown if there is a parse error, such as reading an unexpected character from an input stream.
struct parse_error : std::runtime_error
{
    //!\brief Constructor that forwards the exception string.
    parse_error(std::string const & s) : std::runtime_error{s}
    {}
};

//!\brief Thrown if there is an io error in low level io operations such as in std::basic_streambuf operations.
struct io_error : std::ios_base::failure
{
    //!\brief Constructor that forwards the exception string.
    io_error(std::string const & s,
             std::error_code const & ec) : std::ios_base::failure{s, ec}
    {}

    //!\brief Constructor that forwards the exception string.
    explicit io_error(std::string const & s) : io_error{s, std::error_code{std::io_errc::stream}}
    {}
};

// ----------------------------------------------------------------------------
// parse exceptions
// ----------------------------------------------------------------------------

//!\brief Thrown if I/O was expecting more input (e.g. a delimiter or a new line), but the end of input was reached.
struct unexpected_end_of_input : std::runtime_error
{
    //!\brief Constructor that forwards the exception string.
    unexpected_end_of_input(std::string const & s) : std::runtime_error{s}
    {}
};

// ----------------------------------------------------------------------------
// write exceptions
// ----------------------------------------------------------------------------

//!\brief Thrown if information given to output format didn't match expectations.
struct format_error : std::invalid_argument
{
    //!\brief Constructor that forwards the exception string.
    format_error(std::string const & s) : std::invalid_argument{s}
    {}
};

//!\}

} // namespace seqan3
