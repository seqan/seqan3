// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

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

// ----------------------------------------------------------------------------
// file open exceptions
// ----------------------------------------------------------------------------

//!\brief Thrown if there is no format that accepts a given file extension.
//!\ingroup io
struct unhandled_extension_error : std::invalid_argument
{
    //!\brief Constructor that forwards the exception string.
    unhandled_extension_error(std::string const & s) : std::invalid_argument{s}
    {}
};

//!\brief Thrown if there is an unspecified filesystem or stream error while opening, e.g. permission problem.
//!\ingroup io
struct file_open_error : std::runtime_error
{
    //!\brief Constructor that forwards the exception string.
    file_open_error(std::string const & s) : std::runtime_error{s}
    {}
};

//!\brief Thrown if there is a parse error, such as reading an unexpected character from an input stream.
//!\ingroup io
struct parse_error : std::runtime_error
{
    //!\brief Constructor that forwards the exception string.
    parse_error(std::string const & s) : std::runtime_error{s}
    {}
};

//!\brief Thrown if there is an io error in low level io operations such as in std::basic_streambuf operations.
//!\ingroup io
struct io_error : std::ios_base::failure
{
#if SEQAN3_WORKAROUND_GCC_NO_CXX11_ABI
    // std::ios_base::failure is missing the std::error_code constructor in pre-C++11 ABI
    // see https://en.cppreference.com/w/cpp/io/ios_base/failure
    using base_t = std::ios_base::failure;
    using base_t::base_t;
#else  // ^^^ workaround / no workaround vvv
    //!\brief Constructor that forwards the exception string.
    explicit io_error(std::string const & s, std::error_code const & ec = std::io_errc::stream) :
        std::ios_base::failure{s, ec}
    {}
#endif // SEQAN3_WORKAROUND_GCC_NO_CXX11_ABI
};

// ----------------------------------------------------------------------------
// parse exceptions
// ----------------------------------------------------------------------------

//!\brief Thrown if I/O was expecting more input (e.g. a delimiter or a new line), but the end of input was reached.
//!\ingroup io
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
//!\ingroup io
struct format_error : std::invalid_argument
{
    //!\brief Constructor that forwards the exception string.
    format_error(std::string const & s) : std::invalid_argument{s}
    {}
};

} // namespace seqan3
