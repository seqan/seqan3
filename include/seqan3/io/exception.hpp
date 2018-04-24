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
 * \brief Provides exceptions used in the I/O module.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <stdexcept>

#include <seqan3/core/platform.hpp>

namespace seqan3
{

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

//!\brief Thrown if there there was an unexpected end of input, e.g. the end-of-file marker was reached unexpectedly.
struct unexpected_end_of_error : parse_error
{
    //!\brief Constructor that forwards the exception string.
    unexpected_end_of_error(std::string const & s) : parse_error{s}
    {}
};

}
