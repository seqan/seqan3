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
 * \brief Checks if program is run interactively and retrieves dimensions of
 *        terminal (Transferred from seqan2).
 */

#pragma once

#ifndef _WIN32
#   include <sys/ioctl.h>
#endif

#include <unistd.h>

#include <cstdio>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

// ----------------------------------------------------------------------------
// Function is_terminal()
// ----------------------------------------------------------------------------

/*!\brief Check whether we are printing to a terminal.
 * \return True if code is run in a terminal, false otherwise.
 */
inline bool is_terminal()
{
#ifndef _WIN32
    return isatty(STDOUT_FILENO);
#else
    return false;
#endif
}

// ----------------------------------------------------------------------------
// Function get_terminal_size()
// ----------------------------------------------------------------------------

/*!\brief  Retrieve size of terminal.
 * \return The width of the current terminal in number of characters.
 *
 * \details
 *
 * Note: Only works on Linux/Unix.
 * TIOCGWINSZ is the command (number) to trigger filling the winsize struct.
 * STDOUT_FILENO is the default file descriptor (STDOUT_FILENO == fileno(stdout)).
 */
inline unsigned get_terminal_width()
{
#ifndef _WIN32

    struct winsize w;
    w.ws_row = 0;
    w.ws_col = 0;

    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);

    return w.ws_col;
#else
        return 80;    // ??
#endif
}

}  // namespace seqan::detail
