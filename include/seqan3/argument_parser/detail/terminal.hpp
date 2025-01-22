// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \brief Checks if program is run interactively and retrieves dimensions of
 *        terminal (Transferred from seqan2).
 */

#pragma once

#ifndef _WIN32
#    include <sys/ioctl.h>
#else
#    include <windows.h>
#endif

#include <cstdio>
#include <unistd.h>

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
    return 80; // not implemented in windows
#endif
}

} // namespace seqan3::detail
