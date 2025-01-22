// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::structure_file_output_options.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3
{

/*!\brief The options type defines various option members that influence the behaviour of all or some formats.
 * \ingroup io_structure_file
 *
 * \remark For a complete overview, take a look at \ref io_structure_file
 */
struct structure_file_output_options
{
    /*!\brief The default plain text line-ending is "\n", but on Windows an additional carriage return is
     *        recommended ("\r\n" for line-ending).
     */
    bool add_carriage_return = false;

    //!\brief The precision for writing floating point types.
    int precision = 6;
};

} // namespace seqan3
