// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::structure_file_output_options.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3
{

/*!\brief The options type defines various option members that influence the behaviour of all or some formats.
 * \ingroup structure_file
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
