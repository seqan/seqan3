// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::sam_file_output_options.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3
{

/*!\brief The options type defines various option members that influence the behavior of all or some formats.
 * \ingroup io_sam_file
 *
 * \remark For a complete overview, take a look at \ref io_sam_file
 */
struct sam_file_output_options
{
    /*!\brief The default plain text line-ending is "\n", but on Windows an
     *        additional carriage return is recommended ("\r\n" for line-ending).
     */
    bool add_carriage_return = false;

    /*!\brief Whether to require a header for SAM files.
     *
     * \details
     *
     * In the official SAM format the header is optional but we highly
     * recommend to always specify the header nonetheless to be consistent with
     * BAM files (where the header is always required).
     * If you explicitly want the header not to be written and no related
     * checks to be done (e.g. the record reference name must be present in
     * the reference dictionary of the header) you may set this variable to
     * `false`.
     */
    bool sam_require_header = true;
};

} // namespace seqan3
