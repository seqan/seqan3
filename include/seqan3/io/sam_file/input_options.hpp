// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::sam_file_input_options.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3
{

/*!\brief The options type defines various option members that influence the behaviour of all or some formats.
 * \ingroup io_sam_file
 *
 * \note As of now, there are no specific options for the SAM format. This class may be used in the future for possible
 *       SAM parsing extensions.
 * \remark For a complete overview, take a look at \ref io_sam_file
 */
template <typename sequence_legal_alphabet>
struct sam_file_input_options
{};

} // namespace seqan3
