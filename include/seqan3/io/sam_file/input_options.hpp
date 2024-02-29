// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::sam_file_input_options.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <iostream>

#include <seqan3/core/platform.hpp>

namespace seqan3
{

/*!\brief The options type defines various option members that influence the behaviour of all or some formats.
 * \ingroup io_sam_file
 *
 * \remark For a complete overview, take a look at \ref io_sam_file
 */
template <typename sequence_legal_alphabet>
struct sam_file_input_options
{
    /*!\brief The stream to write warnings to. Defaults to std::cerr.
     * \details
     * ### Example
     * \include test/snippet/io/sam_file/sam_file_input_options.cpp
     * Output to std::cerr:
     * \include test/snippet/io/sam_file/sam_file_input_options.err
     * Output to std::cout:
     * \include test/snippet/io/sam_file/sam_file_input_options.out
     * \experimentalapi{Experimental since version 3.4.}
     */
    std::ostream * stream_warnings_to{std::addressof(std::cerr)};
};

} // namespace seqan3
