// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Includes customized exception types for the \link alignment alignment module \endlink.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <stdexcept>

#include <seqan3/core/platform.hpp>

namespace seqan3
{

//!\brief Thrown in function seqan3::erase_gap, if a position does not contain a gap.
class gap_erase_failure : public std::logic_error
{
public:
    /*!\brief The constructor.
     * \param[in] s The error message.
     */
    gap_erase_failure(std::string const & s) : std::logic_error{s}
    {}
};

//!\brief Thrown if the configuration of the alignment algorithm is invalid.
class invalid_alignment_configuration : public std::invalid_argument
{
public:
    /*!\brief The constructor.
     * \param[in] s The error message.
     */
    invalid_alignment_configuration(std::string const & s) : std::invalid_argument{s}
    {}
};

} // namespace seqan3
