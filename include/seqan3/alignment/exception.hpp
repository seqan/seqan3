// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

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
    gap_erase_failure(std::string const & s) : std::logic_error{s} {}
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
