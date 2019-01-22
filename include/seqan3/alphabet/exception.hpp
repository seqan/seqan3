// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Exceptions thrown by entities in the alphabet module.
 */

#pragma once

#include <stdexcept>
#include <string>

#include <seqan3/io/stream/parse_condition_detail.hpp>

namespace seqan3
{

//!\brief An exception typically thrown by seqan3::alphabet_concept::assign_char_strict.
struct invalid_char_assignment : std::runtime_error
{
    //!\brief Constructor that takes the type name and the failed character as arguments.
    invalid_char_assignment(std::string const & type_name, char const wrong_char) :
        std::runtime_error{std::string{"Assigning "} + detail::make_printable(wrong_char) + " to an alphabet of type " +
                           type_name + " would incur information loss. If you want implicit conversion, use "
                           "seqan3::assign_char instead of seqan3::assign_char_strict."}
    {}
};

} // namespace seqan3
