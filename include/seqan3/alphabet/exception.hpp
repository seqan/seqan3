// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Exceptions thrown by entities in the alphabet module.
 */

#pragma once

#include <stdexcept>
#include <string>

#include <seqan3/core/char_operations/pretty_print.hpp>
#include <seqan3/std/concepts>

namespace seqan3
{

//!\brief An exception typically thrown by seqan3::alphabet::assign_char_strict.
struct invalid_char_assignment : std::runtime_error
{
    //!\brief Constructor that takes the type name and the failed character as arguments.
    invalid_char_assignment(std::string const & type_name, std::string const & wrong_char) :
        std::runtime_error{std::string{"Assigning "} + wrong_char + " to an alphabet of type " +
                           type_name + " would incur information loss. If you want implicit conversion, use "
                           "seqan3::assign_char instead of seqan3::assign_char_strict."}
    {}

    //!\overload
    invalid_char_assignment(std::string const & type_name, char const wrong_char) :
        invalid_char_assignment{type_name, detail::make_printable(wrong_char)}
    {}

    //!\overload
    template <std::convertible_to<char> char_t>
    invalid_char_assignment(std::string const & type_name, char_t const wrong_char) :
        invalid_char_assignment{type_name, static_cast<char>(wrong_char)}
    {}
};

} // namespace seqan3
