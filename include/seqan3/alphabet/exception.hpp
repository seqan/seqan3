// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Exceptions thrown by entities in the alphabet module.
 */

#pragma once

#include <concepts>
#include <stdexcept>
#include <string>

#include <seqan3/utility/char_operations/pretty_print.hpp>

namespace seqan3
{

/*!\brief An exception typically thrown by seqan3::alphabet::assign_char_strict.
 * \ingroup alphabet
 * \details
 * \experimentalapi{Experimental since version 3.1.}
 */
struct invalid_char_assignment : std::runtime_error
{
    /*!\brief Constructor that takes the type name and the failed character as arguments.
     * \details
     * \noapi
     */
    invalid_char_assignment(std::string const & type_name, std::string const & wrong_char) :
        std::runtime_error{std::string{"Assigning "} + wrong_char + " to an alphabet of type " + type_name
                           + " would incur information loss. If you want implicit conversion, use "
                             "seqan3::assign_char instead of seqan3::assign_char_strict."}
    {}

    /*!\overload
     *
     * \noapi
     */
    invalid_char_assignment(std::string const & type_name, char const wrong_char) :
        invalid_char_assignment{type_name, detail::make_printable(wrong_char)}
    {}

    /*!\overload
     *
     * \noapi
     */
    template <std::convertible_to<char> char_t>
    invalid_char_assignment(std::string const & type_name, char_t const wrong_char) :
        invalid_char_assignment{type_name, static_cast<char>(wrong_char)}
    {}
};

} // namespace seqan3
