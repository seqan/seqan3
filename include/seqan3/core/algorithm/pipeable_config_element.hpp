// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::pipeable_config_element.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#pragma once

#include <tuple>
#include <type_traits>

#include <seqan3/core/algorithm/concept.hpp>
#include <seqan3/core/type_traits/basic.hpp>

namespace seqan3
{

/*!\brief Adds pipe interface to configuration elements.
 * \ingroup algorithm
 * \tparam derived_t The type of the derived class.
 * \tparam value_t   The type of the wrapped value (defaults to void).
 */
template <typename derived_t>
struct pipeable_config_element<derived_t, void>
{
private:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    pipeable_config_element() = default; //!< Defaulted.
    pipeable_config_element(pipeable_config_element const &) = default; //!< Defaulted.
    pipeable_config_element(pipeable_config_element &&) = default; //!< Defaulted.
    pipeable_config_element & operator=(pipeable_config_element const &) = default; //!< Defaulted.
    pipeable_config_element & operator=(pipeable_config_element &&) = default; //!< Defaulted.
    ~pipeable_config_element() = default; //!< Defaulted.
    //!\}

    //!\brief Only the derived class is allowed to construct the element
    friend derived_t;
};

//! \overload
template <typename derived_t, typename value_t>
struct pipeable_config_element
{
    //!\brief The stored config value.
    value_t value{};
private:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    pipeable_config_element() = default; //!< Defaulted.
    pipeable_config_element(pipeable_config_element const &) = default; //!< Defaulted.
    pipeable_config_element(pipeable_config_element &&) = default; //!< Defaulted.
    pipeable_config_element & operator=(pipeable_config_element const &) = default; //!< Defaulted.
    pipeable_config_element & operator=(pipeable_config_element &&) = default; //!< Defaulted.
    ~pipeable_config_element() = default; //!< Defaulted.

    //!\brief Construction from the value type.
    constexpr pipeable_config_element(value_t const & val) : value(val) {}
    //!\}

    //!\brief Only the derived class is allowed to construct the element
    friend derived_t;
};

} // namespace seqan3
