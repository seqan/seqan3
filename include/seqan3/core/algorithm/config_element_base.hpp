// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ============================================================================

/*!\file
 * \brief Provides seqan3::detail::config_element_base.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <utility>
#include <type_traits>

#include <seqan3/core/algorithm/concept.hpp>
#include <seqan3/core/algorithm/config_element_access.hpp>
#include <seqan3/std/concept/core_language.hpp>

namespace seqan3::detail
{

/*!\brief Abstract crtp-base class for configurations in use with seqan3::configuration.
 * \ingroup core_algorithm
 *
 * \tparam derived_t The type that should be extended with the config interface.
 *
 * \details
 *
 * This class provides a common interface for config types that are stored in a seqan3::configuration object.
 * It provides getter functions to retrieve the stored state of the config implementation and an additional
 * constructor that takes the seqan3::configuration and copies the state of the corresponding configuration.
 * To allow the base class access to the private members of this config implementation, the class must add a friend declaration for
 * seqan3::detail::config_element_access.
 * The following example demonstrates the usage of this base class for a simple example.
 *
 * ```cpp
 * template <typename t>
 * class my_config = detail::config_element_base<my_config<t>>
 * {
 *     // Grant the base class access to the private member `state`.
 *     friend class detail::config_element_access<my_config<t>>;
 *
 *     t state{};  // Must name the variable `state`.
 * };
 * ```
 *
 * The configuration class must provide a state with the name `state`, which the base class can access via the
 * seqan3::detail::config_element_access struct. This class, then gives access to the underlying data via getter functions.
 * Often, the config is a static type and can be set with an enum value to specify a certain policy for the target
 * algorithm. In case the exact config can also be set at runtime, one can use the seqan3::detail::deferred_config_element_base
 * class to provide functionality of converting the runtime config value to a static config type.
 *
 * \see seqan3::detail::deferred_config_element_base
 */
template <typename derived_t>
class config_element_base
{
protected:

    /*!\name Constructor, destructor and assignment
     * \brief Protected constructor to avoid instantiation of this abstract base class.
     * \{
     */
    config_element_base()                                        = default;
    config_element_base(config_element_base const &)             = default;
    config_element_base(config_element_base &&)                  = default;
    config_element_base & operator=(config_element_base const &) = default;
    config_element_base & operator=(config_element_base &&)      = default;
    ~config_element_base()                                       = default;
    //!\}
public:

    /*!\name Accessor functions
     * \{
     */

    /*!\brief Returns the underlying value associated with this config element.
     *
     * \returns The data of the associated config element.
     *
     * \details
     *
     * Extending classes must provide a data member called `state` and grant access to it by adding
     * seqan3::detail::config_element_access to the type definition.
     *
     * ### Complexity
     *
     * Constant.
     *
     * ### Exception
     *
     * No-throw guarantee.
     *
     * ### Concurrency
     *
     * Thread-safe if data is not written.
     */
    constexpr auto & data() & noexcept
    {
        return config_element_access<derived_t>::_data(static_cast<derived_t &>(*this));
    }

    //!\copydoc config_element_base::data()
    constexpr auto const & data() const & noexcept
    {
        return config_element_access<derived_t>::_data(static_cast<derived_t const &>(*this));
    }

    //!\copydoc config_element_base::data()
    constexpr auto && data() && noexcept
    {
        return config_element_access<derived_t>::_data(std::move(static_cast<derived_t &&>(*this)));
    }

    //!\copydoc config_element_base::data()
    constexpr auto const && data() const && noexcept
    {
        return config_element_access<derived_t>::_data(std::move(static_cast<derived_t const &&>(*this)));
    }
    //!\}
};
} // namespace seqan3::detail
