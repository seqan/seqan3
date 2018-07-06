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
 * \brief Provides seqan3::detail::deferred_config_base.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/algorithm/concept.hpp>
#include <seqan3/core/algorithm/config_base.hpp>
#include <seqan3/core/algorithm/config_access.hpp>

namespace seqan3::detail
{
/*!\brief Abstract crtp-base class for deferred configurations in use with seqan3::configurator.
 * \ingroup core_algorithm
 * \implements seqan3::detail::deferred_config_concept
 *
 * \tparam derived_t The type that should be extended with the config interface.
 *
 * \details
 *
 * This class provides a common interface for deferred config types that are stored in a seqan3::configurator object.
 * It provides getter functions to retrieve the stored state of the config implementation and an additional
 * constructor that takes the seqan3::configurator and copies the state of the corresponding configuration.
 * In addition, this base class requires, that the `derived_t` type provides a member function `invoke`, that
 * transforms the passed configurator to a new configurator replacing the deferred config type with it's static
 * counter part after resolving the runtime information to a static type or value.
 * To allow the base class access to the private members of this config implementation, the class must add a firend declaration for
 * seqan3::detail::config_access.
 * The following example demonstrates the usage of this base class with a simple example:
 *
 * ```cpp
 * template <size_t I>
 * class my_config = detail::config_base<my_config<I>>
 * {
 *     // Grant the base class access to the private member `state`.
 *     friend class detail::config_access<my_config<I>>;
 *
 *     size_t state{I};  // Has to be named `state`.
 * };
 *
 * class my_deferred_config = detail::deferred_config_base<my_deferred_config>
 * {
 *     // Grant the base class access to the private member `state`.
 *     friend class detail::config_access<my_deferred_config>;
 *
 *      template <typename fn_t, detail::configurator_concept configurator_t>
 *      constexpr auto invoke(fn_t && fn, configurator_t && config) const
 *      {
 *          if (state == 0)
 *          {
 *              using new_cfg_t = detail::replace_config_with_t<std::remove_reference_t<configurator_t>,
 *                                                              my_deferred_config,
 *                                                              my_config<0>>;
 *              return fn(new_cfg_t{std::forward<configurator_t>(config)});
 *          }
 *          using new_cfg_t = detail::replace_config_with_t<std::remove_reference_t<configurator_t>,
 *                                                          my_deferred_config,
 *                                                          my_config<1>>;
 *          return fn(new_cfg_t{std::forward<configurator_t>(config)});
 *      }
 *
 *     int state{0};  // Has to be named `state`.
 * };
 * ```
 *
 * The configuration class must provide a member variable with the name `state`, which the base class can access via the
 * seqan3::detail::config_access struct. This class, then gives access to the underlying data via getter functions.
 * For a dynamic dispatching of configurations, that should be translated to a static configuration for the target algorithm,
 * the `invoke` function is triggered by the configuration system before the algorithm is executed. In this function, as
 * can be seen above, the runtime parameter can be translated to a static config and then calls the passed callable `fn`, with
 * the new configurator object. Please make sure, you pass the old configurator to the new one, such that the other configs are
 * correctly passed to the new configurator.
 * One can use the helper meta-function seqan3::detail::replace_config_with to easily create a new configurator with the deferred
 * config replaced with the static one.
 *
 * \see seqan3::detail::config_base
 */
template <typename derived_t>
class deferred_config_base : public config_base<derived_t>
{
protected:

    //!|brief Inherit constructor of base class.
    using config_base<derived_t>::config_base;

    /*!\name Constructor, destructor and assignment
     * \brief Non public constructor to avoid instantiation of this abstract base class.
     * \{
     */
    deferred_config_base()                                         = default;
    deferred_config_base(deferred_config_base const &)             = default;
    deferred_config_base(deferred_config_base &&)                  = default;
    deferred_config_base & operator=(deferred_config_base const &) = default;
    deferred_config_base & operator=(deferred_config_base &&)      = default;
    ~deferred_config_base()                                        = default;
    //!}

public:

    /*!\brief Invokes the actual translation of the configuration within the `derived_t`
     * \param fn  A callable that is invoked with the altered configurator.
     * \param cfg The old configurator containing the deferred_config (`derived_t`) that is currently invoked.
     *
     * \returns The result of invoking `fn` with the altered configurator.
     */
    template <typename fn_t,
              detail::configurator_concept configurator_t>
    constexpr auto operator()(fn_t && fn,
                              configurator_t && cfg) const
    {
        return config_access<derived_t>::_invoke(static_cast<derived_t const &>(*this),
                                                 std::forward<fn_t>(fn),
                                                 std::forward<configurator_t>(cfg));
    }

};

} // namespace seqan3::detail
