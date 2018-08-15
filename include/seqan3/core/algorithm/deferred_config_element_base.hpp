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
 * \brief Provides seqan3::detail::deferred_config_element_base.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/algorithm/concept.hpp>
#include <seqan3/core/metafunction/basic.hpp>

namespace seqan3::detail
{
/*!\brief Abstract crtp-base class for deferred configurations in use with seqan3::configuration.
 * \ingroup algorithm
 *
 * \tparam derived_t The type that should be extended with the config interface.
 *
 * \details
 *
 * This class provides a common interface for deferred config types that are stored in a seqan3::configuration object.
 * It provides a member function `invoke` that must be implemented by `derived_t` as well.
 * Within this `invoke` function the passed configuration is transformed to a new configuration replacing the deferred
 * config type with it's static counter part after resolving the runtime information to a static type or value.
 * The following example demonstrates the usage of this base class with a simple example:
 *
 * ```cpp
 * template <size_t I>
 * struct my_config
 * {
 *     size_t value{I};  // Has to be named `value`.
 * };
 *
 * struct my_deferred_config
 * {
 *      template <typename fn_t, typename configuration_t>
 *      constexpr auto invoke(fn_t && fn, configuration_t && config) const
 *          requires detail::is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
 *      {
 *          if (value == 0)
 *              return fn(std::forward<configuration_t>(cfg).replace_with(my_deferred_config{}, my_config<0>{}));
 *          else
 *              return fn(std::forward<configuration_t>(cfg).replace_with(my_deferred_config{}, my_config<1>{}));
 *      }
 *
 *     int value{0};  // Has to be named `value`.
 * };
 * ```
 *
 * The configuration class must provide a member variable with the name `value`.
 * For a dynamic dispatching of configurations, that should be translated to a static configuration for the
 * target algorithm, the `invoke` function is triggered by the configuration system before the algorithm is executed.
 * In this function, as can be seen above, the runtime parameter can be translated to a static config and then calls the
 * passed callable `fn`, with the new configuration object. Please make sure you pass the old configuration to the new
 * one, such that the other configs are correctly passed to the new configuration.
 * One can use the helper function seqan3::detail::replace_with to easily create a new configuration with the deferred
 * config replaced with the static one.
 *
 * \see seqan3::detail::config_element_base
 */
template <typename derived_t>
class deferred_config_element_base
{
protected:

    //!|brief Give `derived_t` access to protected constructors.
    friend derived_t;

    /*!\name Constructor, destructor and assignment
     * \brief Non public constructor to avoid instantiation of this abstract base class.
     * \{
     */
    deferred_config_element_base()                                                 = default;
    deferred_config_element_base(deferred_config_element_base const &)             = default;
    deferred_config_element_base(deferred_config_element_base &&)                  = default;
    deferred_config_element_base & operator=(deferred_config_element_base const &) = default;
    deferred_config_element_base & operator=(deferred_config_element_base &&)      = default;
    ~deferred_config_element_base()                                                = default;
    //!}

public:

    /*!\brief Invokes the actual translation of the configuration within the `derived_t`
     * \param fn  A callable that is invoked with the altered configuration.
     * \param cfg The old configuration containing the deferred_config (`derived_t`) that is currently invoked.
     *
     * \returns The result of invoking `fn` with the altered configuration.
     *
     * \attention The result is declared `[[nodiscard]]` and cannot be silently ignored.
     */
    template <typename fn_t,
              typename configuration_t>
    [[nodiscard]] constexpr auto operator()(fn_t && fn, configuration_t && cfg) const
    //!\cond
        requires is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
    //!\endcond
    {
        return static_cast<derived_t const &>(*this).invoke(std::forward<fn_t>(fn),
                                                          std::forward<configuration_t>(cfg));
    }
};

} // namespace seqan3::detail
