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
 * \brief Provides seqan3::detail::config_element_access.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <utility>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

// Forward declarations.
//!\cond
template <typename derived_t>
class config_element_base;
//!\endcond

//!\cond
template <typename derived_t>
class deferred_config_element_base;
//!\endcond

/*!\brief An attorney class to grant the crtp-base classes (seqan3::detail::config_element_base and seqan3::detail::deferred_config_element_base)
 *        access to the members of there derived classes.
 * \ingroup core_algorithm
 * \tparam derived_config_t The actual config type which inherits from the respective crtp-base class.
 */
template <typename derived_config_t>
class config_element_access
{
    /*!\brief Grant access to the member variable `state` of the actual config implementation.
     * \param cfg The actual config to get the state from.
     * \returns The state.
     */
    constexpr static auto & _data(derived_config_t & cfg)
    {
        return cfg.state;
    }

    //!\copydoc _data()
    constexpr static auto const & _data(derived_config_t const & cfg)
    {
        return cfg.state;
    }

    //!\copydoc _data()
    constexpr static auto && _data(derived_config_t && cfg)
    {
        return std::move(cfg.state);
    }

    //!\copydoc _data()
    constexpr static auto const && _data(derived_config_t const && cfg)
    {
        return std::move(cfg.state);
    }

    /*!\brief Grant access to the member function `invoke` of the actual config implementation.
     * \param cfg          The actual config to be invoked.
     * \param fn           The callable to be forwarded to the invocation.
     * \param configuration The configuration to be forwarded to the invocation.
     * \returns The result of calling fn at the implementation site.
     */
    template <typename fn_t, typename configuration_t>
    constexpr static auto _invoke(derived_config_t const & cfg,
                                           fn_t && fn,
                                           configuration_t && configuration)
    {
        return cfg.invoke(std::forward<fn_t>(fn), std::forward<configuration_t>(configuration));
    }
    //!\endcond

    //!\brief Grant access to the seqan3::detail::config_element_base class.
    friend class config_element_base<derived_config_t>;
    //!\brief Grant access to the seqan3::detail::deferred_config_element_base class.
    friend class deferred_config_element_base<derived_config_t>;
};

} // namespace seqan3::detail
