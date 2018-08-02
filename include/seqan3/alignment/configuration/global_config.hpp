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
 * \brief Provides seqan3::align::global_config
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/algorithm/all.hpp>
#include <seqan3/core/metafunction/basic.hpp>

namespace seqan3::detail
{

// ----------------------------------------------------------------------------
// global_alignment_cfg
// ----------------------------------------------------------------------------

//!\brief The global_alignment_configuration element to distinguish between alignment policies.
class global_alignment_cfg : public config_element_base<global_alignment_cfg>
{
    //!\brief Grant access to private state.
    friend class config_element_access<global_alignment_cfg>;

    //!\brief State is ignored.
    remove_cvref_t<decltype(std::ignore)> state;
};
}  // namespace seqan3::detail

namespace seqan3::align
{

/*!\brief Global alignment configuration object.
 * \ingroup align
 */
template <seqan3::detail::config_element_concept ... cfg_elements_t>
class global_config : public detail::configuration<cfg_elements_t...>
{
    //!\brief The base type.
    using base_type = detail::configuration<cfg_elements_t...>;

public:
    /*!\name Constructor, destructor and assignment
     * \{
     */
    global_config()                                  = delete;
    global_config(global_config const &)             = default;
    global_config(global_config &&)                  = default;
    global_config & operator=(global_config const &) = default;
    global_config & operator=(global_config &&)      = default;
    ~global_config()                                 = default;

    //TODO(rrahn): Using the intermediate concept form here causes an ICE!
    //!\brief Construction from another seqan3::detail::configuration object.
    template <typename ... other_configs_t>
    //!\cond
        requires (detail::config_element_concept<other_configs_t> && ...)
    //!\endcond
    global_config(detail::configuration<other_configs_t...> const & other) : base_type{other}
    {}

    //TODO(rrahn): Using the intermediate concept form here causes an ICE!
    //!\brief Construction from another seqan3::detail::configuration object.
    template <typename ... other_configs_t>
    //!\cond
        requires (detail::config_element_concept<other_configs_t> && ...)
    //!\endcond
    global_config(detail::configuration<other_configs_t...> && other) : base_type{std::move(other)}
    {}
    //!}
};

/*!\name Deduction guides
 * \relates seqan3::align::global_config
 * \{
 */

//!\brief Deduces template arguments from constructing with another seqan3::detail::configuration.
template <seqan3::detail::config_element_concept ... cfg_elements_t>
global_config(detail::configuration<cfg_elements_t...> const &) ->
    global_config<detail::global_alignment_cfg, cfg_elements_t...>;

//!\brief Deduces template arguments from constructing with another seqan3::detail::configuration.
template <seqan3::detail::config_element_concept ... cfg_elements_t>
global_config(detail::configuration<cfg_elements_t...> &&) ->
    global_config<detail::global_alignment_cfg, cfg_elements_t...>;
//!\}
} // namespace seqan3::align
