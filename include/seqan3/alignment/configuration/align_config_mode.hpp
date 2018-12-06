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
 * \brief Provides global alignment configurations.
 * \author Joshua Kim <joshua.kim AT fu-berlin.de>
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/configuration/detail.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>

namespace seqan3::detail
{
//!\brief Selects the global alignment mode.
//!\ingroup configuration
struct global_alignment_type
{
    //!\privatesection
    //!\brief An internal id used to check for a valid alignment configuration.
    static constexpr detail::align_config_id id{detail::align_config_id::global};
};

} // namespace seqan3::detail

namespace seqan3::align_cfg
{

//!\brief Selects global alignment mode.
//!\ingroup configuration
inline constexpr detail::global_alignment_type global_alignment;

/*!\brief A configuration element for global alignment.
 * \ingroup configuration
 */
template <typename alignment_t>
//!\cond
    requires std::Same<remove_cvref_t<alignment_t>, detail::global_alignment_type>
//!\endcond
class mode : public pipeable_config_element
{
public:
    //!\privatesection
    //!\brief An internal id used to check for a valid alignment configuration inherited from the alignment policy.
    static constexpr detail::align_config_id id{alignment_t::id};

    /*!\name Constructor, destructor and assignment
     * \brief Defaulted all standard constructor.
     * \{
     */
    constexpr mode()                         = default;
    constexpr mode(mode const &)             = default;
    constexpr mode(mode &&)                  = default;
    constexpr mode & operator=(mode const &) = default;
    constexpr mode & operator=(mode &&)      = default;
    ~mode()                                  = default;

    //!\brief Construction from a specific alignment policy.
    constexpr mode(alignment_t)
    {}
    //!}

    //!\publicsection
    //!\brief The value of align_config_global.
    alignment_t value{};
};

/*!\name Type deduction guides
 * \relates seqan3::align_cfg::mode
 * \{
 */
//!\brief Deduces the alignment mode from the given constructor argument.
template <typename alignment_t>
mode(alignment_t) -> mode<alignment_t>;
//!}
} // namespace seqan3::align_cfg
