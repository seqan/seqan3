// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (chr) 2006-2018, Knut Reinert & Freie Universitaet Berlin
// Copyright (chr) 2016-2018, Knut Reinert & MPI Molekulare Genetik
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
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Free function/metafunction wrappers for alphabets with member functions/types.
 *
 * This shall not need be included manually, just include `alphabet/concept.hpp`.
 */

#pragma once

#include <seqan3/alphabet/detail/alphabet_base.hpp>
#include <seqan3/alphabet/detail/convert.hpp>
#include <seqan3/alphabet/aminoacid/concept.hpp>

namespace seqan3
{

/*!\brief A CRTP-base that refines seqan3::alphabet_base and is used by the amino acids.
 * \ingroup aminoacid
 * \tparam derived_type The CRTP parameter type.
 * \tparam size         The size of the alphabet.
 */
template <typename derived_type, auto size>
class aminoacid_base : public alphabet_base<derived_type, size, char>
{
private:
    //!\brief Type of the base class.
    using base_t = alphabet_base<derived_type, size, char>;

    //!\brief Befriend the base class.
    friend base_t;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr aminoacid_base() : base_t{} {}
    constexpr aminoacid_base(aminoacid_base const &) = default;
    constexpr aminoacid_base(aminoacid_base &&) = default;
    constexpr aminoacid_base & operator=(aminoacid_base const &) = default;
    constexpr aminoacid_base & operator=(aminoacid_base &&) = default;
    ~aminoacid_base() = default;
    //!\}

    //!\brief Befriend the derived class so it can instantiate.
    friend derived_type;

public:

    // Import from base:
    using typename base_t::char_type;
    using typename base_t::rank_type;
    using base_t::value_size;
    using base_t::to_rank;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    // This constructor needs to be public, because constructor templates are not inherited otherwise
    //!\brief Allow explicit construction from any other aminoacid type and convert via the character representation.
    template <typename other_aa_type>
    //!\cond
        requires !std::Same<aminoacid_base, other_aa_type> &&
                 !std::Same<derived_type, other_aa_type> &&
                 aminoacid_concept<other_aa_type>
    //!\endcond
    explicit constexpr aminoacid_base(other_aa_type const & other) noexcept
    {
        using seqan3::to_rank;
        static_cast<derived_type &>(*this) =
            detail::convert_through_char_representation<derived_type, other_aa_type>[to_rank(other)];
    }
    //!\}
};

} // namespace seqan3
