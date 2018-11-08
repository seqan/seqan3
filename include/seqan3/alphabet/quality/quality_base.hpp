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
 * \author Marie Hoffmann <marie.hoffmann AT fu-berlin.de>
 * \brief Contains seqan3::phred42 quality scores.
 */

#pragma once

#include <seqan3/alphabet/detail/alphabet_base.hpp>
#include <seqan3/alphabet/detail/convert.hpp>
#include <seqan3/alphabet/quality/concept.hpp>

namespace seqan3
{

/*!\brief A CRTP-base that refines seqan3::alphabet_base and is used by the quality alphabets.
 * \ingroup quality
 * \tparam derived_type The CRTP parameter type.
 * \tparam size         The size of the alphabet.
 */
template <typename derived_type, auto size>
class quality_base : public alphabet_base<derived_type, size, char>
{
public:
    /*!\name Member types
     * \{
     */
    //!\brief The integer representation of a quality score assignable with =operator.
    using phred_type = int8_t;
    //!\}

private:
    //!\brief The base type.
    using base_t = alphabet_base<derived_type, size, char>;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr quality_base() : base_t{} {}
    constexpr quality_base(quality_base const &) = default;
    constexpr quality_base(quality_base &&) = default;
    constexpr quality_base & operator=(quality_base const &) = default;
    constexpr quality_base & operator=(quality_base &&) = default;
    ~quality_base() = default;

    //!\brief Allow construction from the phred value.
    constexpr quality_base(phred_type const p) noexcept
    {
        static_cast<derived_type *>(this)->assign_phred(p);
    }
    //!\}

    //!\brief Befriend the derived_type so it can instantiate.
    friend derived_type;

public:
    // Import from base type:
    using base_t::value_size;
    using base_t::to_rank;
    using base_t::assign_rank;
    using typename base_t::char_type;
    using typename base_t::rank_type;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    // This constructor needs to be public, because constructor templates are not inherited otherwise
    //!\brief Allow explicit construction from any other quality type by means of the phred representation.
    template <typename other_qual_type>
    //!\cond
        requires !std::Same<quality_base, other_qual_type> &&
                 !std::Same<derived_type, other_qual_type> &&
                 quality_concept<other_qual_type>
    //!\endcond
    explicit constexpr quality_base(other_qual_type const & other) noexcept
    {
        using seqan3::assign_phred;
        using seqan3::to_phred;
        assign_phred(static_cast<derived_type &>(*this), to_phred(other));
    }
    //!\}

    /*!\name Read functions
     * \{
     */
    //!\copydoc seqan3::alphabet_base::to_char
    constexpr char_type to_char() const noexcept
    {
        return static_cast<char_type>(to_rank()) + derived_type::offset_char;
    }

    //!\brief Return the alphabet's value in phred representation.
    constexpr phred_type to_phred() const noexcept
    {
        return static_cast<phred_type>(to_rank()) + derived_type::offset_phred;
    }
    //!\}

    /*!\name Write functions
     * \{
     */

    /*!\brief Assign from the numeric phred value.
     *
     * \details
     *
     * Satisfies the seqan3::quality_concept::assign_phred() requirement via the seqan3::assign_rank() wrapper.
     *
     * \par Complexity
     *
     * Constant.
     */
    constexpr derived_type & assign_phred(phred_type const p) noexcept
    {
        return assign_rank(phred_to_rank[static_cast<rank_type>(p)]);
    }
    //!\}

protected:
    //!\brief Phred to rank conversion table.
    static std::array<rank_type, 256> constexpr phred_to_rank
    {
        [] () constexpr
        {
            std::array<rank_type, 256> ret{};

            for (int64_t i = std::numeric_limits<phred_type>::lowest(); i <= std::numeric_limits<phred_type>::max(); ++i)
            {
                if (i < derived_type::offset_phred)                     // map too-small to smallest possible
                    ret[static_cast<rank_type>(i)] = 0;
                else if (i >= derived_type::offset_phred + value_size)  // map too-large to highest possible
                    ret[static_cast<rank_type>(i)] = value_size - 1;
                else                                                    // map valid range to identity
                    ret[static_cast<rank_type>(i)] = i - derived_type::offset_phred;
            }
            return ret;
        }()
    };

    //!\brief Char to rank conversion table.
    static std::array<rank_type, 256> constexpr char_to_rank
    {
        [] () constexpr
        {
            std::array<rank_type, 256> ret{};

            for (int64_t i = std::numeric_limits<char_type>::lowest(); i <= std::numeric_limits<char_type>::max(); ++i)
            {
                if (i < derived_type::offset_char)                     // map too-small to smallest possible
                    ret[static_cast<rank_type>(i)] = 0;
                else if (i >= derived_type::offset_char + value_size)  // map too-large to highest possible
                    ret[static_cast<rank_type>(i)] = value_size - 1;
                else                                                   // map valid range to identity
                    ret[static_cast<rank_type>(i)] = i - derived_type::offset_char;
            }

            return ret;
        }()
    };
};

} // namespace seqan3
