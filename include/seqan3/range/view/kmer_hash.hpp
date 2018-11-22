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
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 * \brief Provides seqan3::view::kmer_hash.
 */

#pragma once

#include <range/v3/view/sliding.hpp>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/metafunction/range.hpp>
#include <seqan3/range/view/detail.hpp>
#include <seqan3/std/view/all.hpp>
#include <seqan3/std/view/transform.hpp>

namespace seqan3::detail
{
// ============================================================================
//  kmer_hash_fn (adaptor definition)
// ============================================================================

//!\brief View adaptor definition for view::kmer_hash.
class kmer_hash_fn : public pipable_adaptor_base<kmer_hash_fn>
{
private:
    //!\brief Type of the CRTP-base.
    using base_t = pipable_adaptor_base<kmer_hash_fn>;

    //!\brief Befriend the base class so it can call impl().
    friend base_t;

    /*!\brief            Call the view's constructor with the underlying view as argument.
     * \param[in] urange The input range to process. Must model std::ranges::ViewableRange and the reference type of the
     *                   range of the range must model seqan3::semi_alphabet_concept.
     * \param[in] k      The k-mer size to construct hashes for.
     * \returns          A range of converted elements.
     */
    template <std::ranges::ViewableRange urng_t>
    //!\cond
        requires semi_alphabet_concept<reference_t<urng_t>>
    //!\endcond
    static auto impl(urng_t && urange, size_t const k)
    {
        return std::forward<urng_t>(urange) | ranges::view::sliding(k) | view::transform(
        [] (auto const in)
        {
            std::hash<decltype(in)> h{};
            return h(in);
        });
    }

public:
    //!\brief Inherit the base class's Constructors.
    using base_t::base_t;
};

} // namespace seqan3::detail

namespace seqan3::view
{
    /*!\name Alphabet related views
     * \{
     */

    /*!\brief               A view that calls std::hash on each substring of length k in the input range.
     * \tparam urng_t       The type of the range being processed. See below for requirements. [template parameter is
     *                      omitted in pipe notation]
     * \param[in] urange    The range being processed. [parameter is omitted in pipe notation]
     * \returns             A range of unsigned integral values where each value is the hash of the resp. k-mer.
     *                      See below for the properties of the returned range.
     * \ingroup view
     *
     * ### View properties
     *
     * | range concepts and reference_t  | `urng_t` (underlying range type)      | `rrng_t` (returned range type)                     |
     * |---------------------------------|:-------------------------------------:|:--------------------------------------------------:|
     * | std::ranges::InputRange         | *required*                            | *preserved*                                        |
     * | std::ranges::ForwardRange       | *required*                            | *preserved*                                        |
     * | std::ranges::BidirectionalRange |                                       | *preserved*                                        |
     * | std::ranges::RandomAccessRange  |                                       | *preserved*                                        |
     * | std::ranges::ContiguousRange    |                                       | *lost*                                             |
     * |                                 |                                       |                                                    |
     * | std::ranges::ViewableRange      | *required*                            | *guaranteed*                                       |
     * | std::ranges::View               |                                       | *guaranteed*                                       |
     * | std::ranges::SizedRange         |                                       | *preserved*                                        |
     * | std::ranges::CommonRange        |                                       | *preserved*                                        |
     * | std::ranges::OutputRange        |                                       | *lost*                                             |
     * | seqan3::const_iterable_concept  |                                       | *preserved*                                        |
     * |                                 |                                       |                                                    |
     * | seqan3::reference_t             | seqan3::semi_alphabet_concept         | std::size_t                                        |
     *
     * See the \link view view submodule documentation \endlink for detailed descriptions of the view properties.
     *
     * ### Example
     * \snippet test/snippet/range/view/kmer_hash.cpp usage
     * \hideinitializer
     */
    inline auto constexpr kmer_hash = detail::kmer_hash_fn{};
} // namespace seqan3::view
