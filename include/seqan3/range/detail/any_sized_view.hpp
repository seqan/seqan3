// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
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
 * \brief Provides the seqan3::detail::any_sized_view class.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \ingroup view
 */

#pragma once

#include <range/v3/view/any_view.hpp>

#include <seqan3/range/view/concept.hpp>

namespace seqan3::detail
{

/*!\brief The equivalent of ranges::any_view<ref_t, c>
 * \tparam ref_t The reference type of the view.
 * \tparam c The ranges::category type that will be supported.
 *
 * \details
 *
 * Allows view type erasure to the category defined by c, but additionally preserves seqan3::sized_range_concept.
 */
template <typename ref_t, ranges::category c>
class any_sized_view : public ranges::any_view<ref_t, c>
{
public:
    //!\brief Publicly accessible size_type.
    using size_type = std::size_t;
    //!\brief Publicly accessible difference_type.
    using difference_type = std::make_signed_t<size_type>;

    //!\brief Construct from another range and save its size.
    //!\tparam irng_t Type of the other range; must satisfy seqan3::input_range_concept and seqan3::sized_range_concept
    template <typename irng_t>
    //!\cond
        requires input_range_concept<irng_t> && sized_range_concept<irng_t>
    //!\endcond
    any_sized_view(irng_t && irange) :
        ranges::any_view<ref_t, c>(std::forward<irng_t>(irange)),
        s{ranges::size(irange)}
    {}

    //!\brief The size of the range. This will also be picked up by the free functions std::size or ranges::size.
    constexpr size_type size() const
    {
        return s;
    }
private:
    //!\brief The member that stores the size.
    size_type s{0};
};

} // namespace seqan3::detail

