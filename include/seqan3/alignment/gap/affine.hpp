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
 * \brief Provides seqan3::align::gap_affine.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/gap/detail.hpp>
#include <seqan3/core/detail/strong_type.hpp>

namespace seqan3
{

/*!\brief Type for gap open costs.
 *\ingroup alignment
 */
template <typename value_t>
struct gap_open_cost : detail::strong_type<value_t, gap_open_cost<value_t>>
{
    //!\brief Inheriting constructors from base class.
    using detail::strong_type<value_t, gap_open_cost<value_t>>::strong_type;
};

/*!\brief Type for gap extension costs.
 *\ingroup alignment
 */
template <typename value_t>
struct gap_extend_cost : detail::strong_type<value_t, gap_extend_cost<value_t>>
{
    //!\brief Inheriting constructors from base class.
    using detail::strong_type<value_t, gap_extend_cost<value_t>>::strong_type;
};

/*!\name Deduction guides
 * \brief Deduces template parameter from the argument.
 * \relates seqan3::gap_open_cost
 * \{
 */

template <typename value_t>
gap_open_cost(value_t) -> gap_open_cost<value_t>;
//!\}


/*!\name Deduction guides
 * \brief Deduces template parameter from the argument.
 * \relates seqan3::gap_extend_cost
 * \{
 */

template <typename value_t>
gap_extend_cost(value_t) -> gap_extend_cost<value_t>;
//!\}


/*!\brief Data structure for affine gaps
 * \ingroup alignment
 * \tparam value_t The value type for the costs.
 */
template <typename value_t>
struct gap_affine
{

    /*!\name Constructor, destructor and assignment
     * \{
     */
    gap_affine()                               = default;
    gap_affine(gap_affine const &)             = default;
    gap_affine(gap_affine &&)                  = default;
    gap_affine & operator=(gap_affine const &) = default;
    gap_affine & operator=(gap_affine &&)      = default;
    ~gap_affine()                              = default;

    //!\brief Construction from seqan3::gap_open_cost and seqan3::gap_extend_cost.
    template <typename inner_value_t>
    constexpr gap_affine(gap_open_cost<inner_value_t> const open_cost,
                         gap_extend_cost<inner_value_t> const extend_cost) noexcept
        : gap_open_cost{std::move(open_cost.get())}, gap_extend_cost{std::move(extend_cost.get())}
    {}
    //!}

    //!\privatesection
    //!\brief The data member storing the gap open cost.
    value_t gap_open_cost{-3};
    //!\brief The data member storing the gap extend cost.
    value_t gap_extend_cost{-1};
};

/*!\name Deduction guides
 * \brief Deduces template parameter from the argument.
 * \relates seqan3::gap_affine
 * \{
 */
template <typename value_t>
gap_affine(gap_open_cost<value_t>, gap_extend_cost<value_t>) -> gap_affine<value_t>;
//!\}

} // namespace seqan3

namespace seqan3::detail
{

//!\cond
template <typename value_t>
struct is_gap_config<gap_affine<value_t>> : public std::true_type
{};
//!\endcond
} // namespace seqan3::detail
