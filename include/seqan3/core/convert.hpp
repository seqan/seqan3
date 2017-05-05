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

/*!\file core/convert.hpp
 * \ingroup core
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::convert().
 */

#pragma once

#include <type_traits>

#include <seqan3/core/concept/core.hpp>

namespace seqan3::detail
{

/*!\brief Default implementation of seqan3::convert(), uses implicit conversion or `static_cast`.
 * \ingroup core
 * \tparam out_t The type to convert to.
 * \tparam in_t The type of the input.
 * \details
 * This class template provides a static function impl() that does the work of seqan3::convert().
 * The implementation is delegated to a class template to be able to use partial template specialisation
 * resolution on the return type.
 *
 * Whenever you wish to specialise / constrain seqan3::convert(), instead specialise / constrain this
 * class template and reimplement impl().
 */
template <typename out_t, typename in_t>
struct convert
{
    static_assert(implicitly_convertible_to_concept<in_t, out_t> || explicitly_convertible_to_concept<in_t, out_t>,
                  "You cannot convert these types.");

    //!\brief The implementation of seqan3::convert().
    //!\sa seqan3::detail::convert
    static constexpr out_t impl(in_t const & in) noexcept
    {
        if constexpr (implicitly_convertible_to_concept<in_t, out_t>)
            return in;
        else
            return static_cast<out_t>(in);
    }
};

} // namespace seqan3::detail

namespace seqan3
{

/*!\brief Convert types to each other explicitly.
 * \ingroup core
 * \tparam out_t The type to convert to (must be given).
 * \tparam in_t The type of the input (deduced from parameter).
 * \param in The input parameter.
 * \returns The input parameter converted to `out_t`.
 * \details
 * The default implementation of this function tries an implicit conversion, then a `static_cast` to
 * convert the parameter. However, other types may provide specialisations with different behaviour.
 * \par Complexity
 * Unless otherwise stated, all specialisations convert in constant time (\f$O(1)\f$).
 * \par Exceptions
 * All specialisations are guaranteed to never throw.
 * \par Thread safety
 * Does not modify data.
 * \par Example
 *
 * ```cpp
 * bool b = convert<bool>(7);  // == true
 * auto i = convert<int>(3.4); // == 3
 *
 * dna4q with_q{dna4::C, 0};
 * auto  wo_q = convert<dna4>(with_q); // == dna4::C
 * ```
 */
template <typename out_t, typename in_t>
constexpr out_t convert(in_t const & in) noexcept
{
    return detail::convert<out_t, in_t>::impl(in);
}

} // namespace seqan3
