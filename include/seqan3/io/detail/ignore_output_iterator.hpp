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
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \brief Provides seqan3::detail::ignore_output_iterator for writing to `null` stream.
 */

#pragma once

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

/*!\brief An output iterator that emulates writing to a `null`-stream in order to dispose the output.
 * \ingroup io
 *
 * \details
 *
 * This helper iterator can be used to explicitly dispose output streams, i.e. the output data is transmitted nowhere.
 * A typical use case is when extracted bytes from an input stream should be ignored, as they contain only structural
 * information, e.g. a newline character.
 */
class ignore_output_iterator
{
public:

    /*!\name Member types
     * \{
     * \brief Associated types are void for output iterators, see also
     * [output iterator concept](http://en.cppreference.com/w/cpp/concept/OutputIterator).
     */
    using value_type        = void;
    using reference         = void;
    using pointer           = void;
    using difference_type   = std::ptrdiff_t;
    using iterator_category = std::output_iterator_tag;
    //!\}

    /*!\name Constructor, destructor and assignment
     * \{
     * \brief All non-user defined constructors are explicitly defaulted.
     */
    ignore_output_iterator() = default;
    ignore_output_iterator(ignore_output_iterator const &) = default;
    ignore_output_iterator(ignore_output_iterator &&) = default;
    ignore_output_iterator & operator= (ignore_output_iterator const &) = default;
    ignore_output_iterator & operator= (ignore_output_iterator &&) = default;
    ~ignore_output_iterator() = default;
    //!\}

    /*!\name Member functions
     * \brief Each function performs no operation. In fact writing to the seqan3::detail::ignore_output_iterator,
     *        is subject to removal by compiler optimizations.
     * \{
     */
    //!\brief Emulates writing the passed value to the `null`-stream.
    template <typename type>
    constexpr ignore_output_iterator & operator= (type const /*v*/) noexcept
    {
        return *this;
    }

    //!\brief This operator performs no function in output iterators.
    constexpr ignore_output_iterator & operator* () noexcept
    {
        return *this;
    }

    //!\brief This operator performs no function in output iterators.
    constexpr ignore_output_iterator & operator++ () noexcept
    {
        return *this;
    }

    //!\brief This operator performs no function in output iterators.
    constexpr ignore_output_iterator & operator++ (int) noexcept
    {
        return *this;
    }
    //!\}
};

} // namespace seqan3::detail
