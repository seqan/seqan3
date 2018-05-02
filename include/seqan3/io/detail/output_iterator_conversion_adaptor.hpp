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
 * \brief Output iterator adaptor that takes care of char to alphabet conversion in the assignment.
 */

#pragma once

#include <utility>

#include <range/v3/range_traits.hpp>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/io/concept.hpp>
#include <seqan3/range/concept.hpp>
#include <seqan3/range/container/concept.hpp>
#include <seqan3/std/concept/core_language.hpp>
#include <seqan3/std/concept/iterator.hpp>

namespace seqan3::detail
{

/*!\brief Output iterator adaptor that converts `char` into the given target type when assigning to it.
 * \tparam oiter_type The adapted output_iterator.
 * \tparam alpha_type The target alphabet type to convert to.
 * \ingroup io
 *
 * \details
 *
 * Adapts any output iterator to a char convertible iterator, such that explicit conversions
 * from char to the `alpha_type` is taken care of during assignment.
 *
 * ### Example
 *
 * ```cpp
 * using insert_iter = ranges::back_insert_iterator<std::vector<dna4>>;
 * using out_iter    = seqan3::detail::output_iterator_conversion_adaptor<insert_iter, dna4>;
 *
 * std::vector<dna4> out_vec;
 * out_iter it{insert_iter{out_vec}};
 *
 * *it = 'A';
 * *it = 'C';
 * *it = 'G';
 * *it = 'T';
 *
 * for (auto && c : out_vec)
 *     std::cout << to_char(c);
 * std::cout << '\n'; // prints ACGT
 * ```
 */
template <typename oiter_type, alphabet_concept alpha_type = char>
//!\cond
    requires output_iterator_concept<oiter_type, alpha_type>
//!\endcond
class output_iterator_conversion_adaptor
{
private:
    //!\brief The adapted output iterator.
    oiter_type oiter;

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
    output_iterator_conversion_adaptor() = default;
    output_iterator_conversion_adaptor(output_iterator_conversion_adaptor const &) = default;
    output_iterator_conversion_adaptor(output_iterator_conversion_adaptor &&) = default;
    output_iterator_conversion_adaptor & operator= (output_iterator_conversion_adaptor const &) = default;
    output_iterator_conversion_adaptor & operator= (output_iterator_conversion_adaptor &&) = default;
    ~output_iterator_conversion_adaptor() = default;

    //!\brief Construction from the adopted iterator.
    output_iterator_conversion_adaptor(oiter_type _oiter) : oiter{_oiter}
    {}
    //!\}

    /*!\name Member functions
     * \{
     */
    //!\brief Inserts an object into the associated writable object.
    constexpr output_iterator_conversion_adaptor & operator= (alpha_type const c)
    {
        *oiter = c;
        ++oiter;
        return *this;
    }

    //!\brief Inserts an object into the associated writable object from char.
    constexpr output_iterator_conversion_adaptor & operator= (char const c)
    //!\cond
        requires !same_concept<alpha_type, char>
    //!\endcond
    {  // Explicit conversion through char conversion.
        *oiter = assign_char(alpha_type{}, c);
        ++oiter;
        return *this;
    }

    //!\brief This operator performs no function in output iterators.
    constexpr output_iterator_conversion_adaptor & operator* ()
    {
        return *this;
    }

    //!\brief This operator performs no function in output iterators.
    constexpr output_iterator_conversion_adaptor & operator++ ()
    {
        return *this;
    }

    //!\brief This operator performs no function in output iterators.
    constexpr output_iterator_conversion_adaptor & operator++ (int)
    {
        return *this;
    }
    //!\}
};

/*!\name Convenience functions
 * \{
 */
/*!\brief     A convenience function template that constructs a seqan3::detail::output_iterator_conversion_adaptor for
 *            the container `c` with the type deduced from the type of the argument.
 * \tparam    container_type Any type that satisfies the seqan3::sequence_concept.
 * \param[in] c The instance to construct a push back iterator for.
 * \returns   seqan3::detail::output_iterator_conversion_adaptor over the output container.
 * \ingroup   io
 * \relates   output_iterator_conversion_adaptor
 */
template <sequence_concept container_type>
//!\cond
    requires !std::is_const_v<container_type>
//!\endcond
auto make_conversion_output_iterator(container_type & c)
{
    using value_type  = typename container_type::value_type;
    using insert_iter = ranges::back_insert_iterator<container_type>;
    using out_iter    = output_iterator_conversion_adaptor<insert_iter, value_type>;

    return out_iter{insert_iter{c}};
}

/*!\brief     A convenience function template that constructs an seqan3::detail::output_iterator_conversion_adaptor for
 *            the ostream `s` with the char type and traits type deduced from the type of the argument.
 * \tparam    ostream_type Any type that satisfies the seqan3::ostream_concept.
 * \param[in] s The stream instance to construct the ostreambuf_iterator for.
 * \returns   seqan3::detail::output_iterator_conversion_adaptor over the output stream.
 * \ingroup   io
 * \relates   output_iterator_conversion_adaptor
 */
template <typename ostream_type>
//!\cond
    requires ostream_concept<std::remove_reference_t<ostream_type>,
                             typename ranges::value_type<std::remove_reference_t<ostream_type>>::type>
//!\endcond
auto make_conversion_output_iterator(ostream_type & s)
{
    using char_type    = typename std::remove_reference_t<ostream_type>::char_type;
    using traits_type  = typename std::remove_reference_t<ostream_type>::traits_type;
    using ostream_iter = ranges::ostreambuf_iterator<char_type, traits_type>;
    using out_iter     = output_iterator_conversion_adaptor<ostream_iter, char_type>;

    return out_iter{ostream_iter{s}};
}
//!}

} // namespace seqan3::detail
