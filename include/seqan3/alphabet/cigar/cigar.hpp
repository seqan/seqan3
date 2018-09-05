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
 * \brief Contains the cigar alphabet and auxiliary functions.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alignment/aligned_sequence/aligned_sequence_concept.hpp>
#include <seqan3/alphabet/cigar/cigar_op.hpp>
#include <seqan3/alphabet/composition/cartesian_composition.hpp>
#include <seqan3/std/concepts>
#include <seqan3/std/ranges>

#include <range/v3/view/zip.hpp>

namespace seqan3
{

/*!\brief The cigar alphabet joins the seqan3::cigar_op alphabet with a length value.
 * \ingroup cigar
 * \tparam length_t Type of the length value; must model std::UnsignedIntegral.
 *
 * This alphabet can be used to represent an alignment in a compressed way using
 * the CIGAR string representation. The seqan3::cigar alphabet implements the
 * CIGAR elements in this context consisting of a CIGAR operator (seqan3::cigar_op)
 * and a length value.
 *
 * Example usage:
 * \snippet test/snippet/alphabet/cigar/cigar_op.cpp general
 *
 * \note Usually you do not want to manipulate cigar elements and vectors on
 *       your own but convert an alignment to a cigar and back. See
 *       seqan3::get_cigar_vector for how to convert two aligned sequences into
 *       a seqan3::cigar_vector.
 */
template<typename length_t = uint32_t>
//!\cond
    requires std::UnsignedIntegral<length_t>
//!\endcond
class cigar :
    public cartesian_composition<cigar<length_t>, cigar_op, length_t>
{
private:
    //!\brief The base type.
    using base_type = cartesian_composition<cigar<length_t>, cigar_op, length_t>;

public:

    //!\brief Equals the char_type of cigar_op type.
    using char_type = underlying_char_t<cigar_op>;

    //!\brief The template parameter as a type member.
    using length_type = length_t;

    //!\brief The type combining both types into a streamable string.
    using string_type = std::string;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    cigar() = default;
    constexpr cigar(cigar const &) = default;
    constexpr cigar(cigar &&) = default;
    constexpr cigar & operator=(cigar const &) = default;
    constexpr cigar & operator=(cigar &&) = default;
    ~cigar() = default;

    using base_type::base_type; // Inherit non-default constructors

    using base_type::operator=; // Inherit non-default assignment operators

    //!\copydoc cartesian_composition::cartesian_composition(component_type const alph)
    SEQAN3_DOXYGEN_ONLY(( constexpr cigar(component_type const alph) {} ))
    //!\copydoc cartesian_composition::cartesian_composition(indirect_component_type const alph)
    SEQAN3_DOXYGEN_ONLY(( constexpr cigar(indirect_component_type const alph) {} ))
    //!\copydoc cartesian_composition::operator=(component_type const alph)
    SEQAN3_DOXYGEN_ONLY(( constexpr cigar & operator=(component_type const alph) {} ))
    //!\copydoc cartesian_composition::operator=(indirect_component_type const alph)
    SEQAN3_DOXYGEN_ONLY(( constexpr cigar & operator=(indirect_component_type const alph) {} ))
    //!\}

    /*!\name Write functions
     * \{
     */
    //!\brief Assign from a character. This modifies the internal seqan3::cigar_op letter.
    constexpr cigar & assign_char(char_type const c)
    {
        seqan3::assign_char(get<0>(*this), c);
        return *this;
    }

    //!\brief Assign from a length value. This modifies the internal length value.
    constexpr cigar & assign_length(length_type const c)
    {
        get<1>(*this) = c;
        return *this;
    }
    //!\}

    /*!\name Read functions
     * \{
     */
    //!\brief Return a character. This reads the internal seqan3::cigar_op letter.
    constexpr char_type to_char() const noexcept
    {
        return seqan3::to_char(get<0>(*this));
    }

    //!\brief Return the length value. This reads the internal length value.
    constexpr length_type to_length() const noexcept
    {
        return get<1>(*this);
    }

    //!\brief Return the cigar element as a string. This reads both, the internal seqan3::cigar_op letter and length value.
    string_type to_string() const noexcept
    {
        std::string tmp{std::to_string(get<1>(*this))};
        tmp.push_back(seqan3::to_char(get<0>(*this)));
        return tmp;
    }
    //!\}
};

//!\brief Type deduction guide enables usage of cigar without specifying template arguments.
//!\relates cigar
template <typename length_type>
cigar(cigar_op &&, length_type &&)
    -> cigar<std::decay_t<length_type>>;

//!\brief Type deduction guide enables usage of cigar without specifying template arguments.
//!\relates cigar
template <typename length_type>
cigar(cigar_op const &, length_type &&)
    -> cigar<std::decay_t<length_type>>;

/*!\brief Specialized implementation of seqan3::alphabet_concept::operator<<()
 *        for the seqan3::cigar alphabet. Delegates to seqan3::cigar::to_string().
 * \relates cigar
 * \tparam length_type Type of the length value; must model std::UnsignedIntegral.
 * \param  os          The output stream you are printing to.
 * \param  cig         The seqan3::cigar element that you wish to print.
 * \returns A reference to the output stream.
 */
template<typename length_type>
//!\cond
    requires std::UnsignedIntegral<length_type>
//!\endcond
std::ostream & operator<<(std::ostream & os, cigar<length_type> const & cig)
{
    os << cig.to_string();
    return os;
}

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

//!\brief Alias for a std::vector over the default seqan3::cigar alphabet.
//!\ingroup cigar
using cigar_vector = std::vector<cigar<>>;

/*!\brief Specialized implementation of the stream operator for a container of
 *        seqan3::cigar elements. Delegates to seqan3::cigar::to_string().
 * \relates cigar
 * \tparam length_type Type of the length value; must model std::UnsignedIntegral.
 * \tparam range_type  Must model std::ranges::ForwardRange and its
 *                     reference type must be seqan3::cigar<length_type>.
 * \param  os          The output stream you are printing to.
 * \param  cigv        The container of cigar elements.
 * \returns A reference to the output stream.
 */
template <typename length_type = uint32_t, std::ranges::ForwardRange range_type>
//!\cond
    requires std::UnsignedIntegral<length_type> &&
             std::Same<reference_t<range_type>, cigar<length_type> &>
//!\endcond
std::ostream & operator<<(std::ostream & os, range_type const & cigv)
{
    std::for_each(cigv.begin(), cigv.end(), [&os] (auto cig) { os << cig; } );
    return os;
}

} // namespace seqan3

// ------------------------------------------------------------------
// std-overloads for the tuple-like interface
// ------------------------------------------------------------------

namespace std
{

/*!\brief std::tuple_size overload.
 * [metafunction specialisation for seqan3::cigar]
 */
template<typename length_type>
struct tuple_size<seqan3::cigar<length_type>>
{
    static constexpr size_t value = 2; //!< Composition of seqan3::cigar_op and length value.
};

/*!\brief std::tuple_element overload.
 * [metafunction specialisation for seqan3::cigar]
 */
template<size_t elem_no,
         typename length_type>
struct tuple_element<elem_no, seqan3::cigar<length_type>>
    : tuple_element<elem_no, std::tuple<seqan3::cigar_op, length_type>>
{};

} // namespace std


// ------------------------------------------------------------------
// create cigar_vector out of alignment
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief Compares two aligned sequence values and returns their seqan3::cigar_op.
 * \ingroup aligned_sequence
 * \relatesalso cigar
 * \tparam reference_char_type Must be equality comparable to seqan3::gap.
 * \tparam query_char_type     Must be equality comparable to seqan3::gap.
 * \param  reference_char      The aligned character of the reference to compare.
 * \param  query_char          The aligned character of the query to compare.
 * \param  extended_cigar      Whether to print the extended cigar alphabet or not. See seqan3::cigar_op.
 * \returns A seqan3::cigar_op representing the alignment operation between the
 *          two values.
 *
 * \details
 *
 * \attention Note that seqan3::cigar elements (respectively by their
 *            seqan3::cigar_op) are always related to on one of the two sequences
 *            in a pairwise alignment. In this case, the resulting
 *            operation is related to the \p query_char.
 *
 * ### Example:
 *
 * The following alignment column shows the reference char ('C') on top and a
 * gap for the query char at the bottom.
 * ```
 * ... C ...
 *     |
 * ... - ...
 * ```
 * In this case, the function seqan3::compare_aligned_values will return
 * seqan3::cigar_op::D, since the query char is "deleted".
 *
 * The next alignment column shows the reference char ('C') on top and a
 * query char ('G') at the bottom.
 * ```
 * ... C ...
 *     |
 * ... G ...
 * ```
 * In this case, the function seqan3::compare_aligned_values will return
 * seqan3::cigar_op::M, for the basic cigar the two bases are aligned, while
 * in the extended CIGAR alphabet (\p extended_cigar = `true`) the function
 * will return an seqan3::cigar_op::X since the bases are aligned but are not
 * equal.
 */
template<typename reference_char_type, typename query_char_type>
//!\cond
    requires std::EqualityComparableWith<std::remove_reference_t<reference_char_type>, gap> &&
             std::EqualityComparableWith<std::remove_reference_t<query_char_type>, gap>
//!\endcond
cigar_op compare_aligned_values(reference_char_type const reference_char,
                                query_char_type const query_char,
                                bool const extended_cigar)
{
    return (reference_char == gap::GAP)
                ? (query_char == gap::GAP)
                    ? cigar_op::P
                    : cigar_op::I
                : (query_char == gap::GAP)
                    ? cigar_op::D
                    : (extended_cigar)
                        ? (query_char == reference_char)
                            ? cigar_op::EQ
                            : cigar_op::X
                        : cigar_op::M;
}

/*!\brief Transforms an alignment represented by two aligned sequences into the
 *        corresponding seqan3::cigar_vector.
 * \ingroup aligned_sequence
 * \relatesalso cigar
 * \tparam ref_seq_type    Must model seqan3::aligned_sequence_concept.
 * \tparam query_seq_type  Must model seqan3::aligned_sequence_concept.
 * \param  alignment       The alignment, represented by a pair of aligned sequences,
 *                         to be transformed into seqan3::cigar_vector based on the
 *                         second (query) sequence.
 * \param  query_start_pos The start position of the alignment in the query
 *                         sequence indicating soft-clipping.
 * \param  query_end_pos   The end position of the alignment in the query
 *                         sequence indicating soft-clipping.
 * \param  extended_cigar  Whether to print the extended cigar alphabet or not. See seqan3::cigar_op.
 * \returns A seqan3::cigar_vector representing the alignment.
 *
 * \details
 *
 * \attention Note that seqan3::cigar elements (respectively by their
 *            seqan3::cigar_op) are always related to on one of the two sequences
 *            in a pairwise alignment. In this case, the resulting cigar_vector
 *            is based on sequence at the second position of the \p alignment pair,
 *            namely the query sequence.
 *
 *
 * ### Theoretical Example:
 *
 * The following alignment reference sequence on top and the query sequence at
 * the bottom.
 * ```
 * ATGG--CGTAGAGC
 * |||X  |||X|  |
 * ATGCCCCGTTG--C
 * ```
 * In this case, the function seqan3::compare_aligned_values will return
 * The following cigar string when printed: "4M2I5M2D1M". The extended cigar
 * string would look like this: "3=1X2I3=1X1=2D1=".
 *
 * ### Code Example:
 *
 * \snippet test/snippet/alphabet/cigar/cigar.cpp get_cigar_vector
 */
template<aligned_sequence_concept ref_seq_type, aligned_sequence_concept query_seq_type>
cigar_vector get_cigar_vector(std::pair<ref_seq_type, query_seq_type> const & alignment,
                              uint32_t query_start_pos = 0,
                              uint32_t query_end_pos = 0,
                              bool extended_cigar = false)
{
    if (get<0>(alignment).size() != get<1>(alignment).size())
        throw std::logic_error("The aligned sequences must have the same length.");

    cigar_vector result{};

    if (!get<0>(alignment).size())
        return result; // return empty cigar vector if sequences are empty

    // Add (S)oft-clipping at the start of the read
    if (query_start_pos)
        result.push_back(cigar<>{cigar_op::S, query_start_pos});

    // Create cigar string from alignment
    // -------------------------------------------------------------------------
    // initialize first operation:
    cigar tmp{compare_aligned_values(get<0>(alignment)[0], get<1>(alignment)[0], extended_cigar), 0u};
    // go through alignment columns
    for (auto column : ranges::view::zip(std::get<0>(alignment), std::get<1>(alignment)))
    {
        cigar_op next_op = compare_aligned_values(get<0>(column), get<1>(column), extended_cigar);

        if (tmp == next_op)
        {
            ++get<1>(tmp);
        }
        else
        {
            result.push_back(tmp);
            tmp = cigar<>{next_op, 1u};
        }
    }
    // append last cigar element
    result.push_back(tmp);

    // Add (S)oft-clipping at the end of the read
    if (query_end_pos)
        result.push_back(cigar<>{cigar_op::S, query_end_pos});

    return result;
}

} // namespace seqan3
