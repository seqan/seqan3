// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the seqan3::cigar alphabet.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/charconv>

#include <seqan3/alphabet/adaptation/uint.hpp>
#include <seqan3/alphabet/cigar/exposition_only/cigar_operation.hpp>
#include <seqan3/alphabet/composite/alphabet_tuple_base.hpp>
#include <seqan3/core/detail/debug_stream_type.hpp>
#include <seqan3/utility/container/small_string.hpp>

// ------------------------------------------------------------------
// cigar
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief The seqan3::cigar semialphabet pairs a counter with a seqan3::cigar::operation letter.
 * \ingroup cigar
 * \implements seqan3::writable_semialphabet
 * \implements seqan3::trivially_copyable
 * \implements seqan3::standard_layout
 * \implements std::regular
 *
 * \details
 *
 * This semialphabet represents a unit in a CIGAR string, typically found in the
 * SAM and BAM formats. It consists of a number and a seqan3::cigar::operation symbol.
 *
 * It has a "visual representation", but since this is a string and not a char,
 * the type only models seqan3::writable_semialphabet and not
 * seqan3::writable_alphabet.
 * Members for reading/writing the string are provided.
 *
 * To avoid confusion between string and char literal, this alphabet has
 * no user defined literal operators. Always assign from a pair of
 * uint32_t and seqan3::cigar::operation.
 *
 * \copydoc seqan3::doxygen::cigar_operation_table
 *
 * ### Example
 *
 * \include test/snippet/alphabet/cigar/cigar.cpp
 *
 * \sa https://samtools.github.io/hts-specs/SAMv1.pdf#page=8
 */
class cigar : public alphabet_tuple_base<cigar, uint32_t, exposition_only::cigar_operation>
{
private:
    //!\brief The base class.
    using base_t = alphabet_tuple_base<cigar, uint32_t, exposition_only::cigar_operation>;

    //!\cond \brief Befriend seqan3::alphabet_tuple_base.
    friend base_t;
    //!\endcond

public:

    /*!\brief The (extended) cigar operation alphabet of M,D,I,H,N,P,S,X,=.
     *
     * \details
     *
     * The CIGAR string can be either basic or extended. The only difference in the extended cigar alphabet is that
     * aligned bases are classified as an actual match ('=') or mismatch ('X'). In contrast, the basic cigar alphabet
     * only indicated the aligned status with an 'M', without further information if the bases are actually equal or
     * not.
     *
     * The main purpose of the seqan3::cigar::operation alphabet is to be used in the seqan3::cigar
     * composition, where a cigar operation is paired with a count value.
     *
     * \copydoc seqan3::doxygen::cigar_operation_table
     *
     * Example usage:
     * \include test/snippet/alphabet/cigar/cigar_operation.cpp
     *
     * \if DEV
     * \note Usually you do not want to manipulate cigar elements and vectors on
     *       your own but convert an alignment to a cigar and back. See
     *       seqan3::detail::get_cigar_vector for how to convert two aligned sequences into
     *       a cigar_vector.
     * \endif
     *
     * \sa https://samtools.github.io/hts-specs/SAMv1.pdf#page=8
     *
     * \stableapi{Since version 3.1.}
     */
    using operation = exposition_only::cigar_operation;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr cigar() noexcept = default; //!< Defaulted.
    constexpr cigar(cigar const &) noexcept = default; //!< Defaulted.
    constexpr cigar(cigar &&) noexcept = default; //!< Defaulted.
    constexpr cigar & operator=(cigar const &) noexcept = default; //!< Defaulted.
    constexpr cigar & operator=(cigar &&) noexcept = default; //!< Defaulted.
    ~cigar() noexcept = default; //!< Defaulted.

    // Inherit constructors from base
    using base_t::base_t;

    /*!\brief Construction via a value of one of the components.
     * \tparam component_type One of the component types; must be uniquely contained in the type list of the composite.
     * \param[in] alph        The value of a component that should be assigned.
     *
     * \include test/snippet/alphabet/cigar/cigar_value_construction.cpp
     *
     * \stableapi{Since version 3.1.}
     */
    SEQAN3_DOXYGEN_ONLY(( constexpr cigar(component_type const alph) noexcept {} ))

    /*!\brief Assignment via a value of one of the components.
     * \tparam component_type One of the component types; must be uniquely contained in the type list of the composite.
     * \param[in] alph        The value of a component that should be assigned.
     *
     * \include test/snippet/alphabet/cigar/cigar_value_assignment.cpp
     *
     * \stableapi{Since version 3.1.}
     */
    SEQAN3_DOXYGEN_ONLY(( constexpr cigar & operator=(component_type const alph) noexcept {} ))

    // Inherit operators from base
    using base_t::operator=;
    //!\}

    /*!\name Read functions
     * \{
     */
    /*!\brief Return the string representation.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    small_string<11> to_string() const noexcept
    {
        small_string<11> ret{}; // maximum number of digits for uint32_t + 1 char for the cigar operation
        ret.resize(11);

        auto [ ptr, errc ] = std::to_chars(ret.data(), ret.data() + 10, get<0>(*this));

        *ptr = seqan3::to_char(get<1>(*this));
        (void)errc;

        ret.resize(ptr - ret.data() + 1);
        return ret;
    }
    //!\}

    /*!\name Write functions
     * \{
     */
    /*!\brief Assign from the string representation.
     * \details
     * \experimentalapi{Experimental since version 3.1.}
     */
    cigar & assign_string(small_string<11> const s) noexcept
    {
        uint32_t num{};
        auto [ ptr, errc ] = std::from_chars(s.data(), s.data() + 10, num);

        if ((errc != std::errc{}) || (!char_is_valid_for<operation>(*ptr)) || (*(ptr + 1) != 0))
        {
            get<0>(*this) = 0;
            assign_char_to('P', get<1>(*this));
        }
        else
        {
            get<0>(*this) = num;
            assign_char_to(*ptr, get<1>(*this));
        }

        return *this;
    }
    //!\}

    /*!\name Get functions
     * \{
     */
    /*!\copydoc alphabet_tuple_base::get(alphabet_tuple_base & l)
     *
     * \include test/snippet/alphabet/cigar/cigar_get_index.cpp
     *
     * \stableapi{Since version 3.1.}
     */
    SEQAN3_DOXYGEN_ONLY(( friend template <size_t index> constexpr auto get(cigar & l) noexcept {} ))

    /*!\copybrief get
     * \tparam type Return the element of specified type; only available if the type is unique in the set of components.
     * \returns A proxy to the contained element that models the same alphabet concepts and supports assignment.
     *
     * \include test/snippet/alphabet/cigar/cigar_get_type.cpp
     *
     * \stableapi{Since version 3.1.}
     */
    SEQAN3_DOXYGEN_ONLY(( friend template <typename type> constexpr auto get(cigar & l) noexcept {} ))
    //!\}
};

//!\brief Overload for the seqan3::debug_stream's operator.
template <typename char_t>
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & s, cigar const c)
{
    s << c.to_string();
    return s;
}

inline namespace literals
{

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

/*!\name Other literals
 * \{
 */

/*!\brief The seqan3::cigar::operation char literal.
 * \relatesalso seqan3::cigar
 * \returns seqan3::cigar::operation
 *
 * You can use this char literal to assign a seqan3::cigar_operation character:
 * \include test/snippet/alphabet/cigar/cigar_operation_char_literal.cpp
 *
 * \stableapi{Since version 3.1.}
 */
inline cigar::operation operator""_cigar_operation(char const c) noexcept
{
    return cigar::operation{}.assign_char(c);
}

//!\deprecated Please use seqan3::""_cigar_operation instead.
SEQAN3_DEPRECATED_310 inline cigar::operation operator""_cigar_op(char const c) noexcept
{
    return cigar::operation{}.assign_char(c);
}
//!\}

} // inline namespace literals

} // namespace seqan3
