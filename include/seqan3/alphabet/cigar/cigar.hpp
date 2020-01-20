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

#include <seqan3/alphabet/adaptation/uint.hpp>
#include <seqan3/alphabet/cigar/cigar_op.hpp>
#include <seqan3/alphabet/composite/alphabet_tuple_base.hpp>
#include <seqan3/core/detail/debug_stream_type.hpp>
#include <seqan3/range/container/small_string.hpp>
#include <seqan3/std/charconv>

// ------------------------------------------------------------------
// cigar
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief The cigar semialphabet pairs a counter with a seqan3::cigar_op letter.
 * \ingroup cigar
 * \implements seqan3::writable_semialphabet
 * \implements seqan3::trivially_copyable
 * \implements seqan3::standard_layout
 * \implements std::regular
 *
 * \details
 *
 * This semialphabet represents a unit in a CIGAR string, typically found in the
 * SAM and BAM formats. It consists of a number and a seqan3::cigar_op symbol.
 *
 * It has a "visual representation", but since this is a string and not a char,
 * the type only models seqan3::writable_semialphabet and not
 * seqan3::writable_alphabet.
 * Members for reading/writing the string are provided.
 *
 * To avoid confusion between string and char literal, this alphabet has
 * no user defined literal operators. Always assign from a pair of
 * uint32_t and seqan3::cigar_op.
 *
 * ### Example
 *
 * \include test/snippet/alphabet/cigar/cigar.cpp
 */
class cigar : public alphabet_tuple_base<cigar, uint32_t, cigar_op>
{
private:
    //!\brief The base class.
    using base_t = alphabet_tuple_base<cigar, uint32_t, cigar_op>;

    //!\cond \brief Befriend seqan3::alphabet_tuple_base.
    friend base_t;
    //!\endcond

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr cigar()                          noexcept = default; //!< Defaulted.
    constexpr cigar(cigar const &)             noexcept = default; //!< Defaulted.
    constexpr cigar(cigar &&)                  noexcept = default; //!< Defaulted.
    constexpr cigar & operator=(cigar const &) noexcept = default; //!< Defaulted.
    constexpr cigar & operator=(cigar &&)      noexcept = default; //!< Defaulted.
    ~cigar()                                   noexcept = default; //!< Defaulted.

    // Inherit constructors from base
    using base_t::base_t;

    /*!\brief Construction via a value of one of the components.
     * \tparam component_type One of the component types; must be uniquely contained in the type list of the composite.
     * \param[in] alph        The value of a component that should be assigned.
     *
     * \include test/snippet/alphabet/cigar/cigar_value_construction.cpp
     */
    SEQAN3_DOXYGEN_ONLY(( constexpr cigar(component_type const alph) noexcept {} ))

    /*!\brief Assignment via a value of one of the components.
     * \tparam component_type One of the component types; must be uniquely contained in the type list of the composite.
     * \param[in] alph        The value of a component that should be assigned.
     *
     * \include test/snippet/alphabet/cigar/cigar_value_assignment.cpp
     */
    SEQAN3_DOXYGEN_ONLY(( constexpr cigar & operator=(component_type const alph) noexcept {} ))
    //!\}

    // Inherit operators from base
    using base_t::operator=;

    /*!\name Read functions
     * \{
     */
    //!\brief Return the string representation.
    small_string<11> to_string() const noexcept
    {
        small_string<11> ret{}; // maximum number of digits for uint32_t + 1 char for the cigar_op
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
    //!\brief Assign from the string representation.
    cigar & assign_string(small_string<11> const s) noexcept
    {
        uint32_t num{};
        auto [ ptr, errc ] = std::from_chars(s.data(), s.data() + 10, num);

        if ((errc != std::errc{}) || (!char_is_valid_for<cigar_op>(*ptr)) || (*(ptr + 1) != 0))
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
     */
    SEQAN3_DOXYGEN_ONLY(( template <size_t index> constexpr auto get(cigar & l) noexcept {} ))

    /*!\copybrief get
     * \tparam type Return the element of specified type; only available if the type is unique in the set of components.
     * \returns A proxy to the contained element that models the same alphabet concepts and supports assignment.
     *
     * \include test/snippet/alphabet/cigar/cigar_get_type.cpp
     */
    SEQAN3_DOXYGEN_ONLY(( template <typename type> constexpr auto get(cigar & l) noexcept {} ))
    //!\}
};

//!\brief Overload for the seqan3::debug_stream's operator.
template <typename char_t>
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & s, cigar const c)
{
    s << c.to_string();
    return s;
}

} // namespace seqan3
