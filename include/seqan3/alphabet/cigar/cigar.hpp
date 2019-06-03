// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the seqan3::cigar alphabet.
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <seqan3/alphabet/cigar/cigar_op.hpp>
#include <seqan3/alphabet/composite/alphabet_tuple_base.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/std/charconv>

// ------------------------------------------------------------------
// cigar
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief The cigar alphabet pairs a counter with a seqan3::cigar_op letter.
 * \ingroup cigar
 * \implements seqan3::WritableAlphabet
 * \implements seqan3::TriviallyCopyable
 * \implements seqan3::StandardLayout
 * \implements std::Regular
 *
 * \details
 *
 * This alphabet represents a unit in a CIGAR string, typically found in the
 * SAM and BAM formats. It consists of a number and a seqan3::cigar_op symbol.
 *
 * The "char representation" of this alphabet is a seqan3::small_string. This
 * has some implications when using it in code that doesn't expect a range
 * of symbols.
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

    //!\copydoc alphabet_tuple_base::alphabet_tuple_base(component_type const alph)
    SEQAN3_DOXYGEN_ONLY(( constexpr cigar(component_type const alph) noexcept {} ))
    //!\copydoc alphabet_tuple_base::alphabet_tuple_base(indirect_component_type const alph)
    SEQAN3_DOXYGEN_ONLY(( constexpr cigar(indirect_component_type const alph) noexcept {} ))
    //!\copydoc alphabet_tuple_base::operator=(component_type const alph)
    SEQAN3_DOXYGEN_ONLY(( constexpr cigar & operator=(component_type const alph) noexcept {} ))
    //!\copydoc alphabet_tuple_base::operator=(indirect_component_type const alph)
    SEQAN3_DOXYGEN_ONLY(( constexpr cigar & operator=(indirect_component_type const alph) noexcept {} ))
    //!\}

    // Inherit operators from base
    using base_t::operator=;
    using base_t::operator==;
    using base_t::operator!=;
    using base_t::operator>=;
    using base_t::operator<=;
    using base_t::operator<;
    using base_t::operator>;

    /*!\name Read functions
     * \{
     */
    //!\brief Return the character representation.
    small_string<11> to_char() const noexcept
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
    //!\brief Assign from the character representation.
    cigar & assign_char(small_string<11> const s) noexcept
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
};

//!\brief Overload for the seqan3::debug_stream's operator.
inline debug_stream_type & operator<<(debug_stream_type & s, cigar const c)
{
    s << to_char(c);
    return s;
}

} // namespace seqan3
