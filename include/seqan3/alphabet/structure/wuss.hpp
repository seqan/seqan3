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
 * \ingroup structure
 * \author Joerg Winkler <j.winkler AT fu-berlin.de>
 * \brief Contains the WUSS format for RNA structure.
 */

#pragma once

#include <cassert>
#include <string>
#include <vector>

#include <seqan3/alphabet/structure/concept.hpp>

// ------------------------------------------------------------------
// wuss
// ------------------------------------------------------------------

namespace seqan3
{

/*!\brief The WUSS structure alphabet of the characters ".<>_,:~()[]{}AaBbCcDd...Ss".
 * \ingroup structure
 *
 * \details
 * The symbols ._,:~ denote unpaired characters, brackets <>()[]{} represent base pair interactions and AaBb...Ss
 * form pseudoknots in the structure. The default alphabet has size 51 (letters until Ss). The size can be varied with
 * the optional template parameter between 13 (no letters for pseudoknots) and 65 (all Aa-Zz for pseudoknots).
 *
 *~~~~~~~~~~~~~~~{.cpp}
 * <<<___>>>,,<<<__>>>
 * <<<<_AAAA____>>>>aaaa
 *~~~~~~~~~~~~~~~
 */

template <uint8_t SIZE = 51>
struct wuss
{

    //!\brief The type of the alphabet when converted to char (e.g. via to_char()).
    using char_type = char;
    //!\brief The type of the alphabet when represented as a number (e.g. via to_rank()).
    using rank_type = uint8_t;

    /*!\name Letter values
     * \brief Static member "letters" that can be assigned to the alphabet or used in aggregate initialization.
     * \details Similar to an Enum interface. *Don't worry about the `internal_type`.*
     */
    //!\{
    static const wuss NP; // not paired .
    static const wuss BL; // bracket left <
    static const wuss BR; // bracket right >
    static const wuss NP1; // not paired :
    static const wuss NP2; // not paired ,
    static const wuss NP3; // not paired _
    static const wuss NP4; // not paired ~
    static const wuss BL1; // bracket left (
    static const wuss BR1; // bracket right )
    static const wuss BL2; // bracket left [
    static const wuss BR2; // bracket right ]
    static const wuss BL3; // bracket left {
    static const wuss BR3; // bracket right }
    // pseudoknots not accessible
    //!\}

    /*!\name Read functions
     * \{
     */
    //!\brief Return the letter as a character of char_type.
    constexpr char_type to_char() const noexcept
    {
        return value_to_char[static_cast<rank_type>(_value)];
    }

    //!\brief Return the letter's numeric value or rank in the alphabet.
    constexpr rank_type to_rank() const noexcept
    {
        return static_cast<rank_type>(_value);
    }
    //!\}

    /*!\name Write functions
     * \{
     */
    //!\brief Assign from a character.
    constexpr wuss<SIZE> & assign_char(char_type const c) noexcept
    {
        _value = char_to_value[c];
        return *this;
    }

    //!\brief Assign from a numeric value.
    constexpr wuss<SIZE> & assign_rank(rank_type const c)
    {
        assert(c < value_size);
        _value = static_cast<internal_type>(c);
        return *this;
    }
    //!\}

    //!\brief The size of the alphabet, i.e. the number of different values it can take.
    static constexpr rank_type value_size{SIZE};

    //!\name Comparison operators
    //!\{
    constexpr bool operator==(wuss<SIZE> const & rhs) const noexcept
    {
        return _value == rhs._value;
    }

    constexpr bool operator!=(wuss<SIZE> const & rhs) const noexcept
    {
        return _value != rhs._value;
    }

    constexpr bool operator<(wuss<SIZE> const & rhs) const noexcept
    {
        return _value < rhs._value;
    }

    constexpr bool operator>(wuss<SIZE> const & rhs) const noexcept
    {
        return _value > rhs._value;
    }

    constexpr bool operator<=(wuss<SIZE> const & rhs) const noexcept
    {
        return _value <= rhs._value;
    }

    constexpr bool operator>=(wuss<SIZE> const & rhs) const noexcept
    {
        return _value >= rhs._value;
    }
    //!\}

protected:
    //!\privatesection
    /*!\brief The internal type is a strictly typed enum.
     *
     * This is done to prevent aggregate initialization from numbers and/or chars.
     * It is has the drawback that it also introduces a scope which in turn makes
     * the static "letter values " members necessary.
     */
    enum struct internal_type : rank_type
    {
        NP, // not paired .
        BL, // bracket left <
        BR, // bracket right >
        NP1, // not paired :
        NP2, // not paired ,
        NP3, // not paired _
        NP4, // not paired ~
        BL1, // bracket left (
        BR1, // bracket right )
        BL2, // bracket left [
        BR2, // bracket right ]
        BL3, // bracket left {
        BR3, // bracket right }
        PKLa, PKRa, PKLb, PKRb, PKLc, PKRc, PKLd, PKRd, // pseudoknots A,a,B,b,...
        PKLe, PKRe, PKLf, PKRf, PKLg, PKRg, PKLh, PKRh,
        PKLi, PKRi, PKLj, PKRj, PKLk, PKRk, PKLl, PKRl,
        PKLm, PKRm, PKLn, PKRn, PKLo, PKRo, PKLp, PKRp,
        PKLq, PKRq, PKLr, PKRr, PKLs, PKRs, PKLt, PKRt,
        PKLu, PKRu, PKLv, PKRv, PKLw, PKRw, PKLx, PKRx,
        PKLy, PKRy, PKLz, PKRz
    };

    //!\brief Value to char conversion table.
    static constexpr std::array<char_type, value_size> value_to_char
    {
        [] () constexpr
        {
            using in_t = internal_type;
            std::array<char_type , value_size> ret{};

            // canonical
            ret[static_cast<rank_type>(in_t::NP)] = '.';
            ret[static_cast<rank_type>(in_t::BL)] = '<';
            ret[static_cast<rank_type>(in_t::BR)] = '>';
            ret[static_cast<rank_type>(in_t::NP1)] = ':';
            ret[static_cast<rank_type>(in_t::NP2)] = ',';
            ret[static_cast<rank_type>(in_t::NP3)] = '_';
            ret[static_cast<rank_type>(in_t::NP4)] = '~';
            ret[static_cast<rank_type>(in_t::BL1)] = '(';
            ret[static_cast<rank_type>(in_t::BR1)] = ')';
            ret[static_cast<rank_type>(in_t::BL2)] = '[';
            ret[static_cast<rank_type>(in_t::BR2)] = ']';
            ret[static_cast<rank_type>(in_t::BL3)] = '{';
            ret[static_cast<rank_type>(in_t::BR3)] = '}';

            // pseudoknot letters
            for (char_type c = 'A'; 2 * (c - 'A') + static_cast<rank_type>(in_t::PKLa) < value_size; ++c)
                ret[2 * (c - 'A') + static_cast<rank_type>(in_t::PKLa)] = c;
            for (char_type c = 'a'; 2 * (c - 'a') + static_cast<rank_type>(in_t::PKRa) < value_size; ++c)
                ret[2 * (c - 'a') + static_cast<rank_type>(in_t::PKRa)] = c;
            return ret;
        } ()
    };

    //!\brief Char to value conversion table.
    static constexpr std::array<internal_type, 256> char_to_value
    {
        [] () constexpr
        {
            using in_t = internal_type;
            std::array<in_t, 256> ret{};

            // initialize with unpaired outside structure (std::array::fill unfortunately not constexpr)
            for (auto & c : ret)
                c = in_t::NP1;

            // set alphabet values
            for (rank_type val = 0; val < value_size; ++val)
                ret[value_to_char[val]] = static_cast<in_t>(val);
            return ret;
        } ()
    };

public:
    //!\privatesection
    //!\brief The data member.
    internal_type _value;
    //!\publicsection
};

template <uint8_t SIZE>
constexpr wuss<SIZE> wuss<SIZE>::NP{internal_type::NP};

template <uint8_t SIZE>
constexpr wuss<SIZE> wuss<SIZE>::BL{internal_type::BL};

template <uint8_t SIZE>
constexpr wuss<SIZE> wuss<SIZE>::BR{internal_type::BR};

template <uint8_t SIZE>
constexpr wuss<SIZE> wuss<SIZE>::NP1{internal_type::NP1};

template <uint8_t SIZE>
constexpr wuss<SIZE> wuss<SIZE>::NP2{internal_type::NP2};

template <uint8_t SIZE>
constexpr wuss<SIZE> wuss<SIZE>::NP3{internal_type::NP3};

template <uint8_t SIZE>
constexpr wuss<SIZE> wuss<SIZE>::NP4{internal_type::NP4};

template <uint8_t SIZE>
constexpr wuss<SIZE> wuss<SIZE>::BL1{internal_type::BL1};

template <uint8_t SIZE>
constexpr wuss<SIZE> wuss<SIZE>::BR1{internal_type::BR1};

template <uint8_t SIZE>
constexpr wuss<SIZE> wuss<SIZE>::BL2{internal_type::BL2};

template <uint8_t SIZE>
constexpr wuss<SIZE> wuss<SIZE>::BR2{internal_type::BR2};

template <uint8_t SIZE>
constexpr wuss<SIZE> wuss<SIZE>::BL3{internal_type::BL3};

template <uint8_t SIZE>
constexpr wuss<SIZE> wuss<SIZE>::BR3{internal_type::BR3};

} // namespace seqan3

namespace seqan3::detail
{

//!\brief seqan3::wuss is defined as being a structure alphabet.
//!\ingroup structure
template <uint8_t SIZE>
struct is_structure<wuss<SIZE>> : public std::true_type
{};

} // namespace seqan3::detail

#ifndef NDEBUG
static_assert(seqan3::alphabet_concept<seqan3::wuss<>>);
static_assert(seqan3::structure_concept<seqan3::wuss<>>);
#endif

// ------------------------------------------------------------------
// containers
// ------------------------------------------------------------------

namespace seqan3
{

//!\brief Alias for an std::vector of seqan3::wuss.
//!\relates wuss
using wuss_vector = std::vector<wuss<>>;


/*!\brief Alias for an std::basic_string of seqan3::wuss.
 * \relates wuss
 *
 * \attention
 * Note that we recommend using seqan3::wuss_vector instead of wuss_string in almost all situations.
 * While the C++ style operations on the string are well supported, you should not access the internal c-string
 * and should not use C-Style operations on it, e.g. the `char_traits::strlen` function will not return the
 * correct length of the string (while the `.size()` returns the correct value).
 */
using wuss_string = std::basic_string<wuss<>, std::char_traits<wuss<>>>;

} // namespace seqan3

// ------------------------------------------------------------------
// literals
// ------------------------------------------------------------------

namespace seqan3::literal
{

/*!\brief wuss literal
 * \relates seqan3::wuss
 * \returns seqan3::wuss_vector
 *
 * You can use this string literal to easily assign to wuss_vector:
 *
 *~~~~~~~~~~~~~~~{.cpp}
 *     using namespace seqan3::literal;
 *     wuss_vector foo{".<..>."_wuss};
 *     wuss_vector bar = ".<..>."_wuss;
 *     auto bax = ".<..>."_wuss;
 *~~~~~~~~~~~~~~~
 *
 * \attention
 * All seqan3 literals are in the namespace seqan3::literal!
 */

inline wuss_vector operator "" _wuss(const char * s, std::size_t n)
{
    wuss_vector r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}

/*!\brief wuss string literal
 * \relates seqan3::wuss
 * \returns seqan3::wuss_string
 *
 * You can use this string literal to easily assign to wuss_vector:
 *
 *~~~~~~~~~~~~~~~{.cpp}
 *     using namespace seqan3::literal;
 *     wuss_string foo{".<..>."_wusss};
 *     wuss_string bar = ".<..>."_wusss;
 *     auto bax = ".<..>."_wusss;
 *~~~~~~~~~~~~~~~
 *
 * Please note the limitations of seqan3::wuss_string and consider using the \link operator""_wuss \endlink instead.
 *
 * \attention
 * All seqan3 literals are in the namespace seqan3::literal!
 */

inline wuss_string operator "" _wusss(const char * s, std::size_t n)
{
    wuss_string r;
    r.resize(n);

    for (size_t i = 0; i < n; ++i)
        r[i].assign_char(s[i]);

    return r;
}

} // namespace seqan3::literal
