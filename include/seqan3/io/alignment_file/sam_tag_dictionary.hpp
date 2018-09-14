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
 * \brief Provides the seqan3::sam_tag_dictionary class and auxiliaries.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include <map>
#include <variant>

#include <seqan3/core/concept/core_language.hpp>
#include <seqan3/core/metafunction/template_inspection.hpp>
#include <seqan3/io/stream/parse_condition.hpp>
#include <seqan3/std/concepts>

namespace seqan3::detail
{
    //!\brief std::variant of allowed types for optional tag fields of the SAM format.
    //!\ingroup alignment_file_io
    using sam_tag_variant = std::variant<char, int32_t, float, std::string,
                                         std::vector<uint8_t>, std::vector<int8_t>,
                                         std::vector<uint16_t>, std::vector<int16_t>,
                                         std::vector<uint32_t>, std::vector<int32_t>,
                                         std::vector<float>>;
}

namespace seqan3
{

/*!\brief The SAM tag literal, such that tags can be used in constant expressions.
 * \ingroup alignment_file_io
 * \tparam char_t The char type. Usually `char`. Parameter pack `...s` must be of
 *                length 2, since SAM tags consist of two letters (char0 and char1).
 * \returns The unique identifier of the SAM tag computed by char0 * 128 + char1.
 *
 * \details
 *
 * A SAM tag consists of two letters, initialized via the string literal ""_tag,
 * which delegate to its unique id.
 * e.g.
 *
 * \snippet test/snippet/io/alignment_file/sam_tag_dictionary/sam_tag_dictionary.cpp tag
 *
 * The purpose of those tags is to fill or query the seqan3::sam_tag_dictionary
 * for a specific key (tag_id) and retrieve the corresponding value.
 *
 * \sa seqan3::sam_tag_dictionary
 */
template <typename char_t, char_t ...s>
constexpr uint16_t operator""_tag()
{
    static_assert(std::Same<char_t, char>, "Illegal SAM tag: Type must be char.");
    static_assert((sizeof...(s) == 2),     "Illegal SAM tag: Exactly two characters must be given.");

    char constexpr str[sizeof...(s)] = {s...};
    char constexpr char0 = str[0];
    char constexpr char1 = str[1];

    static_assert((is_alpha(char0) && is_alnum(char1)),
                  "Illegal SAM tag: a SAM tag must match /[A-Za-z][A-Za-z0-9]/.");

    return static_cast<uint16_t>(char0) * 256 + static_cast<uint16_t>(char1);
}

/*!\brief The generic base class.
 * \ingroup alignment_file_io
 *
 * \attention This is a pure base class that needs to be specialized in order to
 *            be used.
 *
 * ### How to specialize the type for your custom tag
 *
 * All known tags of the SAM specifications already have a pre-defined type.
 * If you want to specify the type of your custom tag (the SAM specifications
 * recommend to use X?, Y? or Z?) you need to overload the seqan3::sam_tag_type
 * struct in the following way: (take tag "XX" as an example)
 *
 * \snippet test/snippet/io/alignment_file/sam_tag_dictionary/sam_tag_dictionary.cpp type_overload
 *
 * Everything else, like the get and set functions and correct SAM output
 * (XX:i:? in this case) is handled by the seqan3::sam_tag_dictionary.
 *
 * The seqan3::sam_tag_type is overloaded the following SAM tags:
 *
 * | Tag Name | SeqAn3 Type Implementation |
 * | -------- | --------------------- |
 * | "AM"_tag | int32_t               |
 * | "AS"_tag | int32_t               |
 * | "BC"_tag | std::string           |
 * | "BQ"_tag | std::string           |
 * | "BZ"_tag | std::string           |
 * | "CB"_tag | std::string           |
 * | "CC"_tag | std::string           |
 * | "CG"_tag | std::vector<int32_t>  |
 * | "CM"_tag | int32_t               |
 * | "CO"_tag | std::string           |
 * | "CP"_tag | int32_t               |
 * | "CQ"_tag | std::string           |
 * | "CR"_tag | std::string           |
 * | "CS"_tag | std::string           |
 * | "CT"_tag | std::string           |
 * | "CY"_tag | std::string           |
 * | "E2"_tag | std::string           |
 * | "FI"_tag | int32_t               |
 * | "FS"_tag | std::string           |
 * | "FZ"_tag | std::vector<uint16_t> |
 * | "H0"_tag | int32_t               |
 * | "H1"_tag | int32_t               |
 * | "H2"_tag | int32_t               |
 * | "HI"_tag | int32_t               |
 * | "IH"_tag | int32_t               |
 * | "LB"_tag | std::string           |
 * | "MC"_tag | std::string           |
 * | "MD"_tag | std::string           |
 * | "MI"_tag | std::string           |
 * | "MQ"_tag | int32_t               |
 * | "NH"_tag | int32_t               |
 * | "NM"_tag | int32_t               |
 * | "OC"_tag | std::string           |
 * | "OP"_tag | int32_t               |
 * | "OQ"_tag | std::string           |
 * | "OX"_tag | std::string           |
 * | "PG"_tag | std::string           |
 * | "PQ"_tag | int32_t               |
 * | "PT"_tag | std::string           |
 * | "PU"_tag | std::string           |
 * | "Q2"_tag | std::string           |
 * | "QT"_tag | std::string           |
 * | "QX"_tag | std::string           |
 * | "R2"_tag | std::string           |
 * | "RG"_tag | std::string           |
 * | "RT"_tag | std::string           |
 * | "RX"_tag | std::string           |
 * | "SA"_tag | std::string           |
 * | "SM"_tag | int32_t               |
 * | "TC"_tag | int32_t               |
 * | "U2"_tag | std::string           |
 * | "UQ"_tag | int32_t               |
 *
 */
template <uint16_t tag_value>
struct sam_tag_type
{
    //!\brief The type for all unknown tags with no extra overload defaults to an std::variant.
    using type = detail::sam_tag_variant;
};

//!\brief Short cut helper for seqan3::sam_tag_type::type.
//!\ingroup alignment_file_io
template <uint16_t tag_value>
using sam_tag_type_t = typename sam_tag_type<tag_value>::type;

//!\cond
template <> struct sam_tag_type<"AM"_tag> { using type = int32_t; };
template <> struct sam_tag_type<"AS"_tag> { using type = int32_t; };
template <> struct sam_tag_type<"BC"_tag> { using type = std::string; };
template <> struct sam_tag_type<"BQ"_tag> { using type = std::string; };
template <> struct sam_tag_type<"BZ"_tag> { using type = std::string; };
template <> struct sam_tag_type<"CB"_tag> { using type = std::string; };
template <> struct sam_tag_type<"CC"_tag> { using type = std::string; };
template <> struct sam_tag_type<"CG"_tag> { using type = std::vector<int32_t>; };
template <> struct sam_tag_type<"CM"_tag> { using type = int32_t; };
template <> struct sam_tag_type<"CO"_tag> { using type = std::string; };
template <> struct sam_tag_type<"CP"_tag> { using type = int32_t; };
template <> struct sam_tag_type<"CQ"_tag> { using type = std::string; };
template <> struct sam_tag_type<"CR"_tag> { using type = std::string; };
template <> struct sam_tag_type<"CS"_tag> { using type = std::string; };
template <> struct sam_tag_type<"CT"_tag> { using type = std::string; };
template <> struct sam_tag_type<"CY"_tag> { using type = std::string; };
template <> struct sam_tag_type<"E2"_tag> { using type = std::string; };
template <> struct sam_tag_type<"FI"_tag> { using type = int32_t; };
template <> struct sam_tag_type<"FS"_tag> { using type = std::string; };
template <> struct sam_tag_type<"FZ"_tag> { using type = std::vector<uint16_t>; };

// template <> struct sam_tag_type<"GC"_tag> {};
// template <> struct sam_tag_type<"GQ"_tag> {};
// template <> struct sam_tag_type<"GS"_tag> {};

template <> struct sam_tag_type<"H0"_tag> { using type = int32_t; };
template <> struct sam_tag_type<"H1"_tag> { using type = int32_t; };
template <> struct sam_tag_type<"H2"_tag> { using type = int32_t; };
template <> struct sam_tag_type<"HI"_tag> { using type = int32_t; };
template <> struct sam_tag_type<"IH"_tag> { using type = int32_t; };
template <> struct sam_tag_type<"LB"_tag> { using type = std::string; };
template <> struct sam_tag_type<"MC"_tag> { using type = std::string; };
template <> struct sam_tag_type<"MD"_tag> { using type = std::string; };

// template <> struct sam_tag_type<"MF"_tag> {};

template <> struct sam_tag_type<"MI"_tag> { using type = std::string; };
template <> struct sam_tag_type<"MQ"_tag> { using type = int32_t; };
template <> struct sam_tag_type<"NH"_tag> { using type = int32_t; };
template <> struct sam_tag_type<"NM"_tag> { using type = int32_t; };
template <> struct sam_tag_type<"OC"_tag> { using type = std::string; };
template <> struct sam_tag_type<"OP"_tag> { using type = int32_t; };
template <> struct sam_tag_type<"OQ"_tag> { using type = std::string; };
template <> struct sam_tag_type<"OX"_tag> { using type = std::string; };
template <> struct sam_tag_type<"PG"_tag> { using type = std::string; };
template <> struct sam_tag_type<"PQ"_tag> { using type = int32_t; };
template <> struct sam_tag_type<"PT"_tag> { using type = std::string; };
template <> struct sam_tag_type<"PU"_tag> { using type = std::string; };
template <> struct sam_tag_type<"Q2"_tag> { using type = std::string; };
template <> struct sam_tag_type<"QT"_tag> { using type = std::string; };
template <> struct sam_tag_type<"QX"_tag> { using type = std::string; };
template <> struct sam_tag_type<"R2"_tag> { using type = std::string; };
template <> struct sam_tag_type<"RG"_tag> { using type = std::string; };
template <> struct sam_tag_type<"RT"_tag> { using type = std::string; };
template <> struct sam_tag_type<"RX"_tag> { using type = std::string; };

// template <> struct sam_tag_type<"S2"_tag> {};

template <> struct sam_tag_type<"SA"_tag> { using type = std::string; };
template <> struct sam_tag_type<"SM"_tag> { using type = int32_t; };

// template <> struct sam_tag_type<"SQ"_tag> {};

template <> struct sam_tag_type<"TC"_tag> { using type = int32_t; };
template <> struct sam_tag_type<"U2"_tag> { using type = std::string; };
template <> struct sam_tag_type<"UQ"_tag> { using type = int32_t; };
//!\endcond

/*!\brief The SAM tag dictionary class that stores all optional SAM fields.
 * \ingroup alignment_file_io
 *
 * \details
 *
 * ### SAM tags
 *
 * A SAM tag consists of two letters, initialized via the string literal ""_tag,
 * which delegates to its unique id (type uint16_t).
 * Example:
 *
 * \snippet test/snippet/io/alignment_file/sam_tag_dictionary/sam_tag_dictionary.cpp tag
 *
 * The purpose of those tags is to fill or query the seqan3::sam_tag_dictionary
 * for a specific key (tag_id) and retrieve the corresponding value.
 *
 * ### SAM tag types
 *
 * Note that a SAM tag is always associated with a specific type.
 * In the SAM format, the type is indicated in the second argument of the
 * TAG:TYPE:VALUE field. For example "NM:i:3" specifies the NM tag of an integer
 * type with value 3.
 * In seqan3, the types for
 * [known](https://samtools.github.io/hts-specs/SAMtags.pdf) SAM tags
 * are pre-defined by a metafunction called seqan3::sam_tag_type. You can access
 * the type via:
 *
 * \snippet test/snippet/io/alignment_file/sam_tag_dictionary/sam_tag_dictionary.cpp tag_type_t
 *
 * which is the short cut for:
 *
 * \snippet test/snippet/io/alignment_file/sam_tag_dictionary/sam_tag_dictionary.cpp tag_type
 *
 * The following types are allowed by the
 * [SAM specifications](https://samtools.github.io/hts-specs/SAMtags.pdf):
 *
 * |Type | Regexp matching VALUE                  | Description                             | SeqAn Type           |
 * |-----|----------------------------------------|-----------------------------------------|----------------------|
 * | A   | [!-~]                                  |  Printable character                    | char                 |
 * | i   | [-+]?[0-9]+                            |  Signed integer                         | int32_t              |
 * | f   | [-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)? |  Single-precision floating number       | float                |
 * | Z   | [ !-~]*                                |  Printable string, including space      | std::string          |
 * | H   | ([0-9A-F][0-9A-F])*                    |  Byte array in the Hex format           | std::vector<uint8_t> |
 * | B   | [cCsSiIf]\(,[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?)+ |  Integer or numeric array | std::vector<T>       |
 *
 * For  an  integer  or  numeric  array  (type  ‘B’),  the second  letter can be
 * one of ‘cCsSiIf’, corresponding to type **T** = int8_t, uint8_t, int16_t,
 * uint16_t, int32_t, uint32_t and float, respectively.
 *
 * ### Using the sam_tag_dictionary
 *
 * The dictionary can be accessed via the functions seqan3::sam_tag_dictionary::get() and
 * seqan3::sam_tag_dictionary::set(). Every time the SAM tag you wish to query
 * for must be given as a template argument to the functions.
 *
 * Example:
 *
 * \snippet test/snippet/io/alignment_file/sam_tag_dictionary/general_usage.cpp all
 *
 * \attention You can get any SAM_tag out of the dictionary, even if the tag is
 *            user defined, but note that for unkown tags the return type is an
 *            [std::variant](https://en.cppreference.com/w/cpp/utility/variant).
 *            If you want specify the return type of your custom tag, you need
 *            to overload the seqan3::sam_tag_type metafunction.
 *
 * Unknown Tag Example:
 *
 * \snippet test/snippet/io/alignment_file/sam_tag_dictionary/unknown_tag.cpp all
 *
 * As mentioned before you can either overload the metafunction seqan3::sam_tag_type
 * for the tag "XZ" or learn more about an std::variant at
 * https://en.cppreference.com/w/cpp/utility/variant.
 *
 * \sa seqan3::sam_tag_type
 * \sa https://en.cppreference.com/w/cpp/utility/variant
 * \sa https://samtools.github.io/hts-specs/SAMv1.pdf
 * \sa https://samtools.github.io/hts-specs/SAMtags.pdf
 */
class sam_tag_dictionary : public std::map<uint16_t, detail::sam_tag_variant>
{
private:
    //!\brief The base type.
    using base_type = std::map<uint16_t, detail::sam_tag_variant>;

public:
    //!\brief The variant type defining all valid SAM tag field types.
    using variant_type = detail::sam_tag_variant;

    /*!\name Getter function for the seqan3::sam_tag_dictionary.
     *\brief Gets the value of known SAM tags by its correct type instead of the std::variant.
     * \tparam tag The unique tag id of a SAM tag.
     * \returns The value corresponding to the key `tag` of type seqan3::sam_tag_type<tag>::type.
     *
     * \details
     *
     * See the seqan3::sam_tag_dictionary detailed documentation below for an example.
     *
     * \attention This function is only available for tags that have an
     *            seqan3::sam_tag_type<tag>::type overload. See the metafunction
     *            documentation for further details.
     * \{
     */

    //!\brief Uses std::map::operator[] for access and default initializes new keys.
    template <uint16_t tag>
    //!\cond
        requires !std::Same<sam_tag_type_t<tag>, variant_type>
    //!\endcond
    auto & get() &
    {
        if ((*this).count(tag) == 0)
            (*this)[tag] = sam_tag_type_t<tag>{}; // set correct type if tag is not set yet on

        return std::get<sam_tag_type_t<tag>>((*this)[tag]);
    }

    //!\brief Uses std::map::operator[] for access and default initializes new keys.
    template <uint16_t tag>
    //!\cond
        requires !std::Same<sam_tag_type_t<tag>, variant_type>
    //!\endcond
    auto && get() &&
    {
        if ((*this).count(tag) == 0)
            (*this)[tag] = sam_tag_type_t<tag>{}; // set correct type if tag is not set yet on

        return std::get<sam_tag_type_t<tag>>(std::move((*this)[tag]));
    }

    //!\brief Uses std::map::at() for access and throws when the key is unknown.
    //!\throws std::out_of_range if map has no key `tag`.
    template <uint16_t tag>
    //!\cond
        requires !std::Same<sam_tag_type_t<tag>, variant_type>
    //!\endcond
    auto const & get() const &
    {
        return std::get<sam_tag_type_t<tag>>((*this).at(tag));
    }

    //!\brief Uses std::map::at() for access and throws when the key is unknown.
    //!\throws std::out_of_range if map has no key `tag`.
    template <uint16_t tag>
    //!\cond
        requires !std::Same<sam_tag_type_t<tag>, variant_type>
    //!\endcond
    auto const && get() const &&
    {
        return std::get<sam_tag_type_t<tag>>(std::move((*this).at(tag)));
    }
    //!\}
};

} // namespace seqan3
