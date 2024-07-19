// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::debug_stream and related types.
 */

#pragma once

#include <ranges>

#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/adaptation/uint.hpp>
#include <seqan3/alphabet/range/sequence.hpp>
#include <seqan3/core/debug_stream/debug_stream_type.hpp>
#include <seqan3/core/range/type_traits.hpp>

namespace seqan3::detail
{
/*!\brief A helper concept definition for ranges that can be streamed to the seqan3::debug_stream.
 * \tparam rng_t The range type to check.
 * \ingroup core_debug_stream
 *
 * \details
 *
 * This concept refines the std::ranges::input_range concept to allow streaming the range object to the debug stream,
 * with the following requirements:
 *
 * * `rng_t` is not the same type as `std::ranges::range_reference_t<rng_t>`,
 * * `rng_t` is not a pointer or c-style array,
 * * `std::ranges::range_reference_t<rng_t>` is not `char`.
 */
template <typename rng_t>
concept debug_stream_range_guard =
    !std::same_as<std::remove_cvref_t<std::ranges::range_reference_t<rng_t>>, std::remove_cvref_t<rng_t>>
    && !(std::is_pointer_v<std::decay_t<rng_t>>
         && std::same_as<std::remove_cvref_t<std::ranges::range_reference_t<rng_t>>, char>);
// prevent recursive instantiation: std::same_as<std::remove_cvref_t<std::ranges::range_reference_t<rng_t>>, char>
// exclude null-terminated strings: std::is_pointer_v<std::decay_t<rng_t>>

/*!\brief Helper template variable that checks if the reference type of a range can be streamed into an instance of
 *        seqan3::debug_stream_type .
 * \tparam rng_t The range type to check.
 * \tparam char_t The char type of the stream.
 * \ingroup core_debug_stream
 *
 * \details
 *
 * Evaluates to `true` if the following expression is valid: `debug_stream << *rng.begin();`, where rng is of  type
 * rng_t. Otherwise false.
 */
template <std::ranges::range rng_t, typename char_t>
constexpr bool reference_type_is_streamable_v = false;

//!\cond
template <std::ranges::range rng_t, typename char_t>
    requires requires (std::ranges::range_reference_t<rng_t> l, debug_stream_type<char_t> s) {
        { s << l };
    }
constexpr bool reference_type_is_streamable_v<rng_t, char_t> = true;
//!\endcond
} // namespace seqan3::detail

namespace seqan3
{
/*!\name Formatted output overloads
 * \{
 */
/*!\brief All input ranges can be printed to the seqan3::debug_stream element-wise (if their elements are printable).
 * \tparam rng_t Type of the range to be printed; must model std::ranges::input_range.
 * \param s The seqan3::debug_stream.
 * \param r The input range.
 * \relates seqan3::debug_stream_type
 *
 * \details
 *
 * If the element type models seqan3::alphabet (and is not an unsigned integer), the range is printed
 * just as if it were a string, i.e. <tt>std::vector<dna4>{'C'_dna4, 'G'_dna4, 'A'_dna4}</tt> is printed as "CGA".
 *
 * In all other cases the elements are comma separated and the range is enclosed in brackets, i.e.
 * `std::vector<int>{3, 1, 33, 7}` is printed as "[3,1,33,7]".
 *
 * \if DEV
 * Note that overloads for range based streaming need to refine the seqan3::detail::debug_stream_range_guard concept
 * to avoid ambiguous function calls.
 * \endif
 */
template <typename char_t, std::ranges::input_range rng_t>
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & s, rng_t && r)
    requires detail::debug_stream_range_guard<rng_t>
{
    static_assert(detail::reference_type_is_streamable_v<rng_t, char_t>,
                  "The reference type of the passed range cannot be streamed into the debug_stream.");

    s << '[';
    auto b = std::ranges::begin(r);
    auto e = std::ranges::end(r);
    if (b != e)
    {
        s << *b;
        ++b;
    }
    while (b != e)
    {
        s << ',';
        s << *b;
        ++b;
    }
    s << ']';

    return s;
}

/*!\brief All biological sequences can be printed to the seqan3::debug_stream.
 * \tparam sequence_t Type of the (biological) sequence to be printed; must model seqan3::sequence.
 * \param s The seqan3::debug_stream.
 * \param sequence The input range.
 * \relates seqan3::debug_stream_type
 *
 * \details
 *
 * The (biological) sequence (except for ranges over unsigned integers) is printed just as if it were a string, i.e.
 * <tt>std::vector<dna4>{'C'_dna4, 'G'_dna4, 'A'_dna4}</tt> is printed as "CGA".
 *
 * \if DEV
 * Note that overloads for range based streaming need to refine the seqan3::detail::debug_stream_range_guard concept
 * to avoid ambiguous function calls.
 * \endif
 */
template <typename char_t, sequence sequence_t>
inline debug_stream_type<char_t> & operator<<(debug_stream_type<char_t> & s, sequence_t && sequence)
    requires detail::debug_stream_range_guard<sequence_t>
          && (!detail::is_uint_adaptation_v<std::remove_cvref_t<std::ranges::range_reference_t<sequence_t>>>)
{
    for (auto && chr : sequence)
        s << chr;
    return s;
}

//!\}

} // namespace seqan3
