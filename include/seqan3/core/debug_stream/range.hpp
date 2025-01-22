// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 * \brief Provides seqan3::debug_stream and related types.
 */

#pragma once

#include <functional>
#include <ranges>

#include <seqan3/alphabet/adaptation/char.hpp>
#include <seqan3/alphabet/adaptation/uint.hpp>
#include <seqan3/alphabet/range/sequence.hpp>
#include <seqan3/core/debug_stream/debug_stream_type.hpp>
#include <seqan3/core/range/type_traits.hpp>

namespace seqan3
{

/*!
 * \interface seqan3::nonrecursive_range <>
 * \brief A concept that checks whether a range is non-recursive.
 *
 * A range is considered non-recursive if it is itself a range whose reference type is not the range type itself.
 *
 * \tparam rng_t The type to check.
 * \ingroup core_debug_stream
 */
template <typename rng_t>
concept nonrecursive_range =
    !std::same_as<std::remove_cvref_t<std::ranges::range_reference_t<rng_t>>, std::remove_cvref_t<rng_t>>;

/*!\brief A printer for arbitrary input ranges.
 *
 * All input ranges can be printed to the seqan3::debug_stream element-wise (if their elements are printable).
 * If the element type models seqan3::alphabet (and is not an unsigned integer), the range is printed
 * just as if it were a string, i.e. <tt>std::vector<dna4>{'C'_dna4, 'G'_dna4, 'A'_dna4}</tt> is printed as "CGA".
 *
 * In all other cases the elements are comma separated and the range is enclosed in brackets, i.e.
 * `std::vector<int>{3, 1, 33, 7}` is printed as "[3,1,33,7]".
 *
 * This printer excludes recursive ranges, such as entities from the std::filesystem library.
 *
 * \tparam rng_t Type of the range to be printed; must model std::ranges::input_range and be non-recursive.
 * \ingroup core_debug_stream
 */
template <typename rng_t>
    requires std::ranges::input_range<rng_t> && nonrecursive_range<rng_t>
struct input_range_printer<rng_t>
{
    /*!\brief Prints the elements of a sequence to an output stream.
     *
     * \tparam stream_t The type of the stream.
     * \tparam arg_t The type of the argument.
     *
     * \param[in,out] stream The output stream to print to.
     * \param[in] arg The range to be printed.
     */
    template <typename stream_t, typename arg_t>
    constexpr void operator()(stream_t & stream, arg_t && arg) const
    {
        stream << '[';
        auto first = std::ranges::begin(arg);
        auto last = std::ranges::end(arg);
        if (first != last)
        {
            stream << *first;
            ++first;
        }
        while (first != last)
        {
            stream << ',';
            stream << *first;
            ++first;
        }
        stream << ']';
    }
};

/*!\brief A printer for (biological) sequences.
 *
 * The (biological) sequence (except for ranges over unsigned integers) is printed just as if it were a string, i.e.
 * <tt>std::vector<dna4>{'C'_dna4, 'G'_dna4, 'A'_dna4}</tt> is printed as "CGA".
 *
 * \tparam sequence_t The type of the sequence to be printed; must model seqan3::sequence.
 * \ingroup core_debug_stream
 */
template <sequence sequence_t>
struct sequence_printer<sequence_t>
{
    /*!\brief Prints the elements of a sequence to an output stream.
     *
     * \tparam stream_t The type of the stream.
     * \tparam arg_t The type of the argument.
     *
     * \param[in,out] stream The output stream to print to.
     * \param[in] arg The sequence to be printed.
     */
    template <typename stream_t, typename arg_t>
    constexpr void operator()(stream_t & stream, arg_t && arg) const
    {
        for (auto && chr : arg)
            stream << chr;
    }
};

/*!\brief A printer for character sequences.
 * \ingroup core_debug_stream
 *
 * This struct provides a printer for character sequences. It is used to print character sequences to a stream.
 * The character sequence must be an input range and the range reference type must be a character type, i.e.
 * seqan3::detail::is_char_adaptation_v evaluates to `true`.
 *
 * \tparam char_sequence_t The type of the character sequence.
 */
template <typename char_sequence_t>
    requires std::ranges::input_range<char_sequence_t>
          && (detail::is_char_adaptation_v<std::remove_cvref_t<std::ranges::range_reference_t<char_sequence_t>>>)
struct char_sequence_printer<char_sequence_t>
{
    /*!\brief Prints the character sequence to the given stream.
     *
     * \tparam stream_t The type of the stream.
     * \tparam arg_t The type of the argument.
     * \param stream The stream to print to.
     * \param arg The character sequence to print.
     */
    template <typename stream_t, typename arg_t>
    constexpr void operator()(stream_t & stream, arg_t && arg) const
    {
        // null-terminated string
        if constexpr (std::is_pointer_v<std::decay_t<char_sequence_t>>)
            return std::invoke(std_printer<char_sequence_t>{}, stream, std::forward<arg_t>(arg));

        for (auto && chr : arg)
            stream << chr;
    }
};

/*!\brief A printer for integer sequences.
 *
 * This struct provides a printer for integer sequences.
 * The integer sequence must be an input range and the range reference type must model std::integral concept.
 *
 * \tparam integer_sequence_t The type of the integer sequence.
 * \ingroup core_debug_stream
 */
template <typename integer_sequence_t>
    requires std::ranges::input_range<integer_sequence_t>
          && std::integral<std::remove_cvref_t<std::ranges::range_reference_t<integer_sequence_t>>>
struct integer_sequence_printer<integer_sequence_t> : public input_range_printer<integer_sequence_t>
{};

} // namespace seqan3
