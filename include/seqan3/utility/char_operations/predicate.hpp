// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides character predicates for tokenisation.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <seqan3/utility/char_operations/predicate_detail.hpp>

// ----------------------------------------------------------------------------
// General Purpose Char predicates
// ----------------------------------------------------------------------------

namespace seqan3
{

/*!\name Char predicates
 * \brief A set of function objects to check if a character from an input source fulfills certain characteristics.
 * \ingroup utility_char_operations
 * \{
 */

/*!\brief Checks whether a given letter is in the specified interval.
 * \tparam interval_first The first character for which to return true.
 * \tparam interval_last  The last character (inclusive) for which to return true.
 * \ingroup utility_char_operations
 *
 * \details
 *
 * This function like object returns true for all characters in the given range, false otherwise.
 *
 * ### Example
 *
 * \include test/snippet/utility/char_operations/char_predicate_is_in_interval.cpp
 * \hideinitializer
 */
template <uint8_t interval_first, uint8_t interval_last>
    requires (interval_first <= interval_last)
inline constexpr auto is_in_interval = detail::is_in_interval_type<interval_first, interval_last>{};

/*!\brief Checks whether a given letter is the same as the template non-type argument.
 * \tparam char_v The letter to compare against.
 * \ingroup utility_char_operations
 *
 * \details
 *
 * This function like object returns true if the argument is the same as the template argument, false otherwise.
 *
 * ### Example
 *
 * \include test/snippet/utility/char_operations/char_predicate_is_char.cpp
 * \hideinitializer
 */
template <int char_v>
inline constexpr auto is_char = detail::is_char_type<char_v>{};

/*!\brief Checks whether a given letter is equal to the EOF constant defined in `<cstdio>`.
 *
 * \details
 *
 * This function like object returns true if the argument is equal to EOF, false otherwise.
 *
 * ### Example
 *
 * \snippet test/snippet/utility/char_operations/char_predicate.cpp is_eof
 */
inline constexpr auto is_eof = is_char<EOF>;

/*!\brief Checks whether `c` is a control character.
 * \ingroup utility_char_operations
 *
 * \details
 *
 * This function like object can be used to check if a character `c` is a control character.
 * For the standard ASCII character set, control characters are those between ASCII codes 0x00 (NUL) and 0x1f (US) and
 * 0x7f (DEL).
 *
 * ### Example
 *
 * \snippet test/snippet/utility/char_operations/char_predicate.cpp is_cntrl
 */
inline constexpr auto is_cntrl = is_in_interval<'\0', static_cast<char>(31)> || is_char<static_cast<char>(127)>;

/*!\brief Checks whether `c` is a printable character.
 * \ingroup utility_char_operations
 *
 * \details
 *
 * This function like object can be used to check if a character `c` is a printable character.
 * For the standard ASCII character set, printable characters are those between ASCII codes 0x20 (space) and 0x7E (`~`).
 *
 * ### Example
 *
 * \snippet test/snippet/utility/char_operations/char_predicate.cpp is_print
 */
inline constexpr auto is_print = is_in_interval<' ', '~'>;

/*!\brief Checks whether `c` is a space character.
 * \ingroup utility_char_operations
 *
 * \details
 *
 * This function like object can be used to check if a character `c` is a space character.
 * For the standard ASCII character set, the following characters are space characters:
 *
 * * horizontal tab ('\\t')
 * * line feed ('\\n')
 * * vertical tab ('\\v')
 * * from feed ('\\f')
 * * carriage return ('\\r')
 * * space (' ')
 *
 * ### Example
 *
 * \snippet test/snippet/utility/char_operations/char_predicate.cpp is_space
 */
inline constexpr auto is_space = is_in_interval<'\t', '\r'> || is_char<' '>;

/*!\brief Checks whether `c` is a blank character.
 * \ingroup utility_char_operations
 *
 * \details
 *
 * This function like object can be used to check if a character `c` is a blank character.
 * For the standard ASCII character set, the following characters are blank characters:
 *
 * * horizontal tab ('\\t')
 * * space (' ')
 *
 * ### Example
 *
 * \snippet test/snippet/utility/char_operations/char_predicate.cpp is_blank
 */
inline constexpr auto is_blank = is_char<'\t'> || is_char<' '>;

/*!\brief Checks whether `c` is a graphic character.
 * \ingroup utility_char_operations
 *
 * \details
 *
 * This function like object can be used to check if a character `c` is a graphic (has a graphical representation)
 * character.
 * For the standard ASCII character set, the following characters are graphic characters:
 *
 * * digits (0123456789)
 * * uppercase letters (ABCDEFGHIJKLMNOPQRSTUVWXYZ)
 * * lowercase letters (abcdefghijklmnopqrstuvwxyz)
 * * punctuation characters (!"\#$%&'()*+,-./:;<=>?\@[\]^_`{|}~)
 *
 * ### Example
 *
 * \snippet test/snippet/utility/char_operations/char_predicate.cpp is_graph
 */
inline constexpr auto is_graph = is_in_interval<'!', '~'>;

/*!\brief Checks whether `c` is a punctuation character.
 * \ingroup utility_char_operations
 *
 * \details
 *
 * This function like object can be used to check if a character `c` is a punctuation character.
 * For the standard ASCII character set, the following characters are punctuation characters:
 *
 * * punctuation characters (!"\#$%&'()*+,-./:;<=>?\@[\]^_`{|}~)
 *
 * ### Example
 *
 * \snippet test/snippet/utility/char_operations/char_predicate.cpp is_punct
 */
inline constexpr auto is_punct =
    is_in_interval<'!', '/'> || is_in_interval<':', '@'> || is_in_interval<'[', '`'> || is_in_interval<'{', '~'>;

/*!\brief Checks whether `c` is a alphanumeric character.
 * \ingroup utility_char_operations
 *
 * \details
 *
 * This function like object can be used to check if a character `c` is a alphanumeric character.
 * For the standard ASCII character set, the following characters are alphanumeric characters:
 *
 * * digits (0123456789)
 * * uppercase letters (ABCDEFGHIJKLMNOPQRSTUVWXYZ)
 * * lowercase letters (abcdefghijklmnopqrstuvwxyz)
 *
 * ### Example
 *
 * \snippet test/snippet/utility/char_operations/char_predicate.cpp is_alnum
 */
inline constexpr auto is_alnum = is_in_interval<'0', '9'> || is_in_interval<'A', 'Z'> || is_in_interval<'a', 'z'>;

/*!\brief Checks whether `c` is a alphabetical character.
 * \ingroup utility_char_operations
 *
 * \details
 *
 * This function like object can be used to check if a character `c` is a alphabetical character.
 * For the standard ASCII character set, the following characters are alphabetical characters:
 *
 * * uppercase letters (ABCDEFGHIJKLMNOPQRSTUVWXYZ)
 * * lowercase letters (abcdefghijklmnopqrstuvwxyz)
 *
 * ### Example
 *
 * \snippet test/snippet/utility/char_operations/char_predicate.cpp is_alpha
 */
inline constexpr auto is_alpha = is_in_interval<'A', 'Z'> || is_in_interval<'a', 'z'>;

/*!\brief Checks whether `c` is a upper case character.
 * \ingroup utility_char_operations
 *
 * \details
 *
 * This function like object can be used to check if a character `c` is a upper case character.
 * For the standard ASCII character set, the following characters are upper case characters:
 *
 * * uppercase letters (ABCDEFGHIJKLMNOPQRSTUVWXYZ)
 *
 * ### Example
 *
 * \snippet test/snippet/utility/char_operations/char_predicate.cpp is_upper
 */
inline constexpr auto is_upper = is_in_interval<'A', 'Z'>;

/*!\brief Checks whether `c` is a lower case character.
 * \ingroup utility_char_operations
 *
 * \details
 *
 * This function like object can be used to check if a character `c` is a lower case character.
 * For the standard ASCII character set, the following characters are lower case characters:
 *
 * * lowercase letters (abcdefghijklmnopqrstuvwxyz)
 *
 * ### Example
 *
 * \snippet test/snippet/utility/char_operations/char_predicate.cpp is_lower
 */
inline constexpr auto is_lower = is_in_interval<'a', 'z'>;

/*!\brief Checks whether `c` is a digital character.
 * \ingroup utility_char_operations
 *
 * \details
 *
 * This function like object can be used to check if a character `c` is a digital character.
 * For the standard ASCII character set, the following characters are digital characters:
 *
 * * digits (0123456789)
 *
 * ### Example
 *
 * \snippet test/snippet/utility/char_operations/char_predicate.cpp is_digit
 */
inline constexpr auto is_digit = is_in_interval<'0', '9'>;

/*!\brief Checks whether `c` is a hexadecimal character.
 * \ingroup utility_char_operations
 *
 * \details
 *
 * This function like object can be used to check if a character `c` is a hexadecimal character.
 * For the standard ASCII character set, the following characters are hexadecimal characters:
 *
 * * digits (0123456789)
 * * uppercase letters (ABCDEF)
 * * lowercase letters (abcdef)
 *
 * ### Example
 *
 * \snippet test/snippet/utility/char_operations/char_predicate.cpp is_xdigit
 */
inline constexpr auto is_xdigit = is_in_interval<'0', '9'> || is_in_interval<'A', 'F'> || is_in_interval<'a', 'f'>;
//!\}

/*!\name Char predicates
 * \ingroup utility_char_operations
 * \details
 *
 * Char predicates are function like objects that can be used to check if a character `c` fulfills certain
 * constraints. SeqAn implements all predicates also available in the standard library and some more.
 *
 * ### Disjunction and Negation
 *
 * In contrast to the standard library (where the checks are implemented as functions), the functors in SeqAn can be
 * joined efficiently, maintaining constant-time evaluation independent of the number of checks. Functors can be
 * combined with the `||-operator` or negated via the `!-operator`:
 *
 * \include test/snippet/utility/char_operations/char_predicate_is_operator.cpp
 *
 * Defining complex combinations and using them in e.g. input/output can increase speed significantly over checking
 * multiple functions: we measured speed-ups of 10x for a single check and speed-ups of
 * over 20x for complex combinations.
 *
 * ### Custom predicates
 *
 * * seqan3::is_in_interval: Checks if the given character is within specified range of ASCII characters.
 * * seqan3::is_char: Checks if the character is equal to the specified ASCII character.
 * * seqan3::is_eof: Checks if a character is the end-of-file marker.
 *
 * ### Standard library predicates
 *
 * SeqAn offers the 12 predicates exactly
 * [as defined in the standard library](https://en.cppreference.com/w/cpp/string/byte) except that we have introduced
 * an underscore in the name to be consistent with our other naming.
 *
 * The following table lists the predefined character predicates and which constraints are associated with them.
 *
 *<table class="wikitable" style="background-color:#ededed;font-size:85%;text-align:center;border: 1px solid black;border-collapse: collapse">
 *
 *<tr>
 *<th colspan="3" style="border: 1px solid black"> ASCII values
 *</th>
 *<th rowspan="2" style="border: 1px solid black"> characters
 *</th>
 *<th rowspan="2" style="font-size:85%;border: 1px solid black;border: 1px solid black;">
 *<p><tt>is_cntrl</tt>
 *</p>
 *</th>
 *<th rowspan="2" style="font-size:85%;border: 1px solid black">
 *<p><tt>is_print</tt>
 *</p>
 *</th>
 *<th rowspan="2" style="font-size:85%;border: 1px solid black">
 *<p><tt>is_space</tt>
 *</p>
 *</th>
 *<th rowspan="2" style="font-size:85%;border: 1px solid black">
 *<p><tt>is_blank</tt>
 *</p>
 *</th>
 *<th rowspan="2" style="font-size:85%;border: 1px solid black">
 *<p><tt>is_graph</tt>
 *</p>
 *</th>
 *<th rowspan="2" style="font-size:85%;border: 1px solid black">
 *<p><tt>is_punct</tt>
 *</p>
 *</th>
 *<th rowspan="2" style="font-size:85%;border: 1px solid black">
 *<p><tt>is_alnum</tt>
 *</p>
 *</th>
 *<th rowspan="2" style="font-size:85%;border: 1px solid black">
 *<p><tt>is_alpha</tt>
 *</p>
 *</th>
 *<th rowspan="2" style="font-size:85%;border: 1px solid black">
 *<p><tt>is_upper</tt>
 *</p>
 *</th>
 *<th rowspan="2" style="font-size:85%;border: 1px solid black">
 *<p><tt>is_lower</tt>
 *</p>
 *</th>
 *<th rowspan="2" style="font-size:85%;border: 1px solid black">
 *<p><tt>is_digit</tt>
 *</p>
 *</th>
 *<th rowspan="2" style="font-size:85%;border: 1px solid black">
 *<p><tt>is_xdigit</tt>
 *</p>
 *</th></tr>
 *<tr>
 *<th style="border: 1px solid black"> decimal
 *</th>
 *<th style="border: 1px solid black"> hexadecimal
 *</th>
 *<th style="border: 1px solid black"> octal
 *</th></tr>
 *<tr>
 *<td style="border: 1px solid black"> 0–8
 *</td>
 *<td style="border: 1px solid black"> <code>\\x0</code>–<code>\\x8</code>
 *</td>
 *<td style="border: 1px solid black"> <code>\0</code>–<code>\10</code>
 *</td>
 *<td style="border: 1px solid black"> control codes (<code>NUL</code>, etc.)
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td></tr>
 *<tr>
 *<td style="border: 1px solid black"> 9
 *</td>
 *<td style="border: 1px solid black"> <code>\\x9</code>
 *</td>
 *<td style="border: 1px solid black"> <code>\11</code>
 *</td>
 *<td style="border: 1px solid black"> tab (<code>\\t</code>)
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td></tr>
 *<tr>
 *<td style="border: 1px solid black"> 10–13
 *</td>
 *<td style="border: 1px solid black"> <code>\\xA</code>–<code>\\xD</code>
 *</td>
 *<td style="border: 1px solid black"> <code>\12</code>–<code>\15</code>
 *</td>
 *<td style="border: 1px solid black"> whitespaces (<code>\\n</code>, <code>\\v</code>, <code>\\f</code>, <code>\\r</code>)
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td></tr>
 *<tr>
 *<td style="border: 1px solid black"> 14–31
 *</td>
 *<td style="border: 1px solid black"> <code>\\xE</code>–<code>\\x1F</code>
 *</td>
 *<td style="border: 1px solid black"> <code>\16</code>–<code>\37</code>
 *</td>
 *<td style="border: 1px solid black"> control codes
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td></tr>
 *<tr>
 *<td style="border: 1px solid black"> 32
 *</td>
 *<td style="border: 1px solid black"> <code>\\x20</code>
 *</td>
 *<td style="border: 1px solid black"> <code>\40</code>
 *</td>
 *<td style="border: 1px solid black"> space
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td></tr>
 *<tr>
 *<td style="border: 1px solid black"> 33–47
 *</td>
 *<td style="border: 1px solid black"> <code>\\x21</code>–<code>\\x2F</code>
 *</td>
 *<td style="border: 1px solid black"> <code>\41</code>–<code>\57</code>
 *</td>
 *<td style="border: 1px solid black"> <code>!\"\#$%&amp;'()*+,-./</code>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td></tr>
 *<tr>
 *<td style="border: 1px solid black"> 48–57
 *</td>
 *<td style="border: 1px solid black"> <code>\\x30</code>–<code>\\x39</code>
 *</td>
 *<td style="border: 1px solid black"> <code>\60</code>–<code>\71</code>
 *</td>
 *<td style="border: 1px solid black"> <code>0123456789</code>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td></tr>
 *<tr>
 *<td style="border: 1px solid black"> 58–64
 *</td>
 *<td style="border: 1px solid black"> <code>\\x3A</code>–<code>\\x40</code>
 *</td>
 *<td style="border: 1px solid black"> <code>\72</code>–<code>\100</code>
 *</td>
 *<td style="border: 1px solid black"> <code>:;&lt;=&gt;?\@</code>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td></tr>
 *<tr>
 *<td style="border: 1px solid black"> 65–70
 *</td>
 *<td style="border: 1px solid black"> <code>\\x41</code>–<code>\\x46</code>
 *</td>
 *<td style="border: 1px solid black"> <code>\101</code>–<code>\106</code>
 *</td>
 *<td style="border: 1px solid black"> <code>ABCDEF</code>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td></tr>
 *<tr>
 *<td style="border: 1px solid black"> 71–90
 *</td>
 *<td style="border: 1px solid black"> <code>\\x47</code>–<code>\\x5A</code>
 *</td>
 *<td style="border: 1px solid black"> <code>\107</code>–<code>\132</code>
 *</td>
 *<td style="border: 1px solid black"> <code>GHIJKLMNOP</code><br /><code>QRSTUVWXYZ</code>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td></tr>
 *<tr>
 *<td style="border: 1px solid black"> 91–96
 *</td>
 *<td style="border: 1px solid black"> <code>\\x5B</code>–<code>\\x60</code>
 *</td>
 *<td style="border: 1px solid black"> <code>\133</code>–<code>\140</code>
 *</td>
 *<td style="border: 1px solid black"> <code>[\]^_`</code>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td></tr>
 *<tr>
 *<td style="border: 1px solid black"> 97–102
 *</td>
 *<td style="border: 1px solid black"> <code>\\x61</code>–<code>\\x66</code>
 *</td>
 *<td style="border: 1px solid black"> <code>\141</code>–<code>\146</code>
 *</td>
 *<td style="border: 1px solid black"> <code>abcdef</code>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td></tr>
 *<tr>
 *<td style="border: 1px solid black"> 103–122
 *</td>
 *<td style="border: 1px solid black"> <code>\\x67</code>–<code>\\x7A</code>
 *</td>
 *<td style="border: 1px solid black"> <code>\147</code>–<code>\172</code>
 *</td>
 *<td style="border: 1px solid black"> <code>ghijklmnop</code><br /><code>qrstuvwxyz</code>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td></tr>
 *<tr>
 *<td style="border: 1px solid black"> 123–126
 *</td>
 *<td style="border: 1px solid black"> <code>\\x7B</code>–<code>\\x7E</code>
 *</td>
 *<td style="border: 1px solid black"> <code>\172</code>–<code>\176</code>
 *</td>
 *<td style="border: 1px solid black"> <code>{|}~</code>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td></tr>
 *<tr>
 *<td style="border: 1px solid black"> 127
 *</td>
 *<td style="border: 1px solid black"> <code>\\x7F</code>
 *</td>
 *<td style="border: 1px solid black"> <code>\177</code>
 *</td>
 *<td style="border: 1px solid black"> backspace character (<code>DEL</code>)
 *</td>
 *<td style="background: #90ff90; color: black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-yes"><b><code>≠0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td>
 *<td style="background:#ff9090; color:black; vertical-align: middle; text-align: center; border: 1px solid black;" class="table-no"> <b><code>0</code></b>
 *</td></tr>
 *</table>
 * <br>
 */

} // namespace seqan3
