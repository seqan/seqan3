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
 * \brief Provides parse conditions for tokenization.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <cctype>
#include <stdexcept>
#include <string>

#include <seqan3/alphabet/all.hpp>
#include <seqan3/core/detail/reflection.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/range/container/constexpr_string.hpp>
#include <seqan3/std/concept/core_language.hpp>
#include <seqan3/std/concept/callable.hpp>

namespace seqan3::detail
{

// ----------------------------------------------------------------------------
// condition_message_v
// ----------------------------------------------------------------------------

/*!\brief Defines a compound seqan3::constexpr_string consisting of all given conditions separated by the
 *        operator-name `op`.
 * \ingroup stream
 * \tparam op               non-type template parameter specifying the separator character, e.g. '|'.
 * \tparam condition_head_t The first condition type in the message. Ensures that there is at least one type.
 * \tparam condition_ts     Remaining list of conditions separated by `op`.
 * \relates seqan3::detail::parse_condition
 */
template <char op, typename condition_head_t, typename ...condition_ts>
constexpr constexpr_string condition_message_v
{
    constexpr_string{"("} +
    (condition_head_t::msg + ... +
        (constexpr_string{" "} + constexpr_string{{op, op, '\0'}} + constexpr_string{" "} + condition_ts::msg)) +
    constexpr_string{")"}
};

// ----------------------------------------------------------------------------
// parse_condition_concept
// ----------------------------------------------------------------------------

//!\cond
template <typename condition_t>
class parse_condition;
//!\endcond

/*!\interface seqan3::detail::parse_condition_concept <>
 * \brief An internal concept to check if an object fulfills the requirements of a seqan3::detail::parse_condition.
 * \ingroup stream
 *
 * \details
 *
 * The must be invocable with an seqan3::char_adaptation_concept type and supply a static constexpr `msg` member of type
 * seqan3::constexpr_string.
 */
//!\cond
template <typename condition_t>
concept bool parse_condition_concept = requires
{
    requires predicate_concept<std::remove_reference_t<condition_t>, char>;
    requires std::is_base_of_v<parse_condition<remove_cvref_t<condition_t>>,
                               remove_cvref_t<condition_t>>;

    std::remove_reference_t<condition_t>::msg;

    //The msg type can be added with a constexpr_string.
    { constexpr_string<0>{} + std::remove_reference_t<condition_t>::msg } ->
        decltype(std::remove_reference_t<condition_t>::msg);
};
//!\endcond

/*!\name Requirements for seqan3::detail::parse_condition_concept
 * \brief You can expect the variable and the predicate function on all types that satisfy seqan3::ostream_concept.
 * \{
 */
/*!\fn      bool operator()(char_type c);
 * \brief   Predicate function to test if `c` satisfies the given condition.
 * \memberof seqan3::detail::parse_condition_concept
 * \param   c The character to be tested.
 * \returns `true` on success, `false` otherwise.
 *
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */

/*!\var static constexpr auto msg
 * \memberof seqan3::detail::parse_condition_concept
 * \brief Defines the condition msg. The type is deduced from the constant expression in the definition of the variable.
 */
//!\}

// ----------------------------------------------------------------------------
// make_printable
// ----------------------------------------------------------------------------

/*!\brief Returns a printable value for the given character `c`.
 * \param[in] c The character to be represented as printable string.
 * \return    a std::string containing a printable version of the given character `c`.
 *
 * \details
 *
 * Some characters, e.g. control commands cannot be printed. This function converts them to a std::string
 * containing the visual representation of this character. For all control commands the value `'CTRL'` is returned.
 *
 * ### Exception
 *
 * Strong exception guarantee is given.
 *
 * ### Complexity
 *
 * Constant.
 *
 * ### Concurrency
 *
 * Thread-safe.
 */
template <typename char_type>
inline std::string
make_printable(char_type const c)
{
    switch (c)
    {
        case '\0':                   return "'\\0'";
        case '\t':                   return "'\\t'";
        case '\n':                   return "'\\n'";
        case '\v':                   return "'\\v'";
        case '\f':                   return "'\\f'";
        case '\r':                   return "'\\r'";
        case static_cast<char>(127): return "'DEL'";
        default:
        {
            if ((c >= static_cast<char>(1) && c <= static_cast<char>(8)) ||
                (c >= static_cast<char>(14) && c <= static_cast<char>(31)))
                return "'CTRL'";
            else
                return {'\'', c, '\''};
        }
    }
}

// ----------------------------------------------------------------------------
// parse_condition
// ----------------------------------------------------------------------------

//!\cond
template <parse_condition_concept... condition_ts>
    requires sizeof...(condition_ts) >= 2
struct parse_condition_combiner;

template <parse_condition_concept condition_t>
struct parse_condition_negator;
//!\endcond

/*!\brief An abstract [CRTP](https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern) base class for
 *        parse conditions to add logical disjunction and negation operator.
 * \ingroup stream
 * \tparam derived_t The parse condition type to be extended with the logical operators.
 *                   Must satisfy the seqan3::detail::parse_condition_concept.
 */
template <typename derived_t>
class parse_condition
{
private:

    //!\brief Friend declaration for the derived type to access the private constructors.
    friend derived_t;

    /*!\name Constructor, destructor and assignment
     * \brief This base class is abstract and all constructors are declared private.
     * \{
     */
    constexpr parse_condition() = default;
    constexpr parse_condition(parse_condition const &) = default;
    constexpr parse_condition(parse_condition &&) = default;
    constexpr parse_condition & operator=(parse_condition const &) = default;
    constexpr parse_condition & operator=(parse_condition &&) = default;
    ~parse_condition() =  default;
    //!\}

public:

    /*!\name Logical operators
     * \brief Adds logical operators to allow logical disjunction, conjunction and negation on parse conditions.
     * \{
     */
    //!\brief Combines the result of two seqan3::detail::parse_condition via logical disjunction.
    template <parse_condition_concept rhs_derived_t>
    constexpr auto operator||(parse_condition<rhs_derived_t> const &) const
    {
        return parse_condition_combiner<derived_t, rhs_derived_t>{};
    }

    //!\brief Negates the result of a seqan3::detail::parse_condition
    constexpr auto operator!() const
    {
        return parse_condition_negator<derived_t>{};
    }
    //!\}

    /*!\name Function call operator
     * \{
     */
    //!\brief Invokes the condition on `val`.
    template <char_adaptation_concept value_t>
    constexpr bool
    operator()(value_t && val) const noexcept(std::is_nothrow_invocable_r_v<bool, derived_t, value_t>)
    {
        return std::invoke(derived_t(), std::forward<value_t>(val));
    }
    //!\}

    /*!\name Output functions
     * \{
     */
    //!\brief Returns the message representing this condition as std::string.
    std::string message() const
        requires parse_condition_concept<derived_t>
    {
        return derived_t::msg.string();
    }
    //!\}
};

// ----------------------------------------------------------------------------
// parse_condition_combiner
// ----------------------------------------------------------------------------

/*!\brief Logical disjunction operator for parse conditions.
 * \implements seqan3::detail::parse_condition_concept
 * \tparam condition_ts Template parameter pack over all parse condition types. Must contain at least 2 template parameters.
 *                      Must satisfy the seqan3::detail::parse_condition_concept.
 * \ingroup stream
 */
template <parse_condition_concept... condition_ts>
//\cond
    requires sizeof...(condition_ts) >= 2
//\endcond
struct parse_condition_combiner : public parse_condition<parse_condition_combiner<condition_ts...>>
{
    //!\brief The message representing the disjunction of the associated conditions.
    static constexpr auto msg = detail::condition_message_v<'|', condition_ts...>;

    //!\brief Inherit constructors from CRTP base class.
    using parse_condition<parse_condition_combiner<condition_ts...>>::parse_condition;

    /*!\brief Invokes the condition check for `val` and returns the logical disjunction of their results.
     * \tparam    char_t The type of the character to test. Must satisfy seqan3::char_adaptation_concept.
     * \param[in] val    The value to test.
     * \return `bool`; `false` if the condition is met, `true` otherwise.
     *
     * \details
     *
     * ### Complexity
     *
     * Linear in the number of conditions.
     *
     * ### Exception
     *
     * No-throw guarantee if the expression `(... && std::is_nothrow_invocable_r_v<bool, condition_ts, char_t &&>)`
     * evaluates to `true`.
     *
     * ### Concurrency
     *
     * Thread-safe.
     */
    template <char_adaptation_concept char_t>
    constexpr bool
    operator()(char_t && val) const noexcept((... && std::is_nothrow_invocable_r_v<bool, condition_ts, char_t &&>))
    {
        return (... || std::invoke(condition_ts(), std::forward<char_t>(val)));
    }
};

/*!\brief Logical not operator for a parse condition.
 * \implements seqan3::detail::parse_condition_concept
 * \tparam condition_t Template parameter to apply the not-operator for.
 *                     Must satisfy the seqan3::detail::parse_condition_concept.
 * \ingroup stream
 */
template <parse_condition_concept condition_t>
struct parse_condition_negator : public parse_condition<parse_condition_negator<condition_t>>
{
    //!\brief The message representing the negation of the associated condition.
    static constexpr auto msg = constexpr_string{'!'} + condition_t::msg;

    //!\brief Inherit constructors from CRTP base class.
    using parse_condition<parse_condition_negator<condition_t>>::parse_condition;

    /*!\brief Invokes the condition check for `val` and returns the negated result.
     * \tparam    char_t The type of the character to test. Must satisfy seqan3::char_adaptation_concept.
     * \param[in] val    The value to test.
     * \return `bool`; `false` if the condition is met, `true` otherwise.
     *
     * \details
     *
     * ### Complexity
     *
     * Constant
     *
     * ### Exception
     *
     * No-throw guarantee if the expression `std::is_nothrow_invocable_r_v<bool, condition_t, value_t &&>`
     * evaluates to `true`.
     *
     * ### Concurrency
     *
     * Thread-safe.
     */
    template <char_adaptation_concept value_t>
    constexpr bool
    operator()(value_t && val) const noexcept(std::is_nothrow_invocable_r_v<bool, condition_t, value_t &&>)
    {
        return !std::invoke(condition_t(), std::forward<value_t>(val));
    }
};

} // namespace seqan3::detail

namespace seqan3
{

// ----------------------------------------------------------------------------
// is_in_interval
// ----------------------------------------------------------------------------

/*!\brief Parse condition that checks if a given value is in the range of `rng_beg` and `interval_last`.
 * \ingroup stream
 * \implements seqan3::detail::parse_condition_concept
 * \tparam interval_first non-type template parameter denoting the begin of the allowed range.
 *                        Must be less than or equal to `interval_last`.
 * \tparam interval_last non-type template parameter denoting the end of the allowed range.
 *                       Must be greater than or equal to `interval_first`.
 */
template <char interval_first, char interval_last>
//\cond
    requires interval_first <= interval_last
//\endcond
struct is_in_interval : public detail::parse_condition<is_in_interval<interval_first, interval_last>>
{
    //!\brief The message representing this condition.
    static constexpr constexpr_string msg = constexpr_string{"is_in_interval<'"} +
                                            constexpr_string{interval_first}         +
                                            constexpr_string{"', '"}          +
                                            constexpr_string{interval_last}         +
                                            constexpr_string{"'>"};

    //!\brief Inherit constructors from CRTP base class.
    using detail::parse_condition<is_in_interval<interval_first, interval_last>>::parse_condition;

    /*!\brief Invokes the condition check for `val`.
     * \tparam    char_t The type of the character to test. Must satisfy seqan3::char_adaptation_concept.
     * \param[in] val    The value to test.
     * \return `bool`; `true` if the condition is met, `false` otherwise.
     *
     * \details
     *
     * ### Complexity
     *
     * Constant
     *
     * ### Exception
     *
     * No-throw guarantee.
     *
     * ### Concurrency
     *
     * Thread-safe.
     */
    template <char_adaptation_concept char_t>
    constexpr bool
    operator()(char_t const & val) const noexcept
    {
        return (static_cast<uint64_t>(interval_first) <= static_cast<uint64_t>(val)) &&
               (static_cast<uint64_t>(val) <= static_cast<uint64_t>(interval_last));
    }
};

// ----------------------------------------------------------------------------
// is_in_alphabet
// ----------------------------------------------------------------------------

/*!\brief Parse condition that checks if a given value is within the given alphabet `alphabet_t`.
 * \ingroup stream
 * \implements seqan3::detail::parse_condition_concept
 * \tparam alphabet_t The alphabet type. Must satisfy the seqan3::alphabet_concept.
 */
template <alphabet_concept alphabet_t>
struct is_in_alphabet : public detail::parse_condition<is_in_alphabet<alphabet_t>>
{
    //!\brief The message representing this condition.
    static constexpr auto msg = constexpr_string{"is_in_alphabet<"} +
                                constexpr_string{detail::get_display_name_v<alphabet_t>} +
                                constexpr_string{">"};

    //!\brief Inherit constructors from CRTP base class.
    using detail::parse_condition<is_in_alphabet<alphabet_t>>::parse_condition;

    /*!\brief Invokes the condition check for `val`.
     * \tparam    char_t The type of the character to test. Must satisfy seqan3::char_adaptation_concept.
     * \param[in] val    The value to test.
     * \return `bool`; `true` if the condition is met, `false` otherwise.
     *
     * \details
     *
     * ### Complexity
     *
     * Constant
     *
     * ### Exception
     *
     * No-throw guarantee.
     *
     * ### Concurrency
     *
     * Thread-safe.
     */
    template <char_adaptation_concept char_t>
    constexpr bool
    operator()(char_t && val) const noexcept(std::is_same_v<char_t, char>)
    {
        if constexpr (!std::is_same_v<char_t, char>)
        {  // Check if alphabet is able to represent the underlying char type of the tested alphabet.
            if (static_cast<uint64_t>(val) >
                static_cast<uint64_t>(std::numeric_limits<underlying_char_t<alphabet_t>>::max()))
                return false;
        }

        return to_char(assign_char(alphabet_t{}, val)) == std::toupper(static_cast<uint8_t>(val));
    }
};

// ----------------------------------------------------------------------------
// is_char
// ----------------------------------------------------------------------------

/*!\brief Parse condition that checks if a given value is equal to `char_v`.
 * \ingroup stream
 * \implements seqan3::detail::parse_condition_concept
 * \tparam alphabet_t non-type template parameter with the value that should be checked against.
 */
template <char char_v>
struct is_char : public detail::parse_condition<is_char<char_v>>
{
    //!\brief The message representing this condition.
    static constexpr auto msg = constexpr_string{"is_char<'"} +
                                constexpr_string{char_v}     +
                                constexpr_string("'>");

    //!\brief Inherit constructors from CRTP base class.
    using detail::parse_condition<is_char<char_v>>::parse_condition;

    /*!\brief Invokes the condition check for `val`.
     * \tparam    char_t The type of the character to test. Must satisfy seqan3::char_adaptation_concept.
     * \param[in] val    The value to test.
     * \return `bool`; `true` if the condition is met, `false` otherwise.
     *
     * \details
     *
     * ### Complexity
     *
     * Constant
     *
     * ### Exception
     *
     * No-throw guarantee.
     *
     * ### Concurrency
     *
     * Thread-safe.
     */
    template <char_adaptation_concept char_t>
    constexpr bool
    operator()(char_t && val) const noexcept
    {
        return val == char_v;
    }
};

// ----------------------------------------------------------------------------
// General Purpose Parse Conditions
// ----------------------------------------------------------------------------

/*!\name Parse conditions
 * \brief A set of function objects to check if a character from an input source fulfills certain characteristics.
 * \ingroup stream
 * \{
 */

/*!\brief Checks wether `c` is a control character.
 * \ingroup stream
 *
 * \details
 *
 * This function like object can be used to check if a character `c` is a control character.
 * For the standard ASCII character set, control characters are those between ASCII codes 0x00 (NUL) and 0x1f (US) and
 * 0x7f (DEL).
 *
 * ### Example
 *
 * ```cpp
 * seqan3::is_cntrl('\0');  // returns true.
 * ```
 */
constexpr detail::parse_condition_combiner<is_in_interval<'\0', static_cast<char>(31)>,
                                           is_char<static_cast<char>(127)>> is_cntrl;

/*!\brief Checks wether `c` is a printable character.
 * \ingroup stream
 *
 * \details
 *
 * This function like object can be used to check if a character `c` is a printable character.
 * For the standard ASCII character set, printable characters are those between ASCII codes 0x20 (space) and 0x7E (`~`).
 *
 * ### Example
 *
 * ```cpp
 * seqan3::is_print(' ');  // returns true.
 * ```
 */
constexpr is_in_interval<' ', '~'> is_print;

/*!\brief Checks wether `c` is a space character.
 * \ingroup stream
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
 * ```cpp
 * seqan3::is_space('\n');  // returns true.
 * ```
 */
constexpr detail::parse_condition_combiner<is_in_interval<'\t', '\r'>, is_char<' '>> is_space;

/*!\brief Checks wether `c` is a blank character.
 * \ingroup stream
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
 * ```cpp
 * seqan3::is_blank('\t');  // returns true.
 * ```
 */
constexpr detail::parse_condition_combiner<is_char<'\t'>, is_char<' '>> is_blank;

/*!\brief Checks wether `c` is a graphic character.
 * \ingroup stream
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
 * * punctuation characters (!"#$%&'()*+,-./:;<=>?\@[\]^_`{|}~)
 *
 * ### Example
 *
 * ```cpp
 * seqan3::is_graph('%');  // returns true.
 * ```
 */
constexpr is_in_interval<'!', '~'> is_graph;

/*!\brief Checks wether `c` is a punctuation character.
 * \ingroup stream
 *
 * \details
 *
 * This function like object can be used to check if a character `c` is a punctuation character.
 * For the standard ASCII character set, the following characters are punctuation characters:
 *
 * * punctuation characters (!"#$%&'()*+,-./:;<=>?\@[\]^_`{|}~)
 *
 * ### Example
 *
 * ```cpp
 * seqan3::is_punct(':');  // returns true.
 * ```
 */
constexpr detail::parse_condition_combiner<is_in_interval<'!', '/'>, is_in_interval<':', '@'>,
                                           is_in_interval<'[', '`'>, is_in_interval<'{', '~'>> is_punct;

/*!\brief Checks wether `c` is a alphanumeric character.
 * \ingroup stream
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
 * ```cpp
 * seqan3::is_alnum('9');  // returns true.
 * ```
 */
constexpr detail::parse_condition_combiner<is_in_interval<'0','9'>, is_in_interval<'A','Z'>, is_in_interval<'a','z'>> is_alnum;

/*!\brief Checks wether `c` is a alphabetical character.
 * \ingroup stream
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
 * ```cpp
 * seqan3::is_alpha('z');  // returns true.
 * ```
 */
constexpr detail::parse_condition_combiner<is_in_interval<'A', 'Z'>, is_in_interval<'a', 'z'>> is_alpha;

/*!\brief Checks wether `c` is a upper case character.
 * \ingroup stream
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
 * ```cpp
 * seqan3::is_upper('K');  // returns true.
 * ```
 */
constexpr is_in_interval<'A', 'Z'> is_upper;

/*!\brief Checks wether `c` is a lower case character.
 * \ingroup stream
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
 * ```cpp
 * seqan3::is_lower('a');  // returns true.
 * ```
 */
constexpr is_in_interval<'a', 'z'> is_lower;

/*!\brief Checks wether `c` is a digital character.
 * \ingroup stream
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
 * ```cpp
 * seqan3::is_digit('1');  // returns true.
 * ```
 */
constexpr is_in_interval<'0', '9'> is_digit;

/*!\brief Checks wether `c` is a hexadecimal character.
 * \ingroup stream
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
 * ```cpp
 * seqan3::is_xdigit('e');  // returns true.
 * ```
 */
constexpr detail::parse_condition_combiner<is_in_interval<'0', '9'>, is_in_interval<'A', 'F'>,
                                           is_in_interval<'a', 'f'>> is_xdigit;
//!\}

/*!\brief A condition checker, that wraps a parse condition and throws a specified exception if the condition was not
 *        met.
 * \ingroup stream
 * \tparam condition_type The wrapped parse condition type to be use for testing
 *
 * \details
 *
 * The `parse_asserter` type is used to enforce that a parsed character satisfies certain conditions.
 * For example, when reading an input fasta file and the expected alphabet is a dna4
 * but the actual data contained in the file is based on amino acids.
 * Thus, the condition would not be satisfied, causing the exception to be thrown.
 *
 * ```cpp
 * using namespace seqan3;
 *
 * std::istringstream istr{"ATZE"};
 *
 * std::istream_iterator<char> it{istr},
 * parse_asserter<is_in_alphabet<dna5>> asserter{};
 *
 * while (it != std::istream_iterator<char>{})
 * {
 *     asserter(*it);  // will throw when reading `Z` from the input stream.
 *     ++it;
 * }
 * ```
 *
 * ### Deduction Guide
 *
 * The `parse_asserter` class itself is stateless. Still, it holds a member of the given parse_condition, in order
 * to allow template deduction from a passed argument.
 * The following listing shows the alternative definition, which generates the same `parse_asserter` type as above.
 *
 * ```cpp
 * parse_asserter asserter{is_alnum};
 * ```
 */
template <typename condition_type>
struct parse_asserter
{
    //!\brief Stores an instance of the stateless condition.
    condition_type cond;

    /*!\brief Checks if the given character satisfies the associated parse condition.
     * \param[in] c The character to be checked. Must satisfy the seqan3::char_adaptation_concept.
     *
     * \details
     *
     * ### Complexity
     *
     * Depends on the underlying condition to be checked. In the regular case it is constant.
     *
     * ### Exception
     *
     * Throws seqan3::parse_error if the associated condition is not met.
     *
     * ### Concurrency
     *
     * Thread-safe.
     */
    template <char_adaptation_concept char_type>
    void operator()(char_type && c) const
    {
        if (!std::invoke(cond, std::forward<char_type>(c)))
        {
            using namespace std::literals;
            // I can not assure that all parameters are convertible to c.
            throw parse_error{"Parsed value <"s + detail::make_printable(to_char(c)) +
                              "> which does not fulfill the following condition: "s  +
                              cond.message()};
        }
    }
};

/*!\name Deduction Guide
 * \brief Deduction guide for parse_asserter.
 * \ingroup stream
 * \relates parse_asserter
 * \{
 */
//!\brief Deduction guide to infer the template argument for the condition type from the constructor argument.
template <detail::parse_condition_concept parse_cond_type>
parse_asserter(parse_cond_type) -> parse_asserter<parse_cond_type>;
//!\}

/*!\name Parse conditions
 *
 * \details
 *
 * Parse conditions are function like objects that can be used to check if a character `c` fulfills certain
 * constraints. There are three basic condition types:
 *
 * * seqan3::is_in_alphabet: Checks if the given character is part of the specified alphabet.
 * * seqan3::is_in_interval: Checks if the given character is within specified range of ASCII characters.
 * * seqan3::is_char: Checks if the character is equal to the specified ASCII character.
 *
 * These checks are necessary when parsing input streams, where the input characters need to be checked for consistency,
 * in order to detect ill-formed input data that could cause otherwise lead to unexpected behavior and run-time errors
 * that difficult to track down.
 *
 * ### Disjunction and Negation
 *
 * All functors can be combined with the `||-operator` or negated via the `!-operator`, such that chains of conditions
 * can be easily created.
 *
 * ```cpp
 * auto my_cond = is_char<'#'>{} || is_char<'%'>{};
 * bool is_comment = my_cond(*stream_it);
 * ```
 *
 * ### General Purpose Parse Conditions
 *
 * There are 12 predefined parse conditions that are all compositions of the previously mentioned parse conditions.
 * The following table lists the predefined parse conditions and which constraints are associated with them.
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
 *<td style="border: 1px solid black"> <code>!"#$%&amp;'()*+,-./</code>
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
