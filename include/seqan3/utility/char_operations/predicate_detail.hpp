// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides parse conditions for tokenization.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <array>
#include <cctype>
#include <concepts>
#include <stdexcept>
#include <string>

#include <seqan3/utility/concept.hpp>
#include <seqan3/utility/detail/type_name_as_string.hpp>
#include <seqan3/utility/type_traits/basic.hpp>

namespace seqan3::detail
{

// ----------------------------------------------------------------------------
// constexpr_pseudo_bitset (could be replaced with constexpr std::bitset or github.com/ClaasBontus/bitset2)
// ----------------------------------------------------------------------------

/*!\brief A data structure that implements a subset of std::bitset as constexpr.
 * \ingroup utility_char_operations
 * \tparam N                The number of bits.
 */
template <size_t N>
class constexpr_pseudo_bitset : public std::array<bool, N>
{
private:
    //!\brief The base type.
    using base_t = std::array<bool, N>;

public:
    //!\brief Inherit constructors.
    using base_t::base_t;

    //!\brief Return a new bitset that is a logical disjunction of the two given ones.
    constexpr constexpr_pseudo_bitset operator|(constexpr_pseudo_bitset rhs) const noexcept
    {
        for (size_t i = 0; i < N; ++i)
            rhs[i] = rhs[i] || base_t::operator[](i);

        return rhs;
    }
    //!\brief Return a new bitset with all bits flipped.
    constexpr constexpr_pseudo_bitset operator~() const noexcept
    {
        constexpr_pseudo_bitset ret{};
        for (size_t i = 0; i < N; ++i)
            ret[i] = !base_t::operator[](i);

        return ret;
    }
};

// ----------------------------------------------------------------------------
// condition_message_v
// ----------------------------------------------------------------------------

/*!\brief Defines a compound std::string consisting of all given conditions separated by the operator-name `op`.
 * \ingroup utility_char_operations
 *
 * \tparam op               non-type template parameter specifying the separator character, e.g. '|'.
 * \tparam condition_head_t The first condition type in the message. Ensures that there is at least one type.
 * \tparam condition_ts     Remaining list of conditions separated by `op`.
 * \relates seqan3::detail::char_predicate
 */
template <char op, typename condition_head_t, typename... condition_ts>
inline std::string const condition_message_v{
    std::string{"("}
    + (condition_head_t::msg + ... + (std::string{" "} + std::string{{op, op}} + std::string{" "} + condition_ts::msg))
    + std::string{")"}};

// ----------------------------------------------------------------------------
// char_predicate
// ----------------------------------------------------------------------------

//!\cond
template <typename condition_t>
struct char_predicate_base;
//!\endcond

/*!\interface seqan3::detail::char_predicate <>
 * \brief An internal concept to check if an object fulfills the requirements of a seqan3::detail::char_predicate.
 * \ingroup utility_char_operations
 *
 * \details
 *
 * An object of the type must be invocable with a std::integral type and supply a static constexpr `msg` member of type
 * std::string.
 */
//!\cond
template <typename condition_t>
concept char_predicate = requires {
    requires std::predicate<std::remove_reference_t<condition_t>, char>;
    requires std::is_base_of_v<char_predicate_base<std::remove_cvref_t<condition_t>>, std::remove_cvref_t<condition_t>>;

    std::remove_reference_t<condition_t>::msg;

    //The msg type can be added with a std::string.
    {
        std::string{} + std::remove_reference_t<condition_t>::msg
    } -> std::convertible_to<decltype(std::remove_reference_t<condition_t>::msg)>;
};
//!\endcond

/*!\name Requirements for seqan3::detail::char_predicate
 * \brief You can expect the variable and the predicate function on all types that satisfy seqan3::output_stream_over.
 * \{
 */
/*!\fn      bool operator()(char_type c);
 * \brief   predicate function to test if `c` satisfies the given condition.
 * \memberof seqan3::detail::char_predicate
 * \param   c The character to be tested.
 * \returns `true` on success, `false` otherwise.
 *
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */

/*!\var static constexpr auto msg
 * \memberof seqan3::detail::char_predicate
 * \brief Defines the condition msg. The type is deduced from the constant expression in the definition of the variable.
 */
//!\}

// ----------------------------------------------------------------------------
// char_predicate
// ----------------------------------------------------------------------------

//!\cond
template <char_predicate... condition_ts>
    requires (sizeof...(condition_ts) >= 2)
struct char_predicate_disjunction;

template <char_predicate condition_t>
struct char_predicate_negator;
//!\endcond

/*!\brief An abstract [CRTP](https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern) base class for
 *        parse conditions to add logical disjunction and negation operator.
 * \ingroup utility_char_operations
 * \tparam derived_t The parse condition type to be extended with the logical operators.
 *                   Must model seqan3::detail::char_predicate.
 */
template <typename derived_t>
struct char_predicate_base
{
    //!\brief Type of the data field; as soon as std::bitset is fully constexpr, use that instead!
    using data_t = constexpr_pseudo_bitset<257>; // sizeof(char) plus EOF

    /*!\name Logical operators
     * \brief Adds logical operators to allow logical disjunction, conjunction and negation on parse conditions.
     * \{
     */
    //!\brief Combines the result of two seqan3::detail::char_predicate via logical disjunction.
    template <char_predicate rhs_t>
    constexpr auto operator||(rhs_t const &) const
    {
        return char_predicate_disjunction<derived_t, rhs_t>{};
    }

    //!\brief Return a new condition with all bits flipped.
    constexpr auto operator!() const
    {
        return char_predicate_negator<derived_t>{};
    }
    //!\}

    /*!\name Function call operator
     * \{
     */
    //!\brief Invokes the condition on `val`.
    template <std::integral value_t>
    constexpr bool operator()(value_t const val) const noexcept
        requires (sizeof(value_t) == 1)
    {
        return derived_t::data[static_cast<unsigned char>(val)];
    }

    //!\overload
    template <std::integral value_t>
    constexpr bool operator()(value_t const val) const noexcept
        requires (sizeof(value_t) != 1)
    {
        // std::char_traits is only guaranteed to be defined for character types.
        // libc++ deprecates other specialisations in llvm-17, and removes them in llvm-18.
        // We map the non-character types to corresponding chracter types.
        // For example, `seqan3::is_eof(EOF)` will call this function with `value_t == int`.
        using char_value_t = std::conditional_t<
            seqan3::builtin_character<value_t>,
            value_t,
            std::conditional_t<
                std::same_as<value_t, std::char_traits<char>::int_type>,
                char,
                std::conditional_t<
                    std::same_as<value_t, std::char_traits<wchar_t>::int_type>,
                    wchar_t,
                    std::conditional_t<
                        std::same_as<value_t, std::char_traits<char8_t>::int_type>,
                        char8_t,
                        std::conditional_t<
                            std::same_as<value_t, std::char_traits<char16_t>::int_type>,
                            char16_t,
                            std::conditional_t<std::same_as<value_t, std::char_traits<char32_t>::int_type>,
                                               char32_t,
                                               void>>>>>>;
        static_assert(!std::same_as<char_value_t, void>, "There is no valid character representation.");
        using char_trait = std::char_traits<char_value_t>;
        return (static_cast<std::make_unsigned_t<value_t>>(val) < 256) ? operator()(static_cast<uint8_t>(val))
             : (char_trait::eq_int_type(val, char_trait::eof()))       ? derived_t::data[256]
                                                                       : false;
    }
    //!\}

    /*!\name Output functions
     * \{
     */
    //!\brief Returns the message representing this condition as std::string.
    std::string message() const
    {
        return derived_t::msg;
    }
    //!\}
};

// ----------------------------------------------------------------------------
// char_predicate_disjunction
// ----------------------------------------------------------------------------

/*!\brief Logical disjunction operator for parse conditions.
 * \implements seqan3::detail::char_predicate
 * \tparam condition_ts Template parameter pack over all parse condition types. Must contain at least 2 template parameters.
 *                      Must model seqan3::detail::char_predicate.
 * \ingroup utility_char_operations
 */
template <char_predicate... condition_ts>
    requires (sizeof...(condition_ts) >= 2)
struct char_predicate_disjunction : public char_predicate_base<char_predicate_disjunction<condition_ts...>>
{
    //!\brief The message representing the disjunction of the associated conditions.
    static inline std::string const msg = detail::condition_message_v<'|', condition_ts...>;

    //!\brief The base type.
    using base_t = char_predicate_base<char_predicate_disjunction<condition_ts...>>;

    //!\brief Import the data type from the base class.
    using typename base_t::data_t;
    //!\brief The look-up table that is used to evaluate the input.
    static constexpr data_t data = (condition_ts::data | ...);
};

/*!\brief Logical not operator for a parse condition.
 * \implements seqan3::detail::char_predicate
 * \tparam condition_t Template parameter to apply the not-operator for.
 *                     Must model seqan3::detail::char_predicate.
 * \ingroup utility_char_operations
 */
template <char_predicate condition_t>
struct char_predicate_negator : public char_predicate_base<char_predicate_negator<condition_t>>
{
    //!\brief The message representing the negation of the associated condition.
    static inline std::string const msg = std::string{'!'} + condition_t::msg;

    //!\brief The base type.
    using base_t = char_predicate_base<char_predicate_negator<condition_t>>;

    //!\brief Import the data type from the base class.
    using typename base_t::data_t;
    //!\brief The look-up table that is used to evaluate the input.
    static constexpr data_t data = ~condition_t::data;
};

// ----------------------------------------------------------------------------
// is_in_interval_type
// ----------------------------------------------------------------------------

/*!\brief Parse condition that checks if a given value is in the range of `rng_beg` and `interval_last`.
 * \ingroup utility_char_operations
 * \implements seqan3::detail::char_predicate
 * \tparam interval_first non-type template parameter denoting the begin of the allowed range.
 *                        Must be less than or equal to `interval_last`.
 * \tparam interval_last non-type template parameter denoting the end of the allowed range.
 *                       Must be greater than or equal to `interval_first`.
 */
template <uint8_t interval_first, uint8_t interval_last>
    requires (interval_first <= interval_last)
struct is_in_interval_type : public char_predicate_base<is_in_interval_type<interval_first, interval_last>>
{
    //!\brief The message representing this condition.
    static inline std::string const msg = std::string{"is_in_interval<'"} + std::string{interval_first}
                                        + std::string{"', '"} + std::string{interval_last} + std::string{"'>"};

    //!\brief The base type.
    using base_t = char_predicate_base<is_in_interval_type<interval_first, interval_last>>;

    //!\brief Import the data type from the base class.
    using typename base_t::data_t;
    //!\brief The look-up table that is used to evaluate the input.
    static constexpr data_t data = []() constexpr
    {
        data_t ret{};

        for (size_t i = interval_first; i <= static_cast<size_t>(interval_last); ++i)
            ret[i] = true;

        return ret;
    }();
};

// ----------------------------------------------------------------------------
// is_char_type
// ----------------------------------------------------------------------------

/*!\brief Parse condition that checks if a given value is equal to `char_v`.
 * \ingroup utility_char_operations
 * \implements seqan3::detail::char_predicate
 * \tparam char_v non-type template parameter with the value that should be checked against.
 */
template <int char_v>
struct is_char_type : public char_predicate_base<is_char_type<char_v>>
{
    static_assert(char_v == EOF || static_cast<uint64_t>(char_v) < 256, "TODO");

    //!\brief The message representing this condition.
    static inline std::string const msg =
        std::string{"is_char<'"} + ((char_v == EOF) ? std::string{"EOF"} : std::string{char_v}) + std::string{"'>"};

    //!\brief The base type.
    using base_t = char_predicate_base<is_char_type<char_v>>;

    //!\brief Import the data type from the base class.
    using typename base_t::data_t;
    //!\brief The look-up table that is used to evaluate the input.
    static constexpr data_t data = []() constexpr
    {
        data_t ret{};

        if (char_v == EOF)
            ret[256] = true;
        else
            ret[static_cast<uint8_t>(char_v)] = true;

        return ret;
    }();
};

} // namespace seqan3::detail
