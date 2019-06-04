// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides parse conditions for tokenization.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 * \author Hannes Hauswedell <hannes.hauswedell AT fu-berlin.de>
 */

#pragma once

#include <array>
#include <cctype>
#include <stdexcept>
#include <string>

#include <seqan3/std/concepts>

#include <seqan3/alphabet/concept.hpp>
#include <seqan3/core/detail/reflection.hpp>
#include <seqan3/core/type_traits/basic.hpp>
#include <seqan3/range/container/small_string.hpp>

namespace seqan3::detail
{

// ----------------------------------------------------------------------------
// constexpr_pseudo_bitset (could be replaced with constexpr std::bitset or github.com/ClaasBontus/bitset2)
// ----------------------------------------------------------------------------

/*!\brief A data structure that implements a subset of std::bitset as constexpr.
 * \ingroup stream
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

/*!\brief Defines a compound seqan3::small_string consisting of all given conditions separated by the
 *        operator-name `op`.
 * \ingroup stream
 * \tparam op               non-type template parameter specifying the separator character, e.g. '|'.
 * \tparam condition_head_t The first condition type in the message. Ensures that there is at least one type.
 * \tparam condition_ts     Remaining list of conditions separated by `op`.
 * \relates seqan3::detail::char_predicate
 */
template <char op, typename condition_head_t, typename ...condition_ts>
inline small_string constexpr condition_message_v
{
    small_string{"("} +
    (condition_head_t::msg + ... +
        (small_string{" "} + small_string{{op, op, '\0'}} + small_string{" "} + condition_ts::msg)) +
    small_string{")"}
};

// ----------------------------------------------------------------------------
// CharPredicate
// ----------------------------------------------------------------------------

//!\cond
template <typename condition_t>
class char_predicate_base;
//!\endcond

/*!\interface seqan3::detail::CharPredicate <>
 * \brief An internal concept to check if an object fulfills the requirements of a seqan3::detail::char_predicate.
 * \ingroup stream
 *
 * \details
 *
 * An object of the type must be invocable with a std::Integral type and supply a static constexpr `msg` member of type
 * seqan3::small_string.
 */
//!\cond
template <typename condition_t>
SEQAN3_CONCEPT CharPredicate = requires
{
    requires std::Predicate<std::remove_reference_t<condition_t>, char>;
    requires std::is_base_of_v<char_predicate_base<remove_cvref_t<condition_t>>,
                               remove_cvref_t<condition_t>>;

    std::remove_reference_t<condition_t>::msg;

    //The msg type can be added with a small_string.
    { small_string<0>{} + std::remove_reference_t<condition_t>::msg } ->
        decltype(std::remove_reference_t<condition_t>::msg);
};
//!\endcond

/*!\name Requirements for seqan3::detail::CharPredicate
 * \brief You can expect the variable and the predicate function on all types that satisfy seqan3::OStream.
 * \{
 */
/*!\fn      bool operator()(char_type c);
 * \brief   Predicate function to test if `c` satisfies the given condition.
 * \memberof seqan3::detail::CharPredicate
 * \param   c The character to be tested.
 * \returns `true` on success, `false` otherwise.
 *
 * \attention This is a concept requirement, not an actual function (however types satisfying this concept
 * will provide an implementation).
 */

/*!\var static constexpr auto msg
 * \memberof seqan3::detail::CharPredicate
 * \brief Defines the condition msg. The type is deduced from the constant expression in the definition of the variable.
 */
//!\}


// ----------------------------------------------------------------------------
// char_predicate
// ----------------------------------------------------------------------------

//!\cond
template <CharPredicate... condition_ts>
    requires sizeof...(condition_ts) >= 2
struct char_predicate_combiner;

template <CharPredicate condition_t>
struct char_predicate_negator;
//!\endcond

/*!\brief An abstract [CRTP](https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern) base class for
 *        parse conditions to add logical disjunction and negation operator.
 * \ingroup stream
 * \tparam derived_t The parse condition type to be extended with the logical operators.
 *                   Must model seqan3::detail::CharPredicate.
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
    //!\brief Combines the result of two seqan3::detail::CharPredicate via logical disjunction.
    template <CharPredicate rhs_t>
    constexpr auto operator||(rhs_t const &) const
    {
        return char_predicate_combiner<derived_t, rhs_t>{};
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
    template <std::Integral value_t>
    constexpr bool operator()(value_t const val) const noexcept
        requires sizeof(value_t) == 1
    {
        return derived_t::data[static_cast<unsigned char>(val)];
    }
    template <std::Integral value_t>
    constexpr bool operator()(value_t const val) const noexcept
        requires sizeof(value_t) != 1
    {
        return (static_cast<std::make_unsigned_t<value_t>>(val) < 256) ? operator()(static_cast<uint8_t>(val)) :
               (static_cast<decltype(EOF)>(val) == EOF)                ? derived_t::data[256]                  : false;
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
// char_predicate_combiner
// ----------------------------------------------------------------------------

/*!\brief Logical disjunction operator for parse conditions.
 * \implements seqan3::detail::CharPredicate
 * \tparam condition_ts Template parameter pack over all parse condition types. Must contain at least 2 template parameters.
 *                      Must model seqan3::detail::CharPredicate.
 * \ingroup stream
 */
template <CharPredicate... condition_ts>
//!\cond
    requires sizeof...(condition_ts) >= 2
//!\endcond
struct char_predicate_combiner : public char_predicate_base<char_predicate_combiner<condition_ts...>>
{
    //!\brief The message representing the disjunction of the associated conditions.
    static constexpr auto msg = detail::condition_message_v<'|', condition_ts...>;

    //!\brief The base type.
    using base_t = char_predicate_base<char_predicate_combiner<condition_ts...>>;

    //!\brief Import the data type from the base class.
    using typename base_t::data_t;
    //!\brief The look-up table that is used to evaluate the input.
    static constexpr data_t data = (condition_ts::data | ...);
};

/*!\brief Logical not operator for a parse condition.
 * \implements seqan3::detail::CharPredicate
 * \tparam condition_t Template parameter to apply the not-operator for.
 *                     Must model seqan3::detail::CharPredicate.
 * \ingroup stream
 */
template <CharPredicate condition_t>
struct char_predicate_negator : public char_predicate_base<char_predicate_negator<condition_t>>
{
    //!\brief The message representing the negation of the associated condition.
    static constexpr auto msg = small_string{'!'} + condition_t::msg;

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
 * \ingroup stream
 * \implements seqan3::detail::CharPredicate
 * \tparam interval_first non-type template parameter denoting the begin of the allowed range.
 *                        Must be less than or equal to `interval_last`.
 * \tparam interval_last non-type template parameter denoting the end of the allowed range.
 *                       Must be greater than or equal to `interval_first`.
 */
template <uint8_t interval_first, uint8_t interval_last>
//!\cond
    requires interval_first <= interval_last
//!\endcond
struct is_in_interval_type : public char_predicate_base<is_in_interval_type<interval_first, interval_last>>
{
    //!\brief The message representing this condition.
    static constexpr small_string msg = small_string{"is_in_interval<'"} +
                                        small_string{interval_first}     +
                                        small_string{"', '"}             +
                                        small_string{interval_last}      +
                                        small_string{"'>"};

    //!\brief The base type.
    using base_t = char_predicate_base<is_in_interval_type<interval_first, interval_last>>;

    //!\brief Import the data type from the base class.
    using typename base_t::data_t;
    //!\brief The look-up table that is used to evaluate the input.
    static constexpr data_t data = [] () constexpr
    {
        data_t ret{};

        for (size_t i = interval_first; i <= static_cast<size_t>(interval_last); ++i)
            ret[i] = true;

        return ret;
    }();
};

// ----------------------------------------------------------------------------
// is_in_alphabet_type
// ----------------------------------------------------------------------------

/*!\brief Parse condition that checks if a given value is within the given alphabet `alphabet_t`.
 * \ingroup stream
 * \implements seqan3::detail::CharPredicate
 * \tparam alphabet_t The alphabet type. Must model seqan3::Alphabet.
 */
template <detail::ConstexprAlphabet alphabet_t>
struct is_in_alphabet_type : public char_predicate_base<is_in_alphabet_type<alphabet_t>>
{
public:
    //!\brief The message representing this condition.
    static constexpr auto msg = small_string{"is_in_alphabet<"} +
                                small_string{detail::get_display_name_v<alphabet_t>} +
                                small_string{">"};

    //!\brief The base type.
    using base_t = char_predicate_base<is_in_alphabet_type<alphabet_t>>;

    //!\brief Import the data type from the base class.
    using typename base_t::data_t;
    //!\brief The look-up table that is used to evaluate the input.
    static constexpr data_t data = [] () constexpr
    {
        data_t ret{};

        for (size_t i = 0; i < 256; ++i)
            ret[i] = char_is_valid_for<alphabet_t>(static_cast<uint8_t>(i));

        return ret;
    }();
};

// ----------------------------------------------------------------------------
// is_char_type
// ----------------------------------------------------------------------------

/*!\brief Parse condition that checks if a given value is equal to `char_v`.
 * \ingroup stream
 * \implements seqan3::detail::CharPredicate
 * \tparam alphabet_t non-type template parameter with the value that should be checked against.
 */
template <int char_v>
struct is_char_type : public char_predicate_base<is_char_type<char_v>>
{
    static_assert(char_v == EOF || static_cast<uint64_t>(char_v) < 256, "TODO");

    //!\brief The message representing this condition.
    static constexpr auto msg = small_string{"is_char<'"} +
                                small_string{char_v}      +
                                small_string("'>");



    //!\brief The base type.
    using base_t = char_predicate_base<is_char_type<char_v>>;

    //!\brief Import the data type from the base class.
    using typename base_t::data_t;
    //!\brief The look-up table that is used to evaluate the input.
    static constexpr data_t data = [] () constexpr
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
