// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides basic data structure for strong types.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <seqan3/core/add_enum_bitwise_operators.hpp>
#include <seqan3/core/type_traits/basic.hpp>
#include <seqan3/std/concepts>

namespace seqan3::detail
{

/*!\brief Enum class for all supported operations that can be added to a seqan3::detail::strong_type.
 * \ingroup core
 */
enum struct strong_type_skill
{
    none           = 0,
    add            = 1 << 0,
    subtract       = 1 << 1,
    multiply       = 1 << 2,
    divide         = 1 << 3,
    modulo         = 1 << 4,
    bitwise_and    = 1 << 5,
    bitwise_or     = 1 << 6,
    bitwise_xor    = 1 << 7,
    bitwise_not    = 1 << 8,
    bitwise_lshift = 1 << 9,
    bitwise_rshift = 1 << 10,
    logical_and    = 1 << 11,
    logical_or     = 1 << 12,
    logical_not    = 1 << 13,
    increment      = 1 << 14,
    decrement      = 1 << 15,
    convert        = 1 << 16,
    additive       = add | subtract,
    multiplicative = multiply | divide | modulo,
    bitwise_logic  = bitwise_and | bitwise_or | bitwise_xor | bitwise_not,
    bitwise_shift  = bitwise_lshift | bitwise_rshift,
    logic          = logical_and | logical_or | logical_not
};
}  //namespace seqan3::detail

namespace seqan3
{
//!\cond
template <>
constexpr bool add_enum_bitwise_operators<seqan3::detail::strong_type_skill> = true;
//!\endcond
}

namespace seqan3::detail
{

/*!\brief CRTP base class to declare a strong typedef for a regular type to avoid ambiguous parameter settings in function
 *        calls.
 * \ingroup core
 * \tparam value_t   The underlying type to create a strong typedef for.
 * \tparam derived_t The derived class inheriting from this base class. [see CRTP](https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern).
 * \tparam skills    A set of skills to be add to the expressive type.
 *
 * \details
 *
 * There are many cases in which interfaces use ambiguous parameters which can be mixed up very easily.
 * For example, if we consider an interface, that expects a window size and a number of maximal errors to search for
 * possible hits in a region of interest, then both values might be given in form of an unsigned integer.
 * The following snippet shows a typical interface:
 *
 * \include test/snippet/core/detail/strong_type_usage.cpp
 *
 * The first parameter is the window size and the last parameter defines the error threshold.
 * But, what happens if the user accidentally switches the `window_size` with the `error` parameter?
 * In a bigger code base this mistake can be so subtle, that it becomes hard to figure out that the
 * search interface was used in a wrong way. Also the compiler cannot detect this.
 * In order to make this interface safe, we have to use a _strong type_.
 * A strong type is expressive in what it actually represents as a value.
 * In our toy example we could define two strong types as follows:
 *
 * \include test/snippet/core/detail/strong_type_error_window.cpp
 * Our interface could now be changed to:
 *
 * \include test/snippet/core/detail/strong_type_new_usage.cpp
 *
 * Now the user is forced to pass the parameters as their named type. If the parameter order is mixed up by accident
 * the compiler would emit an error message, since the `error` type is not convertible to the `window_size` type and
 * vice versa.
 *
 * ### Adding Skills
 *
 * In most cases it is sufficient to protect the user interface from misuse by using strong types. Within the
 * implementation the underlying value can then be extracted using the `getter`-member functions of the
 * seqan3::detail::strong_type class. However, there might be scenarios, where the strong type is exposed to the user
 * to work with, but with a restricted set of operations that can be applied to the type.
 * A familiar use case are the time typedefs of the [std::chrono](https://en.cppreference.com/w/cpp/header/chrono)
 * library. It is convenient for the user to subtract two time values representing, for example seconds, to get the
 * duration between two time points. But not all operations like modulo or multiplication are really useful for this.
 * In order to add skills to the seqan3::detail::strong_type the typedef can be further specialized with
 * operations from the seqan3::detail::strong_type_skill enum.
 * For example, we could further specify our error type to support increment and decrement operations.
 *
 * \include test/snippet/core/detail/strong_type_adding_skills.cpp
 */
template <typename value_t, typename derived_t, strong_type_skill skills = strong_type_skill::none>
class strong_type
{
public:

    //!\brief The underlying value type.
    using value_type = value_t;

    /*!\name Constructor, destructor and assignment.
     * \brief The standard functions are explicitly defaulted.
     * \{
     */
    constexpr strong_type()                                 noexcept = default; //!< Defaulted.
    constexpr strong_type(strong_type const &)              noexcept = default; //!< Defaulted.
    constexpr strong_type(strong_type &&)                   noexcept = default; //!< Defaulted.
    constexpr strong_type & operator= (strong_type const &) noexcept = default; //!< Defaulted.
    constexpr strong_type & operator= (strong_type &&)      noexcept = default; //!< Defaulted.
    ~strong_type()                                          noexcept = default; //!< Defaulted.

    //!\brief Construction from underlying value type.
    constexpr explicit strong_type(value_t _value) : value(std::move(_value))
    {}
    //!\}

    /*!\name Accessor.
    * \{
    */
    //!\brief Returns the underlying value.
    constexpr value_t & get() & noexcept
    {
        return value;
    }

    //!\brief Returns the underlying value.
    constexpr value_t const & get() const & noexcept
    {
        return value;
    }

    //!\brief Returns the underlying value as rvalue.
    constexpr value_t && get() && noexcept
    {
        return std::move(value);
    }

    //!\brief Returns the underlying value as rvalue.
    constexpr value_t const && get() const && noexcept
    {
        return std::move(value);
    }
    //!\}

    /*!\name Arithmetic additive operators.
     * \brief Only available if the corresponding skills from seqan3::detail::strong_type_skill are added and
     *        the underlying type supports this operation.
     * \{
     */
    //!\brief Adds addition operator to the strong type.
    constexpr derived_t operator+(strong_type const & other)
        //!\cond
        requires ((skills & strong_type_skill::add) != strong_type_skill::none)
        //!\endcond
    {
        return derived_t{get() + other.get()};
    }

    //!\brief Adds subtraction operator to the strong type.
    constexpr derived_t operator-(strong_type const & other)
        //!\cond
        requires ((skills & strong_type_skill::subtract) != strong_type_skill::none)
        //!\endcond
    {
        return derived_t{get() - other.get()};
    }
    //!\}

    /*!\name Arithmetic multiplicative operators.
     * \brief Only available if the corresponding skills from seqan3::detail::strong_type_skill are added and
     *        the underlying type supports this operation.
     * \{
     */
    //!\brief Adds multiplication operator to the strong type.
    constexpr derived_t operator*(strong_type const & other)
        //!\cond
        requires ((skills & strong_type_skill::multiply) != strong_type_skill::none)
        //!\endcond
    {
        return derived_t{get() * other.get()};
    }

    //!\brief Adds division operator to the strong type.
    constexpr derived_t operator/(strong_type const & other)
        //!\cond
        requires ((skills & strong_type_skill::divide) != strong_type_skill::none)
        //!\endcond
    {
        return derived_t{get() / other.get()};
    }

    //!\brief Adds modulo operator to the strong type.
    constexpr derived_t operator%(strong_type const & other)
        //!\cond
        requires ((skills & strong_type_skill::modulo) != strong_type_skill::none)
        //!\endcond
    {
        return derived_t{get() % other.get()};
    }
    //!\}

    /*!\name Bitwise logic operators.
     * \brief Only available if the corresponding skills from seqan3::detail::strong_type_skill are added and
     *        the underlying type supports this operation.
     * \{
     */

    //!\brief Adds bitwise and operator to the strong type.
    constexpr derived_t operator&(strong_type const & other)
        //!\cond
        requires ((skills & strong_type_skill::bitwise_and) != strong_type_skill::none)
        //!\endcond
    {
        return derived_t{get() & other.get()};
    }

    //!\brief Adds bitwise or operator to the strong type.
    constexpr derived_t operator|(strong_type const & other)
        //!\cond
        requires ((skills & strong_type_skill::bitwise_or) != strong_type_skill::none)
        //!\endcond
    {
        return derived_t{get() | other.get()};
    }

    //!\brief Adds bitwise xor operator to the strong type.
    constexpr derived_t operator^(strong_type const & other)
        //!\cond
        requires ((skills & strong_type_skill::bitwise_xor) != strong_type_skill::none)
        //!\endcond
    {
        return derived_t{get() ^ other.get()};
    }

    //!\brief Adds bitwise not operator to the strong type.
    constexpr derived_t operator~()
        //!\cond
        requires ((skills & strong_type_skill::bitwise_not) != strong_type_skill::none)
        //!\endcond
    {
        return derived_t{~get()};
    }
    //!\}

    /*!\name Bitwise shift operators.
     * \brief Only available if the corresponding skills from seqan3::detail::strong_type_skill are added and
     *        the underlying type supports this operation.
     * \{
     */

    //!\brief Adds bitwise left shift operator to the strong type.
    constexpr derived_t operator<<(strong_type const & other)
        //!\cond
        requires ((skills & strong_type_skill::bitwise_lshift) != strong_type_skill::none)
        //!\endcond
    {
        return derived_t{get() << other.get()};
    }

    //!\brief Adds bitwise left shift operator to the strong type.
    template <std::integral integral_t>
    constexpr derived_t operator<<(integral_t const shift)
        //!\cond
        requires ((skills & strong_type_skill::bitwise_lshift) != strong_type_skill::none)
        //!\endcond
    {
        return derived_t{get() << shift};
    }

    //!\brief Adds bitwise right shift operator to the strong type.
    constexpr derived_t operator>>(strong_type const & other)
        //!\cond
        requires ((skills & strong_type_skill::bitwise_rshift) != strong_type_skill::none)
        //!\endcond
    {
        return derived_t{get() >> other.get()};
    }

    //!\brief Adds bitwise right shift operator to the strong type.
    template <std::integral integral_t>
    constexpr derived_t operator>>(integral_t const shift)
        //!\cond
        requires ((skills & strong_type_skill::bitwise_rshift) != strong_type_skill::none)
        //!\endcond
    {
        return derived_t{get() >> shift};
    }
    //!\}

    /*!\name Logical operators.
     * \brief Only available if the corresponding skills from seqan3::detail::strong_type_skill are added and
     *        the underlying type supports this operation.
     * \{
     */

    //!\brief Adds logical and operator to the strong type.
    constexpr bool operator&&(strong_type const & other)
        //!\cond
        requires ((skills & strong_type_skill::logical_and) != strong_type_skill::none)
        //!\endcond
    {
        return get() && other.get();
    }

    //!\brief Adds logical or operator to the strong type.
    constexpr bool operator||(strong_type const & other)
        //!\cond
        requires ((skills & strong_type_skill::logical_or) != strong_type_skill::none)
        //!\endcond
    {
        return get() || other.get();
    }

    //!\brief Adds logical not operator to the strong type.
    constexpr bool operator!()
        //!\cond
        requires ((skills & strong_type_skill::logical_not) != strong_type_skill::none)
        //!\endcond
    {
        return !get();
    }
    //!\}

    /*!\name Increment and decrement operators.
     * \brief Only available if the corresponding skills from seqan3::detail::strong_type_skill are added and
     *        the underlying type supports this operation.
     * \{
     */
    //!\brief Adds pre-increment operator to the strong type.
    constexpr derived_t & operator++()
        //!\cond
        requires ((skills & strong_type_skill::increment) != strong_type_skill::none)
        //!\endcond
    {
        ++get();
        return static_cast<derived_t &>(*this);
    }

    //!\brief Adds post-increment operator to the strong type.
    constexpr derived_t operator++(int)
        //!\cond
        requires ((skills & strong_type_skill::increment) != strong_type_skill::none)
        //!\endcond
    {
        derived_t tmp{get()};
        ++get();
        return tmp;
    }

    //!\brief Adds pre-decrement operator to the strong type.
    constexpr derived_t & operator--()
        //!\cond
        requires ((skills & strong_type_skill::decrement) != strong_type_skill::none)
        //!\endcond
    {
        --get();
        return static_cast<derived_t &>(*this);
    }

    //!\brief Adds post-decrement operator to the strong type.
    constexpr derived_t operator--(int)
        //!\cond
        requires ((skills & strong_type_skill::decrement) != strong_type_skill::none)
        //!\endcond
    {
        derived_t tmp{get()};
        --get();
        return tmp;
    }
    //!\}

    /*!\name Conversion operators.
     * \brief Only available if the corresponding skill from seqan3::detail::strong_type_skill is added.
     * \{
     */

    //!\brief Adds explicit conversion to it's underlying type.
    explicit constexpr operator value_t() const
        //!\cond
        requires ((skills & strong_type_skill::convert) != strong_type_skill::none)
        //!\endcond
    {
        return get();
    }
    //!\}
private:
    //!\brief The underlying value, which is wrapped as a strong type.
    value_t value;
};

} // namespace seqan3::detail
