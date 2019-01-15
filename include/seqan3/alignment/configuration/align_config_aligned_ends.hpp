// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::align_cfg::aligned_ends.
 * \author Jörg Winkler <j.winkler AT fu-berlin.de>
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <type_traits>
#include <tuple>

#include <seqan3/alignment/configuration/detail.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/core/metafunction/template_inspection.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>
#include <seqan3/core/algorithm/parameter_pack.hpp>

namespace seqan3::align_cfg
{

/*!\brief A mixin class which can maintain a static or a dynamic bool state.
 * \ingroup configuration
 * \tparam value_t    The value type to be represented. Must be of type std::true_type, std::false_type or `bool`.
 * \tparam _is_static A boolean that evaluates to true only if `value_t` is a std::integral_constant.
 * \tparam _value     A boolean that captures the value of `value_t`, iff `_is_static` is `true`. Otherwise false.
 *
 * \details
 *
 * This mixin base class provides an optional pattern regarding the static state of the represented value.
 * If the mixin is constructed from an std::integral_constant it will hold an static state of the wrapped value.
 * In the other case, when constructing it from a boolean, the state of the value will be dynamic.
 */
template <typename value_t,
          bool _is_static  = std::conditional_t<std::is_same_v<value_t, bool>, std::false_type, std::true_type>::value,
          bool _value      = std::conditional_t<_is_static, value_t, std::false_type>::value>
//!\cond
    requires std::Same<value_t, std::true_type> || std::Same<value_t, std::false_type> || std::Same<value_t, bool>
//!\endcond
struct seq_end_gap_base
{
protected:

    //!\brief Friend of this class to provide access to static constexpr members.
    template <typename ...ends_t>
    friend class end_gaps;

    //!\brief Used to differentiate between static and dynamic state.
    static constexpr bool is_static    = _is_static;
    //!\brief Holds the static value if the state is static.
    static constexpr bool static_value = _value;

public:

    //!\brief Returns the wrapped value.
    constexpr bool operator()() const noexcept
    {
        return value;
    }

    //!\brief The wrapped value.
    value_t value{};
};

// ----------------------------------------------------------------------------
// first_seq_leading
// ----------------------------------------------------------------------------

/*!\brief The penalty configuration for a leading gap in the first sequence of a pairwise sequence alignment.
 * \ingroup configuration
 * \tparam value_t The type of the value to be wrapped. Can be of type std::true_type, std::false_type or bool.
 *
 * \details
 *
 * This configuration element enables/disables penalties for the respective gaps in the pairwise
 * sequence alignment. If the wrapped value evaluates to true, leading gaps in the first sequence will be not penalized.
 * If one constructs this element with a std::integral_constant it will convert to a static type such that
 * compile time optimizations can be enabled. If the type is constructed from a bool it will convert to a dynamic type
 * but will be later during the instantiation of the pairwise alignment converted to a static type.
 * Using a `bool` allows to dynamically set the value if the option is only known at runtime. If the option is
 * already known at compile time the static version will be a better option.
 *
 * \see first_seq_trailing
 * \see second_seq_leading
 * \see second_seq_trailing
 */
template <typename value_t>
struct first_seq_leading : public seq_end_gap_base<value_t>
{
    //!\privatesection
    //!\brief An internal id to allow consistency checks with other gap specifiers.
    static constexpr std::integral_constant<uint8_t, 0> id{};
};

/*!\name Type deduction guides
 * \relates seqan3::align_cfg::first_seq_leading
 * \{
 */
//!\brief Deduces the template argument from the type of the wrapped value.
template <typename value_t>
first_seq_leading(value_t) -> first_seq_leading<value_t>;
//!\}

// ----------------------------------------------------------------------------
// first_seq_trailing
// ----------------------------------------------------------------------------

/*!\brief The penalty configuration for a trailing gap in the first sequence of a pairwise sequence alignment.
 * \ingroup configuration
 * \tparam value_t The type of the value to be wrapped. Can be of type std::true_type, std::false_type or bool.
 *
 * \copydetails seqan3::align_cfg::first_seq_leading
 *
 * \see first_seq_leading
 * \see second_seq_leading
 * \see second_seq_trailing
 */
template <typename value_t>
struct first_seq_trailing : public seq_end_gap_base<value_t>
{
    //!\privatesection
    //!\brief An internal id to allow consistency checks with other gap specifiers.
    static constexpr std::integral_constant<uint8_t, 1> id{};
};

/*!\name Type deduction guides
 * \relates seqan3::align_cfg::first_seq_trailing
 * \{
 */
//!\brief Deduces the template argument from the type of the wrapped value.
template <typename value_t>
first_seq_trailing(value_t) -> first_seq_trailing<value_t>;
//!\}

// ----------------------------------------------------------------------------
// second_seq_leading
// ----------------------------------------------------------------------------

/*!\brief The penalty configuration for a leading gap in the second sequence of a pairwise sequence alignment.
 * \ingroup configuration
 * \tparam value_t The type of the value to be wrapped. Can be of type std::true_type, std::false_type or bool.
 *
 * \copydetails seqan3::align_cfg::first_seq_leading
 *
 * \see first_seq_leading
 * \see first_seq_trailing
 * \see second_seq_trailing
 */
template <typename value_t>
struct second_seq_leading : public seq_end_gap_base<value_t>
{
    //!\privatesection
    //!\brief An internal id to allow consistency checks with other gap specifiers.
    static constexpr std::integral_constant<uint8_t, 2> id{};
};

/*!\name Type deduction guides
 * \relates seqan3::align_cfg::second_seq_leading
 * \{
 */
//!\brief Deduces the template argument from the type of the wrapped value.
template <typename value_t>
second_seq_leading(value_t) -> second_seq_leading<value_t>;
//!\}

// ----------------------------------------------------------------------------
// second_seq_trailing
// ----------------------------------------------------------------------------

/*!\brief The penalty configuration for a trailing gap in the second sequence of a pairwise sequence alignment.
 * \ingroup configuration
 * \tparam value_t The type of the value to be wrapped. Can be of type std::true_type, std::false_type or bool.
 *
 * \copydetails seqan3::align_cfg::first_seq_leading
 *
 * \see first_seq_leading
 * \see first_seq_trailing
 * \see second_seq_leading
 */
template <typename value_t>
struct second_seq_trailing : public seq_end_gap_base<value_t>
{
    //!\privatesection
    //!\brief An internal id to allow consistency checks with other gap specifiers.
    static constexpr std::integral_constant<uint8_t, 3> id{};
};

/*!\name Type deduction guides
 * \relates seqan3::align_cfg::second_seq_trailing
 * \{
 */
//!\brief Deduces the template argument from the type of the wrapped value.
template <typename value_t>
second_seq_trailing(value_t) -> second_seq_trailing<value_t>;
//!\}

// ----------------------------------------------------------------------------
// end_gaps
// ----------------------------------------------------------------------------

/*!\brief Wraps the sequence end-gap specifiers and provides ordered access to the respective values.
 * \ingroup configuration
 * \tparam ends_t A parameter pack containing at most 4 sequence end-gap specifier.
 *
 * \details
 *
 * A wrapper for providing ordered access to the end-gap specifiers independent of the input order.
 * The possible input types can be: seqan3::first_seq_leading, seqan3::first_seq_trailing, seqan3::second_seq_leading and
 * seqan3::second_seq_trailing.
 * The types in the parameter pack `ends_t` are deduced by the corresponding constructor argument.
 * If a specifier is not set it will default to `false` and thus the respective end-gap will be penalized in the
 * pairwise alignment.
 */
template <typename ...ends_t>
//!\cond
    requires sizeof...(ends_t) <= 4 &&
             ((detail::is_type_specialisation_of_v<ends_t, first_seq_leading> ||
               detail::is_type_specialisation_of_v<ends_t, first_seq_trailing> ||
               detail::is_type_specialisation_of_v<ends_t, second_seq_leading>  ||
               detail::is_type_specialisation_of_v<ends_t, second_seq_trailing>) && ...)
//!\endcond
class end_gaps
{
    //!\brief Helper function to check valid end_gaps configuration.
    template <typename ..._ends_t>
    static constexpr bool check_consistency(_ends_t ...ends)
    {
        if constexpr (sizeof...(ends) < 2)
        {
            return true;
        }
        else
        {
            return [] (auto head, auto ...tail) constexpr
            {
                using head_t = decltype(head);
                if constexpr (((head_t::id != decltype(tail)::id) && ...))
                    return check_consistency(tail...);
                else
                    return false;
            }(ends...);
        }
    }

    static_assert(check_consistency(ends_t{}...),
                  "You may not use the same end_gap specifier more than once.");

public:
    /*!\name Constructor, destructor and assignment
     * \{
     */
    constexpr end_gaps() noexcept
    {
        [[maybe_unused]] auto dummy = ((values[std::remove_reference_t<ends_t>::id()] =
            std::remove_reference_t<ends_t>::static_value), ..., 0);
    }

    constexpr end_gaps(end_gaps const &)             noexcept = default;
    constexpr end_gaps(end_gaps &&)                  noexcept = default;
    constexpr end_gaps & operator=(end_gaps const &) noexcept = default;
    constexpr end_gaps & operator=(end_gaps &&)      noexcept = default;
    ~end_gaps()                                      noexcept = default;

    //!\brief Construction from at least one sequence end-gap specifier.
    constexpr end_gaps(ends_t && ...args) noexcept
        requires sizeof...(ends_t) > 0
    {
        detail::for_each_value([this](auto e)
        {
            values[remove_cvref_t<decltype(e)>::id()] = e();
        }, std::forward<ends_t>(args)...);
    }
    //!}

    /*!\name Element access
     * \{
     */
    /*!\brief Returns the value for the specifier at the given position.
     * \param[in] pos The position to get the value for the respective end-gap for.
     *
     * \details
     *
     * The sequence end-gap specifier are stored in an ordered fashion. The following position mapping will be used
     * to access the respective values:
     * seqan3::first_seq_leading &rarr; 0; seqan3::first_seq_trailing &rarr; 1; seqan3::second_seq_leading &rarr; 2;
     * seqan3::second_seq_trailing &rarr; 3.
     *
     * \returns `true` if the respective sequence end-gap is set to be free, `false` otherwise.
     */
    constexpr bool operator[](size_t const pos) const noexcept
    {
        assert(pos < values.size());
        return values[pos];
    }

    /*!\brief Returns the static value for the specifier at the given position.
     * \tparam pos The position to get the value for the respective end-gap for.
     *
     * \details
     *
     * The sequence end-gap specifier are stored in an ordered fashion. The following position mapping will be used
     * to access the respective values:
     * seqan3::first_seq_leading &rarr; 0; seqan3::first_seq_trailing &rarr; 1; seqan3::second_seq_leading &rarr; 2;
     * seqan3::second_seq_trailing &rarr; 3.
     *
     * \returns `true` if the respective sequence end-gap is set to be free, `false` otherwise.
     */
    template <size_t pos>
    static constexpr bool get_static() noexcept
    {
        static_assert(is_static_array[pos],
                      "You may not access an element that was not set in a core constant expression.");
        return get<pos>(static_values);
    }
    //!\}

    /*!\ name Observers
     * \{
     */
    /*!\brief Returns whether a value at the given position was set statically.
     * \tparam pos The position to get the value for the respective end-gap for.
     *
     * \details
     *
     * \copydetails end_gaps::get_static
     *
     * \returns `true` if the respective sequence end-gap was set in a static context, `false` otherwise.
     */
    template <size_t pos>
    static constexpr bool is_static() noexcept
    {
        return get<pos>(is_static_array);
    }
    //!\}

private:

    //!\brief Stores the values.
    std::array<bool, 4> values{false, false, false, false};

    //!\brief Stores whether a value is accessible in a constexpr context.
    static constexpr std::array<bool, 4> is_static_array
    {
        [](auto ...ends) constexpr
        {
            std::array<bool, 4> tmp{false, false, false, false};
            detail::for_each_value([&tmp](auto v)
            {
                tmp[decltype(v)::id()] = decltype(v)::is_static;
            }, ends...);
            return tmp;
        }(ends_t{}...)
    };

    //!\brief Stores the static values.
    static constexpr std::array<bool, 4> static_values
    {
        [](auto ...ends) constexpr
        {
            std::array<bool, 4> tmp{false, false, false, false};
            detail::for_each_value([&tmp](auto v)
            {
                tmp[decltype(v)::id()] = decltype(v)::static_value;
            }, ends...);
            return tmp;
        }(ends_t{}...)
    };
};

/*!\name Type deduction guides
 * \relates seqan3::align_cfg::end_gaps
 * \{
 */

//!\brief Deduces the end-gap specifier from the constructor arguments.
template <typename ... ends_t>
end_gaps(ends_t && ...) -> end_gaps<std::remove_reference_t<ends_t>...>;

//!\}

// ----------------------------------------------------------------------------
// aligned_ends
// ----------------------------------------------------------------------------

/*!\brief The configuration for aligned sequence ends.
 * \ingroup configuration
 * \tparam end_gaps_t The type of the end gaps. Must be of type seqan3::end_gaps.
 */
template <typename end_gaps_t>
//!\cond
    requires detail::is_type_specialisation_of_v<end_gaps_t, end_gaps>
//!\endcond
struct aligned_ends : public pipeable_config_element<aligned_ends<end_gaps_t>, end_gaps_t>
{
    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::align_config_id id{detail::align_config_id::aligned_ends};
};

/*!\name Type deduction guides
 * \relates seqan3::align_cfg::aligned_ends
 * \{
 */
//!\brief Deduces the end gaps object type from the constructor argument.
template <typename end_gaps_t>
aligned_ends(end_gaps_t) -> aligned_ends<std::remove_reference_t<end_gaps_t>>;
//!\}

// ----------------------------------------------------------------------------
// all_ends_free
// ----------------------------------------------------------------------------

/*!\name Fixed end gaps configurations
 * \brief These variables are pre-configured end gaps that are frequently used in pairwise sequence alignments.
 * \relates seqan3::end_gaps
 */

//!\brief Disables all end-gaps, such that they are not penalized in the pairwise alignment.
inline constexpr end_gaps all_ends_free{first_seq_leading<std::true_type>{},
                                        first_seq_trailing<std::true_type>{},
                                        second_seq_leading<std::true_type>{},
                                        second_seq_trailing<std::true_type>{}};

// ----------------------------------------------------------------------------
// none_ends_free
// ----------------------------------------------------------------------------

//!\brief Enables all end-gaps, such that they are penalized in the pairwise alignment.
inline constexpr end_gaps none_ends_free{first_seq_leading<std::false_type>{},
                                         first_seq_trailing<std::false_type>{},
                                         second_seq_leading<std::false_type>{},
                                         second_seq_trailing<std::false_type>{}};

// ----------------------------------------------------------------------------
// seq1_ends_free
// ----------------------------------------------------------------------------

//!\brief Enables end-gaps for the first sequence, to compute a semi-global alignment
//!       with free end gaps in the first sequence
inline constexpr end_gaps seq1_ends_free{first_seq_leading<std::true_type>{},
                                         first_seq_trailing<std::true_type>{},
                                         second_seq_leading<std::false_type>{},
                                         second_seq_trailing<std::false_type>{}};

// ----------------------------------------------------------------------------
// seq2_ends_free
// ----------------------------------------------------------------------------

//!\brief Enables end-gaps for the second sequence, to compute a semi-global alignment
//!       with free end gaps in the second sequence
inline constexpr end_gaps seq2_ends_free{first_seq_leading<std::false_type>{},
                                         first_seq_trailing<std::false_type>{},
                                         second_seq_leading<std::true_type>{},
                                         second_seq_trailing<std::true_type>{}};
//!\}
} // namespace seqan3
