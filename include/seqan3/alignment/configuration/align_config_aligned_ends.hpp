// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
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
#include <seqan3/core/algorithm/pipeable_config_element.hpp>
#include <seqan3/core/detail/pack_algorithm.hpp>
#include <seqan3/core/type_traits/basic.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>

namespace seqan3
{

/*!\brief A mixin class which can maintain a static or a dynamic bool state.
 * \ingroup alignment_configuration
 * \tparam value_t    The value type to be represented. Must be of type std::true_type, std::false_type or `bool`.
 * \tparam _is_static A boolean that evaluates to true only if `value_t` is a std::integral_constant.
 * \tparam _value     A boolean that captures the value of `value_t`, iff `_is_static` is `true`. Otherwise false.
 *
 * \details
 *
 * This mixin base class provides an optional pattern regarding the static state of the represented value.
 * If the mixin is constructed from an std::integral_constant, it will hold an static state of the wrapped value.
 * In the other case, when constructing it from a boolean, the state of the value will be dynamic.
 */
template <typename value_t,
          bool _is_static  = std::conditional_t<std::is_same_v<value_t, bool>, std::false_type, std::true_type>::value,
          bool _value      = std::conditional_t<_is_static, value_t, std::false_type>::value>
//!\cond
    requires std::same_as<value_t, std::true_type> || std::same_as<value_t, std::false_type> || std::same_as<value_t, bool>
//!\endcond
struct sequence_end_gap_specifier_base
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
// front_end_first
// ----------------------------------------------------------------------------

/*!\brief The penalty configuration for aligning the front of the first sequence with a gap.
 * \ingroup alignment_configuration
 * \tparam value_t The type of the value to be wrapped. Can be of type std::true_type, std::false_type or bool.
 *
 * \details
 *
 * This strong type enables (`false`) or disables (`true`) penalties for aligning the respective sequence end with gaps.
 * If one constructs this element with a std::integral_constant it will convert to a static type such that
 * compile time optimisations can be used. If the type is constructed from a bool it will convert to a dynamic type
 * but will be converted to a static type during the configuration of the pairwise alignment algorithm.
 * Using a `bool` allows to dynamically set the value if the option is only known at runtime. If the option is
 * already known at compile time the static version will be the preferred option.
 *
 * \see seqan3::back_end_first
 * \see seqan3::front_end_second
 * \see seqan3::back_end_second
 */
template <typename value_t>
struct front_end_first : public sequence_end_gap_specifier_base<value_t>
{
    //!\privatesection
    //!\brief An internal id to allow consistency checks with other gap specifiers.
    static constexpr std::integral_constant<uint8_t, 0> id{};
};

/*!\name Type deduction guides
 * \relates seqan3::front_end_first
 * \{
 */
//!\brief Deduces the template argument from the type of the wrapped value.
template <typename value_t>
front_end_first(value_t) -> front_end_first<value_t>;
//!\}

// ----------------------------------------------------------------------------
// back_end_first
// ----------------------------------------------------------------------------

/*!\brief The penalty configuration for aligning the back of the first sequence with a gap.
 * \ingroup alignment_configuration
 * \tparam value_t The type of the value to be wrapped. Can be of type std::true_type, std::false_type or bool.
 *
 * \copydetails seqan3::front_end_first
 *
 * \see front_end_first
 * \see front_end_second
 * \see back_end_second
 */
template <typename value_t>
struct back_end_first : public sequence_end_gap_specifier_base<value_t>
{
    //!\privatesection
    //!\brief An internal id to allow consistency checks with other gap specifiers.
    static constexpr std::integral_constant<uint8_t, 1> id{};
};

/*!\name Type deduction guides
 * \relates seqan3::back_end_first
 * \{
 */
//!\brief Deduces the template argument from the type of the wrapped value.
template <typename value_t>
back_end_first(value_t) -> back_end_first<value_t>;
//!\}

// ----------------------------------------------------------------------------
// front_end_second
// ----------------------------------------------------------------------------

/*!\brief The penalty configuration for aligning the front of the second sequence with a gap.
 * \ingroup alignment_configuration
 * \tparam value_t The type of the value to be wrapped. Can be of type std::true_type, std::false_type or bool.
 *
 * \copydetails seqan3::front_end_first
 *
 * \see front_end_first
 * \see back_end_first
 * \see back_end_second
 */
template <typename value_t>
struct front_end_second : public sequence_end_gap_specifier_base<value_t>
{
    //!\privatesection
    //!\brief An internal id to allow consistency checks with other gap specifiers.
    static constexpr std::integral_constant<uint8_t, 2> id{};
};

/*!\name Type deduction guides
 * \relates seqan3::front_end_second
 * \{
 */
//!\brief Deduces the template argument from the type of the wrapped value.
template <typename value_t>
front_end_second(value_t) -> front_end_second<value_t>;
//!\}

// ----------------------------------------------------------------------------
// back_end_second
// ----------------------------------------------------------------------------

/*!\brief The penalty configuration for aligning the back of the second sequence with a gap.
 * \ingroup alignment_configuration
 * \tparam value_t The type of the value to be wrapped. Can be of type std::true_type, std::false_type or bool.
 *
 * \copydetails seqan3::front_end_first
 *
 * \see front_end_first
 * \see back_end_first
 * \see front_end_second
 */
template <typename value_t>
struct back_end_second : public sequence_end_gap_specifier_base<value_t>
{
    //!\privatesection
    //!\brief An internal id to allow consistency checks with other gap specifiers.
    static constexpr std::integral_constant<uint8_t, 3> id{};
};

/*!\name Type deduction guides
 * \relates seqan3::back_end_second
 * \{
 */
//!\brief Deduces the template argument from the type of the wrapped value.
template <typename value_t>
back_end_second(value_t) -> back_end_second<value_t>;
//!\}

// ----------------------------------------------------------------------------
// end_gaps
// ----------------------------------------------------------------------------

/*!\brief Wraps the sequence end-gap specifiers and provides ordered access to the respective values.
 * \ingroup alignment_configuration
 * \tparam ends_t A parameter pack containing at most 4 sequence end-gap specifier.
 *
 * \details
 *
 * A wrapper for providing ordered access to the end-gap specifiers independent of the input order.
 * The possible input types can be: seqan3::front_end_first, seqan3::back_end_first,
 * seqan3::front_end_second and seqan3::back_end_second.
 * The types in the parameter pack `ends_t` are deduced by the corresponding constructor argument.
 * If a specifier is not set it will default to `false` and thus the respective end-gap will be penalised in the
 * pairwise alignment.
 *
 * ### Static vs runtime configuration
 *
 * The end_gaps class preserves the static/non-static property of the respective end-gap specifier. Those specifiers
 * can, depending on how they are constructed, contain a static information or a runtime information whether or not a
 * specific end-gap is enabled. To check whether the information was static the function end_gaps::is_static can be
 * used. If it was static the function end_gaps::get_static can be used to obtain the respective value at compile time.
 * Note that both functions are static and thus need to be called on the type and not on an instance of the type.
 *
 * \include snippet/alignment/configuration/align_cfg_aligned_ends_access.cpp
 *
 * To get the respective value at runtime use the \ref operator[]() "[]-operator". This function returns always the
 * respective value independent of whether the value was provided by a static variable or a runtime variable.
 * Also note that static and non-static end-gap specifier can be mixed within the end_gaps class
 * template.
 *
 * It is strongly recommended to use the static information if possible and only make those specifiers depend on
 * runtime parameters that cannot be resolved at compile time. This will reduce the compile time while configuring
 * the alignment algorithm since the runtime information need to be translated into static information for the alignment
 * algorithm.
 *
 * \see seqan3::front_end_first
 * \see seqan3::back_end_first
 * \see seqan3::front_end_second
 * \see seqan3::back_end_second
 */
template <typename ...ends_t>
//!\cond
    requires sizeof...(ends_t) <= 4 &&
             ((detail::is_type_specialisation_of_v<ends_t, front_end_first> ||
               detail::is_type_specialisation_of_v<ends_t, back_end_first> ||
               detail::is_type_specialisation_of_v<ends_t, front_end_second>  ||
               detail::is_type_specialisation_of_v<ends_t, back_end_second>) && ...)
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
    /*!\name Constructors, destructor and assignment
     * \{
     */

    //! \brief Default constructor.
    constexpr end_gaps() noexcept
    {
        [[maybe_unused]] auto dummy = ((values[std::remove_reference_t<ends_t>::id()] =
            std::remove_reference_t<ends_t>::static_value), ..., 0);
    }

    constexpr end_gaps(end_gaps const &)             noexcept = default; //!< Defaulted
    constexpr end_gaps(end_gaps &&)                  noexcept = default; //!< Defaulted
    constexpr end_gaps & operator=(end_gaps const &) noexcept = default; //!< Defaulted
    constexpr end_gaps & operator=(end_gaps &&)      noexcept = default; //!< Defaulted
    ~end_gaps()                                      noexcept = default; //!< Defaulted

    //!\brief Construction from at least one sequence end-gap specifier.
    constexpr end_gaps(ends_t const ...args) noexcept
        requires sizeof...(ends_t) > 0
    {
        detail::for_each([this](auto e)
        {
            values[remove_cvref_t<decltype(e)>::id()] = e();
        }, args...);
    }
    //\!}

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
     * seqan3::front_end_first &rarr; 0; seqan3::back_end_first &rarr; 1;
     * seqan3::front_end_second &rarr; 2; seqan3::back_end_second &rarr; 3.
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
     * seqan3::front_end_first &rarr; 0; seqan3::back_end_first &rarr; 1;
     * seqan3::front_end_second &rarr; 2; seqan3::back_end_second &rarr; 3.
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
            detail::for_each([&tmp](auto v)
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
            detail::for_each([&tmp](auto v)
            {
                tmp[decltype(v)::id()] = decltype(v)::static_value;
            }, ends...);
            return tmp;
        }(ends_t{}...)
    };
};

/*!\name Type deduction guides
 * \{
 */

/*!\brief Deduces the end-gap specifier from the constructor arguments.
 * \relates seqan3::end_gaps
 * \tparam ends_t A template parameter pack containing at most 4 sequence end-gap specifiers.
 */
template <typename ...ends_t>
end_gaps(ends_t const & ...) -> end_gaps<ends_t...>;
//!\}

// ----------------------------------------------------------------------------
// free_ends_all
// ----------------------------------------------------------------------------

/*!\name Predefined end-gaps configurations
 * \relates seqan3::end_gaps
 * \anchor predefined_end_gap_configurations
 * \brief These variables are pre-configured end-gaps that are frequently used in pairwise sequence alignments.
 * \{
 */

/*!\brief All ends are free.
 *
 * \details
 *
 * Computes an overlap alignment where the end of one sequence can overlap with the end of the other sequence.
 * In the following example the gaps at the ends are not penalised and the sequences are aligned such that the prefix
 * of the first sequences matches the suffix of the second sequence.
 *
 * ```
 * -----ACGTAAAACGT
 *      |||||
 * TTTTTACGTA------
 * ```
 */
inline constexpr end_gaps free_ends_all{front_end_first<std::true_type>{},
                                        back_end_first<std::true_type>{},
                                        front_end_second<std::true_type>{},
                                        back_end_second<std::true_type>{}};

// ----------------------------------------------------------------------------
// free_ends_none
// ----------------------------------------------------------------------------

/*!\brief All ends are penalised.
 *
 * \details
 *
 * Computes a global alignment where all end-gaps are penalised. For example in the following alignment, the
 * alignment is forced to cover the entire sequences and the leading gaps will be penalised.
 *
 * ```
 * ---ACG--TAAAACGT
 *    |||  || | |||
 * AAAACGTATAGACCGT
 * ```
 */
inline constexpr end_gaps free_ends_none{front_end_first<std::false_type>{},
                                         back_end_first<std::false_type>{},
                                         front_end_second<std::false_type>{},
                                         back_end_second<std::false_type>{}};

// ----------------------------------------------------------------------------
// free_ends_first
// ----------------------------------------------------------------------------

/*!\brief Ends of the first sequence are free.
 *
 * \details
 *
 * Computes a semi-global alignment where the ends of the first sequence can align to gaps without additional costs.
 * For example in the following alignment, the leading and trailing gaps are not penalised and the smaller sequence
 * can be aligned such that it matches the middle part of the longer sequence.
 *
 * ```
 * TTTTTACGT---ATGTCCCCC
 *      ||||   | ||
 * -----ACGTAAAACGT-----
 * ```
 */
inline constexpr end_gaps free_ends_first{front_end_first<std::true_type>{},
                                          back_end_first<std::true_type>{},
                                          front_end_second<std::false_type>{},
                                          back_end_second<std::false_type>{}};

// ----------------------------------------------------------------------------
// free_ends_second
// ----------------------------------------------------------------------------

/*!\brief Ends for the second sequence are free.
 *
 * \details
 *
 * Computes a semi-global alignment where the ends of the second sequence can align to gaps without additional costs.
 * For example in the following alignment, the leading and trailing gaps are not penalised and the smaller sequence
 * can be aligned such that it matches the middle part of the longer sequence.
 *
 * ```
 * -----ACGTAAAACGT-----
 *      ||||   | ||
 * TTTTTACGT---ATGTCCCCC
 * ```
 */
inline constexpr end_gaps free_ends_second{front_end_first<std::false_type>{},
                                           back_end_first<std::false_type>{},
                                           front_end_second<std::true_type>{},
                                           back_end_second<std::true_type>{}};
//!\}
} // namespace seqan3

namespace seqan3::align_cfg
{

// ----------------------------------------------------------------------------
// aligned_ends
// ----------------------------------------------------------------------------

/*!\brief The configuration for aligned sequence ends.
 * \ingroup alignment_configuration
 * \tparam end_gaps_t The type of the end-gaps. Must be a specialisation of seqan3::end_gaps.
 *
 * \details
 *
 * This configuration element configures the aligned ends to further refine the global alignment algorithm.
 * Particularly, the ends of the alignment can be penalised with gap costs or not. For example, the semi-global
 * alignment does not penalise the leading and trailing gaps of one sequence while it does for the other sequence.
 *
 * The class is instantiated with an object of seqan3::end_gaps. The user can configure each of the
 * gap specifier separately allowing for maximal flexibility when configuring the alignment algorithm.
 * However, there are also predefined \ref predefined_end_gap_configurations "configurations" which should be
 * preferred whenever possible.
 *
 * If this configuration element is not specified for the alignment algorithm, it will automatically default to
 * \ref seqan3::end_gaps::free_ends_none "free_ends_none" which computes a global alignment.
 *
 * ### Example
 *
 * \include snippet/alignment/configuration/align_cfg_aligned_ends.cpp
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
//!\brief Deduces the end-gaps object type from the constructor argument.
template <typename end_gaps_t>
aligned_ends(end_gaps_t) -> aligned_ends<std::remove_reference_t<end_gaps_t>>;
//!\}

} // namespace seqan3::align_cfg
