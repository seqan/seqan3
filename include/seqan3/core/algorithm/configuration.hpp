// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::configuration and utility functions.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <tuple>

#include <meta/meta.hpp>

#include <seqan3/core/algorithm/concept.hpp>
#include <seqan3/core/algorithm/configuration_utility.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>
#include <seqan3/core/tuple_utility.hpp>
#include <seqan3/core/type_list/type_list.hpp>
#include <seqan3/core/type_list/traits.hpp>
#include <seqan3/core/type_traits/template_inspection.hpp>
#include <seqan3/std/concepts>

namespace seqan3
{

// ----------------------------------------------------------------------------
// configuration
// ----------------------------------------------------------------------------

/*!\brief Collection of elements to configure an algorithm.
 * \ingroup algorithm
 *
 * \tparam configs_t Template parameter pack containing all configuration elements; Must model
 *                   seqan3::detail::config_element_specialisation
 *
 * \details
 *
 * This class provides a unified interface to create and query such
 * configurations for a specific algorithm. It extends the standard tuple interface with some useful functions to modify
 * and query the user configurations.
 */
template <detail::config_element_specialisation ... configs_t>
class configuration : public std::tuple<configs_t...>
{
    //!\brief Friend declaration for other instances of the configuration.
    template <detail::config_element_specialisation ... _configs_t>
    friend class configuration;

public:
    //!\privatesection
    //!\brief A type alias for the base class.
    using base_type = std::tuple<configs_t...>;

    //!\publicsection
    /*!\name Constructor, destructor and assignment
     * \{
     */
    constexpr configuration()                                  = default; //!< Defaulted.
    constexpr configuration(configuration const &)             = default; //!< Defaulted.
    constexpr configuration(configuration &&)                  = default; //!< Defaulted.
    constexpr configuration & operator=(configuration const &) = default; //!< Defaulted.
    constexpr configuration & operator=(configuration &&)      = default; //!< Defaulted.
    ~configuration()                                           = default; //!< Defaulted.

    /*!\brief Constructs a configuration from a single configuration element.
     * \tparam config_element_t The configuration element to add; must model
     *                          seqan3::detail::config_element_specialisation.
     * \param[in] config_element The configuration element to construct the configuration from.
     */
    template <typename config_element_t>
    //!\cond
        requires (!std::same_as<std::remove_cvref_t<config_element_t>, configuration>) &&
                 detail::config_element_specialisation<std::remove_cvref_t<config_element_t>>
    //!\endcond
    constexpr configuration(config_element_t && config_element) :
        base_type{std::forward<config_element_t>(config_element)}
    {}
    //!\}

    /*!\name Capacity
     * \{
     */

    //!\brief Returns the number of contained config elements.
    constexpr size_t size() const noexcept
    {
        return std::tuple_size_v<base_type>;
    }

    /*!\name Observers
     * \{
     */

    /*!\brief Returns the stored configuration element if present otherwise the given alternative.
     * \tparam alternative_t The type of the configuration element that is queried.
     *
     * \param[in] alternative The alternative whose type is used to check for an existing configuration element.
     *
     * \details
     *
     * Uses the type `alternative_t` of the given alternative to check if such an configuration element was already
     * stored inside of the configuration. If no suitable candidate can be found the passed value `alternative` will
     * be returned. If `alternative_t` is a class template, then any specialisation of this alternative type will be
     * searched and returned if present.
     *
     * \returns The stored configuration element identified by `alternative_t` or the alternative if not present.
     *
     * ### Example
     *
     * \include test/snippet/core/algorithm/configuration_get_or.cpp
     *
     * ### Exception
     *
     * no-throw guarantee.
     *
     * ### Complexity
     *
     * Constant time.
     */
    template <typename alternative_t>
    constexpr decltype(auto) get_or(alternative_t && alternative) & noexcept
    {
        return get_or_impl(*this, alternative, std::forward<alternative_t>(alternative));
    }

    //!\overload
    template <typename alternative_t>
    constexpr decltype(auto) get_or(alternative_t && alternative) const & noexcept
    {
        return get_or_impl(*this, alternative, std::forward<alternative_t>(alternative));
    }

    //!\overload
    template <typename alternative_t>
    constexpr decltype(auto) get_or(alternative_t && alternative) && noexcept
    {
        return get_or_impl(std::move(*this), alternative, std::forward<alternative_t>(alternative));
    }

    //!\overload
    template <typename alternative_t>
    constexpr decltype(auto) get_or(alternative_t && alternative) const && noexcept
    {
        return get_or_impl(std::move(*this), alternative, std::forward<alternative_t>(alternative));
    }

    //!\brief Checks if the given type exists in the tuple.
    template <typename query_t>
    static constexpr bool exists() noexcept
    {
        return pack_traits::contains<query_t, configs_t...>;
    }
    //!\brief Checks if the given type exists in the tuple.
    template <template <typename ...> typename query_t>
    static constexpr bool exists() noexcept
    {
        return (pack_traits::find_if<detail::is_same_configuration_f<query_t>::template invoke, configs_t...> > -1);
    }
    //!\}

    /*!\name Modifiers
     * \{
     */
    /*!\brief Remove a config element from the configuration.
     * \tparam query_t The config element type to remove from the configuration.
     * \returns A new configuration object without the config element identified by `query_t`.
     */
    template <typename query_t>
    [[nodiscard]] constexpr auto remove() const
    //!\cond
#if !SEQAN3_WORKAROUND_GCC_95371
        requires (exists<query_t>())
#endif // !SEQAN3_WORKAROUND_GCC_95371
    //!\endcond
    {
        constexpr int index = pack_traits::find<query_t, configs_t...>;
        return remove_at<index>();
    }

    //!\overload
    template <template <typename ...> typename query_t>
    [[nodiscard]] constexpr auto remove() const
    //!\cond
#if !SEQAN3_WORKAROUND_GCC_95371
        requires (exists<query_t>())
#endif // !SEQAN3_WORKAROUND_GCC_95371
    //!\endcond
    {
        constexpr int index = pack_traits::find_if<detail::is_same_configuration_f<query_t>::template invoke,
                                                   configs_t...>;
        return remove_at<index>();
    }
    //!\}

    /*!\brief Returns a new configuration by appending the given configuration to the current one.
     *
     * \tparam other_configuration_t Another configuration type or configuration element type; each configuration
     *                               element must model seqan3::detail::config_element_pipeable_with each of the
     *                               configurations elements of the current configuration.
     *
     * \param[in] other_config The other configuration to append to the current one.
     *
     * \returns A new configuration containing the appended configuration elements.
     *
     * \details
     *
     * This function generates a new configuration object containing the appended configuration elements. The current
     * configuration will not be modified.
     */
    template <typename other_configuration_t>
    //!\cond
        requires (is_config_element_combineable_v<configs_t, std::remove_cvref_t<other_configuration_t>> && ...)
    //!\endcond
    constexpr auto append(other_configuration_t && other_config) const
    {
        if constexpr (detail::config_element_specialisation<std::remove_cvref_t<other_configuration_t>>)
        {
            return configuration<configs_t..., std::remove_cvref_t<other_configuration_t>>
            {
                std::tuple_cat(static_cast<base_type>(*this),
                               std::tuple{std::forward<other_configuration_t>(other_config)})
            };
        }
        else
        {
            using other_base_t = typename std::remove_cvref_t<other_configuration_t>::base_type;
            using other_base_maybe_const_t =
                std::conditional_t<std::is_const_v<std::remove_reference_t<other_configuration_t>>,
                                   std::add_const_t<other_base_t>,
                                   other_base_t>;
            using fwd_other_base_t = std::conditional_t<std::is_lvalue_reference_v<other_configuration_t>,
                                                        std::add_lvalue_reference_t<other_base_maybe_const_t>,
                                                        std::add_rvalue_reference_t<other_base_maybe_const_t>>;

            using other_configs_list_t = detail::transfer_template_args_onto_t<other_base_t, type_list>;
            using appended_configuration_t =
                    detail::transfer_template_args_onto_t<list_traits::concat<type_list<configs_t...>,
                                                                              other_configs_list_t>,
                                                          configuration>;

            return appended_configuration_t{std::tuple_cat(static_cast<base_type>(*this),
                                                           static_cast<fwd_other_base_t>(other_config))};
        }
    }

private:
    /*!\name Internal constructor
     * \{
     */
    //!\brief Constructs from std::tuple.
    template <typename ..._configs_t>
    explicit constexpr configuration(std::tuple<_configs_t...> const & cfg) : base_type{cfg}
    {}

    //!\brief Constructs from std::tuple.
    template <typename ..._configs_t>
    explicit constexpr configuration(std::tuple<_configs_t...> && cfg) : base_type{std::move(cfg)}
    {}
    //!\}

    /*!\name Modifiers
     * \brief Note that modifications return new configurations and do not modify `this`.
     * \{
     */

    /*!\brief Remove a config element from the configuration.
     * \tparam index The config element at `index` is removed from the config.
     * \returns A new configuration object without the config element at `index`.
     */
    template <int index>
    [[nodiscard]] constexpr auto remove_at() const
    {
        static_assert((index >= 0) && (index < sizeof...(configs_t)), "Index to remove from config is out of bounds.");

        auto [head, middle] = tuple_split<index>(static_cast<base_type>(*this));
        auto tail = tuple_pop_front(middle);

        using head_list_t = detail::transfer_template_args_onto_t<decltype(head), type_list>;
        using tail_list_t = detail::transfer_template_args_onto_t<decltype(tail), type_list>;
        using concat_list_t = list_traits::concat<head_list_t, tail_list_t>;
        using new_configuration_t = detail::transfer_template_args_onto_t<concat_list_t, configuration>;

        return new_configuration_t{std::tuple_cat(std::move(head), std::move(tail))};
    }
    //!\}

    /*!\brief Internal implementation of the get_or interace.
     *
     * \tparam this_t The type of this.
     * \tparam query_t The type of the configuration element to query.
     * \tparam alternative_t The type of the alternative.
     *
     * \param[in] me The perfectly forwarded instance of `*this`.
     * \param[in] query The queried configuration element [only the type is needed].
     * \param[in] alternative The alternative configuration element to return if the query_t is not present.
     *
     * \details
     *
     * Use the type `query_t` to check if such a configuration element is stored in `me`. If this is `true` then
     * the stored configuration element is returned using perfect forwarding. If this evaluates to `false` the
     * given alternative is returned. If `query_t` is a class template then it is checked if any
     * specialisation of this class template is stored.
     */
    template <typename this_t, typename query_t, typename alternative_t>
    static constexpr decltype(auto) get_or_impl(this_t && me,
                                                query_t const & SEQAN3_DOXYGEN_ONLY(query),
                                                alternative_t && alternative) noexcept
    {
        if constexpr (exists<query_t>())
        {
            return get<query_t>(std::forward<this_t>(me));
        }
        else
        {
            using ret_type = remove_rvalue_reference_t<decltype(alternative)>;
            return static_cast<ret_type>(alternative);
        }
    }

    //!\overload
    template <typename this_t,
              template <typename ...> typename query_template_t, typename ...parameters_t,
              typename alternative_t>
    static constexpr decltype(auto) get_or_impl(this_t && me,
                                                query_template_t<parameters_t...> const &,
                                                alternative_t && alternative) noexcept
    {
        if constexpr (exists<query_template_t>())
        {
            return get<query_template_t>(std::forward<this_t>(me));
        }
        else
        {
            using ret_type = remove_rvalue_reference_t<decltype(alternative)>;
            return static_cast<ret_type>(alternative);
        }
    }
};

/*!\brief Combines configuration and configuration elements to create a new configuration object.
 * \relates seqan3::configuration
 *
 * \tparam lhs_config_t The type of the left hand side operand; the seqan3::is_config_element_combineable_v variable
 *                      template must evaluate to true for both operand types.
 * \tparam rhs_config_t The type of the right hand side operand; the seqan3::is_config_element_combineable_v variable
 *                      template must evaluate to true for both operand types.
 *
 * \param[in] lhs The left hand operand.
 * \param[in] rhs The right hand operand.
 *
 * \returns A new seqan3::configuration containing `lhs` and `rhs`.
 *
 * \details
 *
 * The the two operands can be either a seqan3::configuration object or a
 * \ref seqan3::detail::config_element_specialisation "configuration element". The right hand side operand is then
 * appended to the left hand side operand by creating a new seqan3::configuration object. Neither `lhs` nor `rhs` will
 * be modified.
 */
template <typename lhs_config_t, typename rhs_config_t>
//!\cond
    requires (is_config_element_combineable_v<std::remove_cvref_t<lhs_config_t>, std::remove_cvref_t<rhs_config_t>>)
//!\endcond
constexpr auto operator|(lhs_config_t && lhs, rhs_config_t && rhs)
{
    using pure_lhs_t = std::remove_cvref_t<lhs_config_t>;
    using pure_rhs_t = std::remove_cvref_t<rhs_config_t>;

    constexpr bool lhs_is_config_element = detail::config_element_specialisation<pure_lhs_t>;
    constexpr bool rhs_is_config_element = detail::config_element_specialisation<pure_rhs_t>;

    if constexpr (lhs_is_config_element && rhs_is_config_element)
        return configuration{std::forward<lhs_config_t>(lhs)}.append(std::forward<rhs_config_t>(rhs));
    else if constexpr (lhs_is_config_element && !rhs_is_config_element)
        return configuration{std::forward<lhs_config_t>(lhs)}.append(std::forward<rhs_config_t>(rhs));
    else if constexpr (!lhs_is_config_element && rhs_is_config_element)
        return std::forward<lhs_config_t>(lhs).append(std::forward<rhs_config_t>(rhs));
    else
        return lhs.append(rhs);
}

/*!\name Type deduction guides
 * \{
 */

/*!\brief Deduces the correct configuration element type from the passed seqan3::pipeable_config_element.
 * \relates seqan3::configuration
 */
template <typename config_t>
//!\cond
    requires detail::config_element_specialisation<std::remove_cvref_t<config_t>>
//!\endcond
configuration(config_t &&) -> configuration<std::remove_cvref_t<config_t>>;
//!\}

/*!\name Tuple interface
 * \{
 */

/*!\brief Returns the stored element.
 * \ingroup algorithm
 * \relates seqan3::configuration
 *
 * \tparam    query_t A template template.
 * \param[in] config  The configuration to get the element for.
 *
 * \details
 *
 * Extends the position-based and type based `get` interface for the configuration type, with a version
 * that also accepts template-template types (types that are itself templates), such that the exact template definition
 * must not be known.
 *
 * ### Example
 *
 * The following snippet demonstrates the various versions of get that can be used.
 *
 * \include test/snippet/core/algorithm/configuration_get.cpp
 *
 * ### Exception
 *
 * no-throw guarantee.
 *
 * ### Complexity
 *
 * Constant time.
 */
template <template <typename ...> class query_t, typename ...configs_t>
constexpr auto & get(configuration<configs_t...> & config) noexcept
{
    constexpr auto pos = pack_traits::find_if<detail::is_same_configuration_f<query_t>::template invoke, configs_t...>;
    static_assert(pos > -1, "Access error: The requested type is not contained.");

    return get<pos>(config);
}

//!\overload
template <template <typename ...> class query_t, typename ...configs_t>
constexpr auto const & get(configuration<configs_t...> const & config) noexcept
{
    constexpr auto pos = pack_traits::find_if<detail::is_same_configuration_f<query_t>::template invoke, configs_t...>;
    static_assert(pos > -1, "Access error: The requested type is not contained.");

    return get<pos>(config);
}

//!\overload
template <template <typename ...> class query_t, typename ...configs_t>
constexpr auto && get(configuration<configs_t...> && config) noexcept
{
    constexpr auto pos = pack_traits::find_if<detail::is_same_configuration_f<query_t>::template invoke, configs_t...>;
    static_assert(pos > -1, "Access error: The requested type is not contained.");

    return get<pos>(std::move(config));
}

//!\overload
template <template <typename ...> class query_t, typename ...configs_t>
constexpr auto const && get(configuration<configs_t...> const && config) noexcept
{
    constexpr auto pos = pack_traits::find_if<detail::is_same_configuration_f<query_t>::template invoke, configs_t...>;
    static_assert(pos > -1, "Access error: The requested type is not contained.");

    // TODO: change after GCC-7 bug with const && version of get in std::tuple is fixed.
    // return get<pos>(std::move(config));
    return std::move(get<pos>(config));
}
//!\}

} // namespace seqan3::detail

namespace std
{
//!\cond DEV

/*!\brief Returns the number of elements stored in seqan3::detail::configuration.
 * \implements seqan3::unary_type_trait
 * \see std::tuple_size_v
 * \ingroup algorithm
 */
template <seqan3::detail::config_element_specialisation ... configs_t>
struct tuple_size<seqan3::configuration<configs_t...>>
{
    //!\brief The number of elements.
    static constexpr size_t value = std::tuple_size_v<typename seqan3::configuration<configs_t...>::base_type>;
};

/*!\brief Returns the type of the element at the specified position within seqan3::configuration.
 * \implements seqan3::transformation_trait
 * \see [std::tuple_element](https://en.cppreference.com/w/cpp/utility/tuple/tuple_element)
 * \ingroup algorithm
 */
template <size_t pos, seqan3::detail::config_element_specialisation ... configs_t>
struct tuple_element<pos, seqan3::configuration<configs_t...>>
{
    //!\brief The type of the config at position `pos`
    using type = std::tuple_element_t<pos, typename seqan3::configuration<configs_t...>::base_type>;
};
//!\endcond
}  //namespace std
