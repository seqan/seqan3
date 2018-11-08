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
 * \brief Provides seqan3::detail::configuration and utility functions.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <tuple>

#include <meta/meta.hpp>

#include <seqan3/core/algorithm/concept.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/core/metafunction/template_inspection.hpp>
#include <seqan3/core/tuple_utility.hpp>
#include <seqan3/core/type_list.hpp>
#include <seqan3/std/concepts>

namespace seqan3::detail
{

//!\cond
// Forward declaration.
template <typename derived_t, typename ... args_t>
struct configuration_fn_proxy;

template <typename derived_fn_t>
class configuration_fn_base;
//!\endcond

// ----------------------------------------------------------------------------
// configuration
// ----------------------------------------------------------------------------

/*!\brief Collection of configurations objects used to specify the runtime behavior of algorithms.
 * \ingroup algorithm
 *
 * \tparam configs_t Template parameter pack containing all the inherited configs. Each must satisfy the
 *                   seqan3::detail::config_element_concept
 *
 * \details
 *
 * This class provides a unified interface and additional helper functions to create and query configurations
 * for a specific algorithm. Certain bioinformatics algorithms, e.g. alignment or find interfaces, can contain
 * many different configurations and policies that alter the execution of the algorithm. These configurations can be
 * orthogonal or might be mutually exclusive. Using this configuration the interface for the user becomes much more
 * easier, and incompatible configurations can be checked at compile time.
 *
 * ### Usage
 *
 * The entire configuration system is designed to work completely in the background of any algorithm. The implementor
 * of an algorithm only needs to provide the configurations for the algorithm.
 * In general the type of a configuration is static (see seqan3::detail::config_element_base) within the context of an
 * algorithm to select the correct code branches based on a policy-driven design, albeit the stored value might still be
 * a runtime parameter.
 * However, in some cases a specific configuration might be known only at runtime, and needs to be converted to a static
 * type for the algorithm in use. To enable a transparent conversion from the runtime parameters to a static type
 * configuration one can use deferred configs (seqan3::detail::deferred_config_element_base), which are invocable
 * configurations, that store the runtime parameter and on invocation translate this runtime parameter to a static
 * configuration. The implementor of an algorithm can achieve this by using the function
 * seqan3::detail::apply_deferred_configs, which iterates through all configurations and in case
 * it is a deferred configuration it will invoke the translation function and continue with the modified configuration,
 * which contains now the static type for the specific configuration.
 *
 * ### Combining Configurations
 *
 * To enable simple extension of configurations the configuration supports a logical-or interface for the different
 * configurations and configuration adaptors. A configuration adaptor is a functor, which provides the logical-or
 * interface. Adaptors need to implement the abstract base seqan3::detail::configuration_fn_base,
 * to add all the necessary interfaces.
 * Using this, a configuration class can be easily constructed by chaining together different config elements.
 * Consider the following example, where we assume that the config element `bar` and the configuration adaptor
 * `with_foo` are already given:
 *
 * \snippet test/snippet/core/algorithm/configuration.cpp combine
 * The `with_foo` adaptor captures the `configuration<bar>` object and returns a new configuration with the added
 * `foo`-config.

 * ### Access the data
 *
 * The configuration inherits from a std::tuple and exposes a tuple like interface. To access a specific element one
 * can use `std::get` to query the config element at the specified position or for the respective type
 * within the configuration.
 * Considering the example from above one can get the value of the bar-config as:
 *
 * \snippet test/snippet/core/algorithm/configuration.cpp access
 * Note, that the type based get should be used with caution. In the typical usage scenario, the user of the
 * configuration works with config adaptor types. Only these adaptors know what exact type was added to the
 * configuration.
 * To allow algorithm implementors to still access data for specific properties independent of their type and position,
 *  every algorithm will expose a enum based get interface, that allows to access the element associated with the
 * enum identifier. The details can be found in the respective algorithm description.
 */
template <config_element_concept ... configs_t>
class configuration : public std::tuple<configs_t...>
{
    //!\brief Friend declaration for other instances of the configuration.
    template <config_element_concept ... _configs_t>
    friend class configuration;

public:
    /*!\name Member types
     * \{
     */
    //!\brief A type alias for the base class.
    using base_type = std::tuple<configs_t...>;
    //!\}

    /*!\name Constructor, destructor and assignment
     * \{
     */
    constexpr configuration()                                  = default;
    constexpr configuration(configuration const &)             = default;
    constexpr configuration(configuration &&)                  = default;
    constexpr configuration & operator=(configuration const &) = default;
    constexpr configuration & operator=(configuration &&)      = default;
    ~configuration()                                           = default;

    //!\brief Constructs from std::tuple.
    constexpr configuration(std::tuple<configs_t...> const & cfg) : base_type{cfg}
    {}

    //!\brief Constructs from std::tuple.
    constexpr configuration(std::tuple<configs_t...> && cfg) : base_type{std::move(cfg)}
    {}

    /*!\brief Constructs a configuration from a variable of type seqan3::detail::configuration_fn_proxy or
     *        seqan3::detail::configuration_fn_base
     * \tparam cfg_fn_t Either seqan3::detail::configuration_fn_proxy or seqan3::detail::configuration_fn_base.
     * \param cfg_fn    A variable of type `cfg_fn_t`.
     *
     * \details
     *
     * In some scenarios it might be possible, that only a single configuration element is used. In this case, the
     * configuration must be constructed from a variable of type seqan3::detail::configuration_fn_proxy
     * or seqan3::detail::configuration_fn_base. In order to construct a configuration object from it, this constructor
     * invokes the functor with an empty configuration.
     *
     * ### Example
     *
     * \snippet test/snippet/core/algorithm/configuration.cpp constructor
     * In the above example, we assume, that `my_cfg` is a variable of type seqan3::detail::configuration_fn_base.
     */
    template <typename cfg_fn_t>
    //!\cond
        requires is_type_specialisation_of_v<std::remove_reference_t<cfg_fn_t>, configuration_fn_proxy> ||
                 std::is_base_of_v<configuration_fn_base<remove_cvref_t<cfg_fn_t>>, remove_cvref_t<cfg_fn_t>>
    //!\endcond
    constexpr configuration(cfg_fn_t && cfg_fn) :
        configuration{std::invoke(std::forward<cfg_fn_t>(cfg_fn), configuration<>{})}
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

    /*!\name Modifiers
     * \brief Note that modifications return new configurations and do not modify `this`.
     * \{
     */

    /*!\brief Adds a new config element to the beginning of the configuration.
     *
     * \param[in] cfg_element The element to add.
     *
     * \returns A new seqan3::detail::configuration containing the added element.
     *
     * \details
     *
     * Creates a new seqan3::detail::configuration from `this` and adds the passed config element.
     * Note, that `this` is not modified by this operation.
     *
     * Complexity
     *
     * Linear in the number of elements.
     */
    template <config_element_concept config_element_t>
    constexpr auto push_front(config_element_t && cfg_element) const &
    {
        return detail::configuration{std::tuple_cat(std::tuple{std::forward<config_element_t>(cfg_element)},
                                                    static_cast<base_type>(*this))};
    }

    //!\copydoc push_front
    template <config_element_concept config_element_t>
    constexpr auto push_front(config_element_t && cfg_element) &&
    {
        return detail::configuration{std::tuple_cat(std::tuple{std::forward<config_element_t>(cfg_element)},
                                                    std::move(static_cast<base_type>(*this)))};
    }

    /*!\brief Replaces the old config element with the new one.
     *
     * \param[in] old_element The element to replace.
     * \param[in] new_element The new element to insert.
     *
     * \returns A new seqan3::detail::configuration containing the replaced element.
     *
     * \details
     *
     * Splits the `this` at the position of the `old_config_element_t` and replaces it with `new_element` at
     * the same position and constructs a new seqan3::detail::configuration from it. Note, that `this` is not
     * modified by this operation.
     *
     * Complexity
     *
     * Linear in the number of elements.
     */
    template <config_element_concept old_config_element_t,
              config_element_concept new_config_element_t>
    constexpr auto replace_with(old_config_element_t const & SEQAN3_DOXYGEN_ONLY(old_element),
                                new_config_element_t && new_element) const &
    {
        static_assert(std::tuple_size_v<base_type> > 0, "The configuration cannot be empty.");
        static_assert(meta::find_index<transfer_template_args_onto_t<base_type, type_list>,
                                       old_config_element_t>::value != meta::npos::value,
                      "The element to be replaced is not contained in the passed configuration.");

        auto && [prefix, remainder] = seqan3::tuple_split<old_config_element_t>(static_cast<base_type>(*this));

        return detail::configuration{std::tuple_cat(std::move(prefix),
                                                    std::tuple{std::forward<new_config_element_t>(new_element)},
                                                    tuple_pop_front(std::move(remainder)))};
    }

    //!\copydoc replace_with
    template <config_element_concept old_config_element_t,
              config_element_concept new_config_element_t>
    constexpr auto replace_with(old_config_element_t const & SEQAN3_DOXYGEN_ONLY(old_element),
                                new_config_element_t && new_element) &&
    {
        static_assert(std::tuple_size_v<base_type> > 0, "The configuration cannot be empty.");
        static_assert(meta::find_index<transfer_template_args_onto_t<base_type, type_list>,
                                       old_config_element_t>::value != meta::npos::value,
                      "The element to be replaced is not contained in the passed configuration.");

        auto && [prefix, remainder] =
            seqan3::tuple_split<old_config_element_t>(std::move(static_cast<base_type>(*this)));

        return detail::configuration{std::tuple_cat(std::move(prefix),
                                                    std::tuple{std::forward<new_config_element_t>(new_element)},
                                                    tuple_pop_front(std::move(remainder)))};
    }
    //!\}
};

/*!\name Type deduction guides
 * \relates seqan3::detail::configuration
 * \{
 */

//!\brief Deduces the correct templates from a seqan3::detail::configuration_fn_proxy
template <typename cfg_fn_t>
//!\cond
    requires is_type_specialisation_of_v<std::remove_reference_t<cfg_fn_t>, configuration_fn_proxy> ||
             std::is_base_of_v<configuration_fn_base<remove_cvref_t<cfg_fn_t>>, remove_cvref_t<cfg_fn_t>>
//!\endcond
configuration(cfg_fn_t &&) ->
    configuration<meta::front<detail::tuple_type_list_t<typename std::invoke_result_t<std::remove_reference_t<cfg_fn_t>,
                                                                                      configuration<>>::base_type>>>;
//!\}

// ----------------------------------------------------------------------------
// Metafunction is_configuration_combinable_with
// ----------------------------------------------------------------------------

/*!\brief Variable template that checks if two adaptor types can be combined with the logical-or operator.
 * \ingroup algorithm
 * \tparam target_t The type of the left operand.
 * \tparam query_t The type of the right operand.
 *
 * \returns `true` if `target_t` and `query_t` are either a derived type of seqan3::detail::configuration_fn_base or a
 *           seqan3::detail::configuration_fn_proxy and `target_t` and `query_t` are not the same.
 *           Otherwise it returns `false`.
 *
 * \see seqan3::detail::is_configuration_combinable_with_v.
 */
template <typename target_t, typename query_t>
//!\cond
    requires !is_algorithm_configuration_v<remove_cvref_t<target_t>>
//!\endcond
struct is_configuration_combinable_with
{
    //!\brief The result of the test expression.
    static constexpr bool value =
        (std::is_base_of_v<configuration_fn_base<remove_cvref_t<target_t>>, remove_cvref_t<target_t>> ||
         is_type_specialisation_of_v<remove_cvref_t<target_t>, configuration_fn_proxy>) &&
        (is_type_specialisation_of_v<remove_cvref_t<query_t>, configuration_fn_proxy> ||
        (std::is_base_of_v<configuration_fn_base<remove_cvref_t<query_t>>, remove_cvref_t<query_t>> &&
            !std::is_same_v<remove_cvref_t<target_t>, remove_cvref_t<query_t>>));
};

//!\brief Helper variable template for seqan3::detail::is_configuration_combinable_with.
//!\ingroup algorithm
template <typename target_t, typename query_t>
//!\cond
    requires !is_algorithm_configuration_v<remove_cvref_t<target_t>>
//!\endcond
inline constexpr bool is_configuration_combinable_with_v = is_configuration_combinable_with<target_t, query_t>::value;

// ----------------------------------------------------------------------------
// configuration_fn_base
// ----------------------------------------------------------------------------

/*!\brief An abstract crtp-base class to add pipeable concept for configuration functors in combination with
          seqan3::detail::configuration.
 * \ingroup algorithm
 * \tparam derived_t The derived type to be extended with functionality from this base class.
 *
 * \details
 *
 * This abstract crtp-base class provides the configuration adaptor interface and the pipeable interface
 * for configurations and configuration elements.
 */
template <typename derived_fn_t>
class configuration_fn_base
{
protected:

    /*!\name Constructor, destructor and assignment
     * \brief Abstract base class has no public constructor.
     * \{
     */
    configuration_fn_base()                                          = default;
    configuration_fn_base(configuration_fn_base const &)             = default;
    configuration_fn_base(configuration_fn_base &&)                  = default;
    configuration_fn_base & operator=(configuration_fn_base const &) = default;
    configuration_fn_base & operator=(configuration_fn_base &&)      = default;
    ~configuration_fn_base()                                         = default;
    //!}

public:

    /*!\brief Invokes the configuration specific functor to extend the seqan3::detail::configuration with the associated
     *        config.
     * \tparam configuration_t The configuration to be extended with the new configuration type. Must be of type
     *                        seqan3::detail::configuration.
     * \tparam args_t A template parameter pack with the arguments applied to the configuration functor call.
     *
     * \param configuration The configuration to be extended.
     * \param args         The argument pack to be passed to the functor.
     * \returns A new seqan3::detail::configuration with the extended configuration.
     */
    template <typename configuration_t,
              typename ... args_t>
    constexpr auto operator()(configuration_t && configuration,
                              args_t && ... args) const
    //!\cond
        requires is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
    //!\endcond
    {
        return static_cast<derived_fn_t const &>(*this).invoke(std::forward<configuration_t>(configuration),
                                                               std::forward<args_t>(args)...);
    }

    /*!\brief Creates a proxy caching the arguments that should be applied when invoking the configuration specific
     *        configuration.
     * \tparam args_t A template parameter pack with the arguments applied to the configuration functor call.
     *
     * \param args A pack containing the arguments.
     * \returns A proxy type used to defer configuration invocation.
     */
    template <typename ... args_t>
    constexpr configuration_fn_proxy<derived_fn_t, args_t...> operator()(args_t && ... args) const
    {
        return {std::forward<args_t>(args)...};
    }
};

/*!\name Pipe Interface
 * \relates seqan3::detail::configuration_fn_base
 * \{
 */

/*!\brief Combines a seqan3::detail::configuration with a configuration adaptor.
 * \tparam configuration_t The configuration type to be extended. Must be of type seqan3::detail::configuration.
 * \tparam fn_t            The type of the right operand. Must be of type seqan3::detail::configuration_fn_proxy or
 *                         seqan3::detail::configuration_fn_base
 *
 * \param[in] cfg The configuration to be extended.
 * \param[in] fn  The right operand.
 * \returns The result of invoking the configuration adaptor with the passed `cfg` object.
 */
template <typename configuration_t,
          typename fn_t>
constexpr auto operator|(configuration_t && cfg,
                         fn_t && fn)
//!\cond
    requires is_algorithm_configuration_v<remove_cvref_t<configuration_t>> &&
             (is_type_specialisation_of_v<remove_cvref_t<fn_t>, configuration_fn_proxy> ||
              std::is_base_of_v<configuration_fn_base<remove_cvref_t<fn_t>>, remove_cvref_t<fn_t>>)
//!\endcond
{
    return fn(std::forward<configuration_t>(cfg));
}

/*!\brief Combines a configuration adaptor with another configuration adaptor.
 * \tparam lhs_fn_t The type of the left operand. Must be of type seqan3::detail::configuration_fn_base or
 *                  seqan3::detail::configuration_fn_proxy.
 * \tparam rhs_fn_t The type of the right operand. Must be of type seqan3::detail::configuration_fn_base or
 *                  seqan3::detail::configuration_fn_proxy.
 *
 * \param[in] lhs_fn A configuration adaptor or a proxy there of.
 * \param[in] rhs_fn A configuration adaptor or a proxy there of.
 * \returns The result of invoking the right operand with the result of invoking the left operand.
 *
 * \details
 *
 * Allows any configuration adaptor or a proxy there of to be at the beginning of a configuration
 * declaration.
 *
 * ### Example
 *
 * \snippet test/snippet/core/algorithm/configuration.cpp combine_2
 */
template <typename lhs_fn_t,
          typename rhs_fn_t>
constexpr auto operator|(lhs_fn_t && lhs_fn,
                         rhs_fn_t && rhs_fn)
//!\cond
    requires is_configuration_combinable_with_v<lhs_fn_t, rhs_fn_t>
//!\endcond
{
    return rhs_fn(std::forward<lhs_fn_t>(lhs_fn)(configuration<>{}));
}
//!\}

// ----------------------------------------------------------------------------
// configuration_fn_proxy
// ----------------------------------------------------------------------------

/*!\brief A proxy class used to defer invocation of the actual functor.
 * \ingroup algorithm
 * \tparam derived_t The template parameter of seqan3::detail::configuration_fn_base.
 * \tparam args_t    Template parameter pack with intermediate arguments that should be applied on invocation.
 *
 * \details
 *
 * This class is a helper proxy class for some invocations of seqan3::detail::configuration_fn_base, i.e.
 * the associated functor is invoked with the arguments passed to the configuration element that should be
 * created. In this case, the functor returns a proxy, which caches the arguments to the configuration element.
 * This proxy can only be constructed from within seqan3::detail::configuration_fn_base class.
 * There are special pipe-operator overloads, that work in combination with this proxy implementation.
 */
template <typename derived_t, typename ... args_t>
struct configuration_fn_proxy
{
protected:

    //!\brief Friend declaration to allow only seqan3::detail::config_element_base access to the constructor.
    template <typename t>
    friend class configuration_fn_base;

    //!\brief Constructs this proxy while forwarding arguments.
    constexpr configuration_fn_proxy(args_t && ... args) : args_cache{std::forward<args_t>(args)...}
    {}

    /*!\name Auxiliary functions
     * \brief Expand the cached elements.
     * \{
     */
    template <typename configuration_t, std::size_t ... Is>
    constexpr auto explode(configuration_t && cfg, std::index_sequence<Is...> /*unused*/) &&
    {
        // Move out the elements from cache.
        return derived_t{}(std::forward<configuration_t>(cfg), std::move(std::get<Is>(args_cache))...);
    }

    template <typename configuration_t, std::size_t ... Is>
    constexpr auto explode(configuration_t && cfg, std::index_sequence<Is...> /*unused*/) const &
    {
        // Copy the elements from cache.
        return derived_t{}(std::forward<configuration_t>(cfg), std::get<Is>(args_cache)...);
    }
    //!\}

    //!\brief The cached data.
    std::tuple<args_t...> args_cache;

public:
    /*!\name Invocable interface
     * \{
     */
    template <typename configuration_t>
    constexpr auto operator()(configuration_t && cfg) &&
    {
        return explode(std::forward<configuration_t>(cfg), std::make_index_sequence<sizeof...(args_t)>{});
    }

    template <typename configuration_t>
    constexpr auto operator()(configuration_t && cfg) const &
    {
        return explode(std::forward<configuration_t>(cfg), std::make_index_sequence<sizeof...(args_t)>{});
    }
    //!\}
};

// ============================================================================
// Utility functions
// ============================================================================

// ----------------------------------------------------------------------------
// apply_deferred_configs
// ----------------------------------------------------------------------------

//!\cond
// TODO Document me!
template <std::size_t N, typename fn_t, typename config_t>
constexpr auto apply_deferred_configs(fn_t && fn,
                                      config_t && config)
    requires is_algorithm_configuration_v<remove_cvref_t<config_t>>
{
    if constexpr (N == 0)
    {
        return fn(std::forward<config_t>(config));
    }
    else
    {
        using config_as_list = detail::transfer_template_args_onto_t<remove_cvref_t<config_t>, type_list>;
        using current_config_t = meta::at_c<config_as_list, meta::size<config_as_list>::value - N>;

        auto delegate = [&fn](auto && new_config)
        {
            return apply_deferred_configs<N-1>(fn, std::forward<decltype(new_config)>(new_config));
        };

        if constexpr (std::is_invocable_v<current_config_t, decltype(delegate), remove_cvref_t<config_t>>)
        {
            return current_config_t{}(delegate, std::forward<config_t>(config));
        }
        else
        {  // Perfect forwarding of the underlying config.
            return apply_deferred_configs<N - 1>(fn, std::forward<config_t>(config));
        }
    }
}

template <typename fn_t, typename config_t>
constexpr auto apply_deferred_configs(fn_t & fn,
                                      config_t && config)
    requires is_algorithm_configuration_v<remove_cvref_t<config_t>> &&
             std::Invocable<std::remove_reference_t<fn_t>, std::remove_reference_t<config_t>>
{
    using type_list_t = detail::tuple_type_list_t<typename std::remove_reference_t<config_t>::base_type>;
    return apply_deferred_configs<meta::size<type_list_t>::value>(std::forward<fn_t>(fn),
                                                                  std::forward<config_t>(config));
}
//!\endcond
} // namespace seqan3::detail

namespace std
{
//!\cond DEV

/*!\brief Returns the number of elements stored in seqan3::detail::configuration.
 * \ingroup algorithm
 */
template <seqan3::detail::config_element_concept ... configs_t>
struct tuple_size<seqan3::detail::configuration<configs_t...>>
{
    //!\brief The number of elements.
    static constexpr size_t value = std::tuple_size_v<typename seqan3::detail::configuration<configs_t...>::base_type>;
};

/*!\brief Returns the type of the element at the specified position within seqan3::detail::configuration.
 * \ingroup algorithm
 */
template <size_t pos, seqan3::detail::config_element_concept ... configs_t>
struct tuple_element<pos, seqan3::detail::configuration<configs_t...>>
{
    //!\brief The type of the config at position `pos`
    using type = std::tuple_element_t<pos, typename seqan3::detail::configuration<configs_t...>::base_type>;
};
//!\endcond
}  //namespace std
