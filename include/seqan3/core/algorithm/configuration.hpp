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
#include <seqan3/core/type_list.hpp>
#include <seqan3/std/concept/callable.hpp>
#include <seqan3/std/concept/core_language.hpp>

// Forward declaration for providing get interface for top-level namespace.
namespace seqan3
{
//!\cond
template <size_t elem_no, typename ... _configs_t>
constexpr auto & get(detail::configuration<_configs_t...> &) noexcept;

template <size_t elem_no, typename ... _configs_t>
constexpr auto const & get(detail::configuration<_configs_t...> const &) noexcept;

template <size_t elem_no, typename ... _configs_t>
constexpr auto && get(detail::configuration<_configs_t...> &&) noexcept;

template <size_t elem_no, typename ... _configs_t>
constexpr auto const && get(detail::configuration<_configs_t...> const &&) noexcept;

template <typename target_t, typename ... _configs_t>
constexpr auto & get(detail::configuration<_configs_t...> &) noexcept;

template <typename target_t, typename ... _configs_t>
constexpr auto const & get(detail::configuration<_configs_t...> const &) noexcept;

template <typename target_t, typename ... _configs_t>
constexpr auto && get(detail::configuration<_configs_t...> &&) noexcept;

template <typename target_t, typename ... _configs_t>
constexpr auto const && get(detail::configuration<_configs_t...> const &&) noexcept;
//!\endcond
}

namespace seqan3::detail
{

//!\cond
// Forward declaration.
template <typename ... args_t>
struct configuration_fn_proxy;

template <typename derived_fn_t>
class configuration_fn_base;
//!\endcond

// ----------------------------------------------------------------------------
// configuration
// ----------------------------------------------------------------------------

/*!\brief Collection of configurations objects used to specify the runtime behavior of algorithms.
 * \ingroup core_algorithm
 *
 * \tparam configs_t Template parameter pack containing all the inherited configs. Must satisfy the
 *                   seqan3::detail::config_element_concept
 *
 * \details
 *
 * This class provides a unified interface and additional helper functions to create and query configurations
 * for a specific algorithm. Certain bioinformatics algorithms, e.g. alignment or find interfaces, support a various set
 * of different configurations and policies that alter the execution of the algorithm. These configurations can be
 * orthogonal or might be mutual exclusive. Using this configuration the interface for the user becomes much more easier,
 * and incompatible configurations can be checked at compile time.
 *
 * \attention
 *
 * This class uses multiple inheritance to inherit from the actual configuration types. This also means, that the same
 * configuration type cannot occur more than once in the type definition of the configuration class.
 *
 * ### Usage
 *
 * The entire configuration system is designed to work completely in the background of any algorithm. The implementor
 * of an algorithm only needs to provide the configurations for the algorithm.
 * In general the type of a configuration is static (see seqan3::detail::config_element_base) within the context of an algorithm
 * to select the correct code branches based on a policy-driven design, albeit the stored state might still be a runtime
 * parameter.
 * However, in some cases a specific configuration might be known first at runtime but needs to be converted to a static
 * type for the algorithm in use. To enable a transparent conversion from the runtime parameters to a static type
 * configuration one can use deferred configs (seqan3::detail::deferred_config_element_base), which are invocable
 * configurations, that store the runtime parameter and on invocation translate this runtime parameter to a static
 * configuration. The implementor of an algorithm can achieve this by using the function
 * seqan3::detail::apply_deferred_configs, which iterates through all configurations and in case
 * it is a deferred configuration it will invoke the translation function and continue with the modified configuration,
 * which contains now the static type for the specific configuration.
 *
 * ### Pipe Notation
 *
 * To enable simple extension of configurations the configuration provides a generic pipe interface for the different
 * configurations. Thus, a config class can be easily constructed by chaining together different properties.
 * Consider the following example, where we assume that the config `bar` and the config functor `with_foo` are already
 * given:
 *
 * ```cpp
 * auto my_cfg = configuration<bar>{} | with_foo;  // my_cfg is now of type configuration<foo, bar>
 * ```
 * The `with_foo` functor captures the `configuration<bar>` object and returns a new configuration with the added
 * `foo`-config.

 * ### Tuple interface
 *
 * The configuration exposes a tuple interface. The type_trait functions std::tuple_size and std::tuple_element are
 * overloaded for this class. To access a specific element one can either use `seqan3::get` or `std::get` to query
 * the state of the configuration. Considering the example from above one can get the state of
 * the bar-config either as:
 *
 * ```cpp
 * auto bar_state = get<1>(my_cfg);
 * ```
 *
 * or as:
 *
 * ```cpp
 * auto bar_state = get<bar>(my_cfg);
 * ```
 *
 * In addition every algorithm specific config implementation can overload the get-interface with a enum based
 * get-interface, such that the correct config is extracted based on some type-independent enum value.
 * This is necessary for the algorithm implementation, as the algorithm neither knows the order nor the exact type of
 * all configs (the config type could be a template class).
 */
template <config_element_concept ... configs_t>
class configuration
{
    //!\brief Friend declaration for other instances of the configuration.
    template <typename ... _configs_t>
    friend class configuration;

    //!\brief Helper function for constructing a configuration from the predecessor prior it's modification.
    template <typename ... old_configs_t,
              size_t ... Is>
    constexpr void transfer_config_elements(configuration<old_configs_t...> const & other,
                                            std::index_sequence<Is...> const & /*idx*/)
    {
        ((std::get<Is + 1>(this->config_elements) =
            std::get<std::tuple_element_t<Is + 1, std::tuple<configs_t...>>>(other.config_elements)), ...);
    }

    //!\copydoc configuration::transfer_config_elements
    template <typename ... old_configs_t,
              size_t ... Is>
    constexpr void transfer_config_elements(configuration<old_configs_t...> && other,
                                            std::index_sequence<Is...> const & /*idx*/)
    {
        ((std::get<Is + 1>(this->config_elements) =
            std::get<std::tuple_element_t<Is + 1, std::tuple<configs_t...>>>(std::move(other.config_elements))), ...);
    }

public:
    /*!\name Member types
     * \{
     */
    //!\brief Exposes the configurations as a seqan3::type_list to allow meta operations on this list.
    using type_list_type = type_list<configs_t...>;
    //!\}

    /*!\name Constructor, destructor and assignment
     * \brief Defaulted all standard constructor.
     * \{
     */
    configuration()                                  = default;
    configuration(configuration const &)             = default;
    configuration(configuration &&)                  = default;
    configuration & operator=(configuration const &) = default;
    configuration & operator=(configuration &&)      = default;
    ~configuration()                                 = default;

    /*!\brief Constructs a configuration from the configuration prior to the modification with the new config.
     * \tparam other_configs_t The other configs.
     * \param cfg The configuration immediately prior to it's extension by another configuration.
     *
     * \details
     *
     * Note, that this constructor is only enabled for types that are in direct relation to this class, namely for
     * types for which the following invariant holds: `configuration<configs...>` corresponds to
     * `configuration<head, configs...>`, with `decltype(*this) = configuration<head, configs...>`.
     */
    template <typename ... other_configs_t>
    constexpr configuration(configuration<other_configs_t ...> const & cfg) : configuration{}
    {
        if constexpr (sizeof...(configs_t) > 0)
            transfer_config_elements(cfg, std::make_index_sequence<sizeof...(configs_t) - 1>{});
    }

    //!\copydoc configuration::configuration(configuration<other_configs_t ...> const & cfg)
    template <typename ... other_configs_t>
    constexpr configuration(configuration<other_configs_t ...> && cfg) : configuration{}
    {
        if constexpr (sizeof...(configs_t) > 0)
            transfer_config_elements(std::move(cfg), std::make_index_sequence<sizeof...(configs_t) - 1>{});
    }

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
     * ```cpp
     * detail::configuration cfg = my_cfg(1, 2, 3);
     * ```
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

    //!\cond
    template <size_t elem_no, typename ... _configs_t>
        requires (std::is_same_v<configs_t, _configs_t> && ...)
    friend constexpr auto & seqan3::get(configuration<_configs_t...> & cfg) noexcept;

    //!\copydoc seqan3::get()
    template <size_t elem_no, typename ... _configs_t>
        requires (std::is_same_v<configs_t, _configs_t> && ...)
    friend constexpr auto const & seqan3::get(configuration<_configs_t...> const &) noexcept;

    //!\copydoc seqan3::get<elem_no>(seqan3::detail::configuration)
    template <size_t elem_no, typename ... _configs_t>
        requires (std::is_same_v<configs_t, _configs_t> && ...)
    friend constexpr auto && seqan3::get(configuration<_configs_t...> &&) noexcept;

    //!\copydoc seqan3::get<elem_no>(seqan3::detail::configuration)
    template <size_t elem_no, typename ... _configs_t>
        requires (std::is_same_v<configs_t, _configs_t> && ...)
    friend constexpr auto const && seqan3::get(configuration<_configs_t...> const &&) noexcept;

    template <typename target_t, typename ... _configs_t>
        requires (std::is_same_v<configs_t, _configs_t> && ...)
    friend constexpr auto & seqan3::get(configuration<_configs_t...> &) noexcept;

    //!\copydoc configuration::get()
    template <typename target_t, typename ... _configs_t>
        requires (std::is_same_v<configs_t, _configs_t> && ...)
    friend constexpr auto const & seqan3::get(configuration<_configs_t...> const &) noexcept;

    //!\copydoc configuration::get()
    template <typename target_t, typename ... _configs_t>
        requires (std::is_same_v<configs_t, _configs_t> && ...)
    friend constexpr auto && seqan3::get(configuration<_configs_t...> &&) noexcept;

    //!\copydoc configuration::get()
    template <typename target_t, typename ... _configs_t>
        requires (std::is_same_v<configs_t, _configs_t> && ...)
    friend constexpr auto const && seqan3::get(configuration<_configs_t...> const &&) noexcept;
    //!\endcond

private:

    //!\brief The cached config elements.
    std::tuple<configs_t...> config_elements;
};

/*!\name Deduction guide
 * \relates seqan3::detail::configuration
 * \{
 */

//!\brief Deduces the correct templates from a seqan3::detail::configuration_fn_proxy
template <typename cfg_proxy_t>
//!\cond
    requires is_type_specialisation_of_v<std::remove_reference_t<cfg_proxy_t>, configuration_fn_proxy>
//!\endcond
configuration(cfg_proxy_t &&) ->
    configuration<meta::front<typename std::invoke_result_t<std::remove_reference_t<cfg_proxy_t>,
                                                            configuration<>>::type_list_type>>;

//!\brief Deduces the correct templates from a seqan3::detail::configuration_fn_base
template <typename cfg_fn_t>
//!\cond
    requires std::is_base_of_v<configuration_fn_base<remove_cvref_t<cfg_fn_t>>, remove_cvref_t<cfg_fn_t>>
//!\endcond
configuration(cfg_fn_t &&) ->
    configuration<meta::front<typename std::invoke_result_t<std::remove_reference_t<cfg_fn_t>,
                                                                configuration<>>::type_list_type>>;
//!\}
// ----------------------------------------------------------------------------
// Metafunction replace_config_with
// ----------------------------------------------------------------------------

//!\cond
template <typename configuration_t, detail::config_element_concept old_config_t, detail::config_element_concept new_config_t>
struct replace_config_with;
//!\endcond

/*!\brief Helper traits meta-function to replace one configuration type with another by pushing it to the front of
 *        the configurations for a specific seqan3::detail::configuration.
 * \ingroup core_algorithm
 * \relates seqan3::detail::configuration
 *
 * \tparam configuration_t The configuration type to be altered. Must be a template class.
 * \tparam old_config_t   The old_config_t to be removed.
 * \tparam new_config_t   The new_config_t to be added at the front of the configurations.
 *
 * \details
 *
 * This meta-function parses the list of config types contained in the given seqan3::detail::configuration and
 * removes the `old_config_t` from the list if it exists, and adds `new_config_t` at the front of
 * the type list. Subsequently a new `configuration_t` type is defined containing the new config type list.
 *
 * \see seqan3::detail::push_front_config
 * \see seqan3::detail::replace_config_with_t
 */
template <template <typename ...> class configuration_t,
          detail::config_element_concept old_config_t,
          detail::config_element_concept new_config_t,
          detail::config_element_concept ... configs_t>
struct replace_config_with<configuration_t<configs_t...>, old_config_t, new_config_t>
{
protected:

    //!\cond
    template <typename target_t>
    struct invoke_exclude
    {
        /*!\brief A template type which is used by meta::invoke.
         * \tparam query_t The type to query.
         * \returns `std::true_type` if `query_t` and `target_t` are not the same, `std::false_type` otherwise.
         */
        template <typename query_t>
        using invoke = std::negation<std::is_same<query_t, target_t>>;
    };

    //!\brief Helper typedef with the modified type_list.
    using _new_list = meta::push_front<meta::filter<seqan3::type_list<configs_t...>,
                                                    invoke_exclude<old_config_t>>,
                                       new_config_t>;
    //!\endcond
public:

    //!\brief Typedef with the modified configuration type.
    using type = detail::transfer_template_args_onto_t<_new_list, configuration_t>;
};

/*!\brief Shortcut typedef for seqan3::detail::replace_config_with.
 * \ingroup core_algorithm
 * \see seqan3::detail::replace_config_with
 */
template <typename configuration_t, typename old_config_t, typename new_config_t>
using replace_config_with_t = typename replace_config_with<configuration_t, old_config_t, new_config_t>::type;

//!\cond
template <typename configuration_t, detail::config_element_concept new_config_t>
struct push_front_config;
//!\endcond

// ----------------------------------------------------------------------------
// Metafunction push_front_config
// ----------------------------------------------------------------------------

/*!\brief Helper traits meta-function to replace one configuration type with another by pushing it to the front of
 *        the configurations for a specific seqan3::detail::configuration.
 * \ingroup core_algorithm
 *
 * \tparam configuration_t The configuration type to be altered. Must be a template class.
 * \tparam new_config_t   The new_config_t to be added at the front of the configurations.
 *
 * \details
 *
 * This meta-function adds `new_config_t` at the front of the configuration.
 *
 * \see seqan3::detail::push_front_config_t
 * \see seqan3::detail::replace_config_with
 */
template <template <typename ...> class configuration_t,
          detail::config_element_concept new_config_t,
          detail::config_element_concept ... configs_t>
struct push_front_config<configuration_t<configs_t...>, new_config_t>
{
private:
    //!\brief Helper typedef with the modified type_list.
    using _new_list = meta::push_front<seqan3::type_list<configs_t...>, new_config_t>;
public:

    //!\brief Typedef with the modified configuration type.
    using type = detail::transfer_template_args_onto_t<_new_list, configuration_t>;
};

/*!\brief Shortcut typedef for seqan3::detail::push_front_config.
 * \ingroup core_algorithm
 * \see seqan3::detail::push_front_config
 */
template <typename configuration_t, typename new_config_t>
using push_front_config_t = typename push_front_config<configuration_t, new_config_t>::type;

// ----------------------------------------------------------------------------
// configuration_fn_base
// ----------------------------------------------------------------------------

/*!\brief An abstract crtp-base class to add pipeable concept for configuration functors in combination with
          seqan3::detail::configuration.
 * \ingroup core_algorithm
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
        requires is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
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
        return {{std::forward<args_t>(args)...}};
    }
};

/*!\name Pipe Interface
 * \relates seqan3::detail::configuration_fn_base
 * \{
 */

/*!\brief Combines a seqan3::detail::configuration with an configuration adaptor.
 * \tparam configuration_t The configuration type to be extended. Must be of type seqan3::detail::configuration.
 * \tparam fn_t            Either an instance of seqan3::detail::configuration_fn_proxy or
 *                         seqan3::detail::configuration_fn_base
 *
 * \param[in] cfg The configuration to be extended.
 * \param[in] fn  The functor implementation, invoking the associated configuration adaptor.
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
 * \tparam lhs_fn_t The configuration adaptor to be combined.
 *                  Must be an instance of an instance of seqan3::detail::configuration_fn_base.
 * \tparam rhs_fn_t Either an instance of seqan3::detail::configuration_fn_proxy or
 *                  seqan3::detail::configuration_fn_base.
 *
 * \param[in] lhs_fn A configuration adaptor.
 * \param[in] rhs_fn The functor, invoking the associated configuration adaptor.
 * \returns The result of invoking the configuration adaptor with the passed `lhs_fn` object.
 */
template <typename lhs_fn_t,
          typename rhs_fn_t>
constexpr auto operator|(lhs_fn_t && lhs_fn,
                         rhs_fn_t && rhs_fn)
//!\cond
    requires std::is_base_of_v<configuration_fn_base<remove_cvref_t<lhs_fn_t>>, remove_cvref_t<lhs_fn_t>> &&
             (is_type_specialisation_of_v<remove_cvref_t<rhs_fn_t>, configuration_fn_proxy> ||
             (std::is_base_of_v<configuration_fn_base<remove_cvref_t<rhs_fn_t>>, remove_cvref_t<rhs_fn_t>> &&
              !std::is_same_v<remove_cvref_t<lhs_fn_t>, remove_cvref_t<rhs_fn_t>>))
//!\endcond
{
    return rhs_fn(std::invoke(std::forward<lhs_fn_t>(lhs_fn), configuration<>{}));
}

/*!\brief Combines a seqan3::detail::configuration_fn_proxy adaptor with another configuration adaptor.
 * \tparam proxy_fn_t The configuration adaptor proxy to be combined.
 *                    Must be an instance of an instance of seqan3::detail::configuration_fn_proxy.
 * \tparam rhs_fn_t   Either an instance of seqan3::detail::configuration_fn_proxy or
 *                    seqan3::detail::configuration_fn_base.
 *
 * \param[in] proxy_fn A configuration adaptor proxy.
 * \param[in] rhs_fn   The functor, invoking the associated configuration adaptor.
 * \returns The result of invoking the configuration adaptor with the passed `proxy_fn` object.
 */
template <typename proxy_fn_t,
          typename rhs_fn_t>
constexpr auto operator|(proxy_fn_t && proxy_fn,
                         rhs_fn_t && rhs_fn)
//!\cond
    requires is_type_specialisation_of_v<remove_cvref_t<proxy_fn_t>, configuration_fn_proxy> &&
             (is_type_specialisation_of_v<remove_cvref_t<rhs_fn_t>, configuration_fn_proxy> ||
             (std::is_base_of_v<configuration_fn_base<remove_cvref_t<rhs_fn_t>>, remove_cvref_t<rhs_fn_t>> &&
              !std::is_same_v<remove_cvref_t<proxy_fn_t>, remove_cvref_t<rhs_fn_t>>))
//!\endcond
{
    return rhs_fn(std::invoke(std::forward<proxy_fn_t>(proxy_fn), configuration<>{}));
}
//!\}

// ----------------------------------------------------------------------------
// configuration_fn_proxy
// ----------------------------------------------------------------------------

/*!\brief A proxy class used to defer invocation of the actual functor.
 * \ingroup core_algorithm
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
struct configuration_fn_proxy<derived_t, args_t...>
{
protected:

    //!\brief Friend declaration to allow only seqan3::detail::config_element_base access to the constructor.
    template <typename t>
    friend class configuration_fn_base;

    //!\brief Constructs this proxy while forwarding arguments.
    configuration_fn_proxy(args_t && ... args) : args_cache{std::forward<args_t>(args)...}
    {}

    /*!\name Helper functions
     * \{
     */
    template <typename configuration_t, std::size_t ... Is>
    auto explode(configuration_t && cfg, std::index_sequence<Is...> const & /*unsued*/) &&
    {
        return derived_t{}(std::forward<configuration_t>(cfg),
                           std::forward<args_t>(std::get<Is>(args_cache))...);
    }

    // Copy the temporaries at this place.
    template <typename configuration_t, std::size_t ... Is>
    auto explode(configuration_t && cfg, std::index_sequence<Is...> const & /*unsued*/) &
    {
        // TODO: copy the prvalues instead of forwarding them.
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
    auto operator()(configuration_t && cfg) &&
    {
        return explode(std::forward<configuration_t>(cfg), std::make_index_sequence<sizeof...(args_t)>{});
    }

    template <typename configuration_t>
    auto operator()(configuration_t && cfg) &
    {
        return explode(std::forward<configuration_t>(cfg), std::make_index_sequence<sizeof...(args_t)>{});
    }
    //!\}
};

// ----------------------------------------------------------------------------
// apply_deferred_configs
// ----------------------------------------------------------------------------

//!\cond
// TODO Document me!
template <std::size_t N, typename fn_t, typename config_t>
constexpr auto apply_deferred_configs(fn_t & fn,
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

        if constexpr (std::is_invocable_v<current_config_t, decltype(delegate), remove_cvref_t<config_t>>)  // deferred state.
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
             invocable_concept<std::remove_reference_t<fn_t>, std::remove_reference_t<config_t>>
{
    using type_list_t = typename std::remove_reference_t<config_t>::type_list_type;
    return apply_deferred_configs<meta::size<type_list_t>::value>(std::forward<fn_t>(fn),
                                                                  std::forward<config_t>(config));
}
//!\endcond
} // namespace seqan3::detail

namespace seqan3
{
/*!\name Tuple-like get interface
 * \relates seqan3::detail::configuration
 * \{
 */

/*!\brief Returns the state of the corresponding config at the specified position.
 * \tparam     elem_no   Non-type template parameter specifying the position of the config to get.
 * \tparam     configs_t Template parameter pack specifying the seqan3::detail::config_elements.
 * \param[in]  cfg       The configuration to query for it's config.
 *
 * \returns The state of the specified configuration.
 *
 * \details
 *
 * This function has the same semantics as std::get.
 *
 * ### Complexity
 *
 * Constant.
 *
 * ### Exception
 *
 * No-throw guarantee.
 *
 * ### Concurrency
 *
 * Thread-safe as long as the referenced data is not modified.
 */
template <size_t elem_no, typename ... configs_t>
constexpr auto & get(detail::configuration<configs_t...> & cfg) noexcept
{
    using type_list_type = typename detail::configuration<configs_t...>::type_list_type;

    static_assert(elem_no < meta::size<type_list_type>::value,
                  "Not requested position is to large.");
    return std::get<elem_no>(cfg.config_elements).data();
}

/*!\brief Returns the state of the corresponding config at the specified position.
 * \tparam     elem_no   Non-type template parameter specifying the position of the config to get.
 * \tparam     configs_t Template parameter pack specifying the seqan3::detail::config_elements.
 * \param[in]  cfg       The configuration to query for it's config.
 *
 * \returns The state of the specified configuration.
 *
 * \details
 *
 * This function has the same semantics as std::get.
 *
 * ### Complexity
 *
 * Constant.
 *
 * ### Exception
 *
 * No-throw guarantee.
 *
 * ### Concurrency
 *
 * Thread-safe.
 */
template <size_t elem_no, typename ... configs_t>
constexpr auto const & get(detail::configuration<configs_t...> const & cfg) noexcept
{
    using type_list_type = typename detail::configuration<configs_t...>::type_list_type;

    static_assert(elem_no < meta::size<type_list_type>::value,
                  "Not requested position is to large.");
    return std::get<elem_no>(cfg.config_elements).data();
}

/*!\brief Returns the state of the corresponding config at the specified position.
 * \tparam     elem_no   Non-type template parameter specifying the position of the config to get.
 * \tparam     configs_t Template parameter pack specifying the seqan3::detail::config_elements.
 * \param[in]  cfg       The configuration to query for it's config.
 *
 * \returns The state of the specified configuration.
 *
 * \details
 *
 * This function has the same semantics as std::get.
 *
 * ### Complexity
 *
 * Constant.
 *
 * ### Exception
 *
 * No-throw guarantee.
 *
 * ### Concurrency
 *
 * Thread-safe.
 */
template <size_t elem_no, typename ... configs_t>
constexpr auto && get(detail::configuration<configs_t...> && cfg) noexcept
{
    using type_list_type = typename detail::configuration<configs_t...>::type_list_type;

    static_assert(elem_no < meta::size<type_list_type>::value,
                  "Not requested position is to large.");
    return std::get<elem_no>(std::move(cfg.config_elements)).data();
}

/*!\brief Returns the state of the corresponding config at the specified position.
 * \tparam     elem_no   Non-type template parameter specifying the position of the config to get.
 * \tparam     configs_t Template parameter pack specifying the seqan3::detail::config_elements.
 * \param[in]  cfg       The configuration to query for it's config.
 *
 * \returns The state of the specified configuration.
 *
 * \details
 *
 * This function has the same semantics as std::get.
 *
 * ### Complexity
 *
 * Constant.
 *
 * ### Exception
 *
 * No-throw guarantee.
 *
 * ### Concurrency
 *
 * Thread-safe.
 */
template <size_t elem_no, typename ... configs_t>
constexpr auto const && get(detail::configuration<configs_t...> const && cfg) noexcept
{
    using type_list_type = typename detail::configuration<configs_t...>::type_list_type;

    static_assert(elem_no < meta::size<type_list_type>::value,
                  "Not requested position is to large.");

    // TODO: It is unclear why return must be wrapped in std::move here.
    return std::move(std::get<elem_no>(std::move(cfg.config_elements)).data());
}

/*!\brief Returns the state of the corresponding config at the specified position.
 * \tparam     target_t  Template parameter specifying the type of the config to get.
 * \tparam     configs_t Template parameter pack specifying the seqan3::detail::config_elements.
 * \param[in]  cfg       The configuration to query for it's config.
 *
 * \returns The state of the specified configuration.
 *
 * \details
 *
 * This function has the same semantics as std::get.
 *
 * ### Complexity
 *
 * Constant.
 *
 * ### Exception
 *
 * No-throw guarantee.
 *
 * ### Concurrency
 *
 * Thread-safe as long as the referenced data is not modified.
 */
template <typename target_t, typename ... configs_t>
constexpr auto & get(detail::configuration<configs_t...> & cfg) noexcept
{
    using type_list_type = typename detail::configuration<configs_t...>::type_list_type;
    constexpr size_t pos = meta::find_index<type_list_type, target_t>::value;
    static_assert(pos != meta::npos::value,
                  "The requested type does not exists in the list.");

    return get<pos>(cfg);
}

/*!\brief Returns the state of the corresponding config at the specified position.
 * \tparam     target_t  Template parameter specifying the type of the config to get.
 * \tparam     configs_t Template parameter pack specifying the seqan3::detail::config_elements.
 * \param[in]  cfg       The configuration to query for it's config.
 *
 * \returns The state of the specified configuration.
 *
 * \details
 *
 * This function has the same semantics as std::get.
 *
 * ### Complexity
 *
 * Constant.
 *
 * ### Exception
 *
 * No-throw guarantee.
 *
 * ### Concurrency
 *
 * Thread-safe.
 */
template <typename target_t, typename ... configs_t>
constexpr auto const & get(seqan3::detail::configuration<configs_t...> const & cfg) noexcept
{
    using type_list_type = typename detail::configuration<configs_t...>::type_list_type;
    constexpr size_t pos = meta::find_index<type_list_type, target_t>::value;
    static_assert(pos != meta::npos::value,
                  "The requested type does not exists in the list.");

    return get<pos>(cfg);
}

/*!\brief Returns the state of the corresponding config at the specified position.
 * \tparam     target_t  Template parameter specifying the type of the config to get.
 * \tparam     configs_t Template parameter pack specifying the seqan3::detail::config_elements.
 * \param[in]  cfg       The configuration to query for it's config.
 *
 * \returns The state of the specified configuration.
 *
 * \details
 *
 * This function has the same semantics as std::get.
 *
 * ### Complexity
 *
 * Constant.
 *
 * ### Exception
 *
 * No-throw guarantee.
 *
 * ### Concurrency
 *
 * Thread-safe.
 */
template <typename target_t, typename ... configs_t>
constexpr auto && get(detail::configuration<configs_t...> && cfg) noexcept
{
    using type_list_type = typename detail::configuration<configs_t...>::type_list_type;
    constexpr size_t pos = meta::find_index<type_list_type, target_t>::value;
    static_assert(pos != meta::npos::value,
                  "The requested type does not exists in the list.");

    return get<pos>(std::move(cfg));
}

/*!\brief Returns the state of the corresponding config at the specified position.
 * \tparam     target_t  Template parameter specifying the type of the config to get.
 * \tparam     configs_t Template parameter pack specifying the seqan3::detail::config_elements.
 * \param[in]  cfg       The configuration to query for it's config.
 *
 * \returns The state of the specified configuration.
 *
 * \details
 *
 * This function has the same semantics as std::get.
 *
 * ### Complexity
 *
 * Constant.
 *
 * ### Exception
 *
 * No-throw guarantee.
 *
 * ### Concurrency
 *
 * Thread-safe.
 */
template <typename target_t, typename ... configs_t>
constexpr auto const && get(detail::configuration<configs_t...> const && cfg) noexcept
{
    using type_list_type = typename detail::configuration<configs_t...>::type_list_type;
    constexpr size_t pos = meta::find_index<type_list_type, target_t>::value;
    static_assert(pos != meta::npos::value,
                  "The requested type does not exists in the list.");

    return get<pos>(std::move(cfg));
}
} // namespace seqan3

namespace std
{

/*!\brief Value metafunction specialisation for seqan3::detail::configuration
 * \ingroup core_algorithm
 * \relates seqan3::detail::configuration
 * \returns The number of configurations contained in seqan3::detail::configuration.
 */
template <typename ... configs_t>
struct tuple_size<seqan3::detail::configuration<configs_t...>>
{
    //!\brief The number of elements.
    static constexpr size_t value = sizeof...(configs_t);
};

/*!\brief Value metafunction specialisation for seqan3::detail::configuration.
 * \ingroup core_algorithm
 * \relates seqan3::detail::configuration
 * \returns The the type of configuration for the requested position.
 */
template <size_t elem_no, typename ... configs_t>
struct tuple_element<elem_no, seqan3::detail::configuration<configs_t...>>
{
    //!\brief The type of the element at position `elem_no`.
    using type = meta::at_c<typename seqan3::detail::configuration<configs_t...>::type_list_type, elem_no>;
};

//!\cond
//!\brief Provides std::get interface.
template <size_t elem_no, typename ... configs_t>
constexpr auto & get(seqan3::detail::configuration<configs_t...> & cfg) noexcept
{
    return seqan3::get<elem_no>(cfg);
}

//!\brief Provides std::get interface.
template <size_t elem_no, typename ... configs_t>
constexpr auto const & get(seqan3::detail::configuration<configs_t...> const & cfg) noexcept
{
    return  seqan3::get<elem_no>(cfg);
}

//!\brief Provides std::get interface.
template <size_t elem_no, typename ... configs_t>
constexpr auto && get(seqan3::detail::configuration<configs_t...> && cfg) noexcept
{
    return seqan3::get<elem_no>(std::move(cfg));
}

//!\brief Provides std::get interface.
template <size_t elem_no, typename ... configs_t>
constexpr auto const && get(seqan3::detail::configuration<configs_t...> const && cfg) noexcept
{
    return seqan3::get<elem_no>(std::move(cfg));
}

//!\brief Provides std::get interface.
template <typename target_t, typename ... configs_t>
constexpr auto & get(seqan3::detail::configuration<configs_t...> & cfg) noexcept
{
    return seqan3::get<target_t>(cfg);
}

//!\brief Provides std::get interface.
template <typename target_t, typename ... configs_t>
constexpr auto const & get(seqan3::detail::configuration<configs_t...> const & cfg) noexcept
{
    return  seqan3::get<target_t>(cfg);
}

//!\brief Provides std::get interface.
template <typename target_t, typename ... configs_t>
constexpr auto && get(seqan3::detail::configuration<configs_t...> && cfg) noexcept
{
    return seqan3::get<target_t>(std::move(cfg));
}

//!\brief Provides std::get interface.
template <typename target_t, typename ... configs_t>
constexpr auto const && get(seqan3::detail::configuration<configs_t...> const && cfg) noexcept
{
    return seqan3::get<target_t>(std::move(cfg));
}
//!\endcond
} // namespace std
