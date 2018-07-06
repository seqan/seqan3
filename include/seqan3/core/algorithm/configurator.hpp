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
 * \brief Provides seqan3::detail::configurator.
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

namespace seqan3::detail
{
// ----------------------------------------------------------------------------
// configurator_base
// ----------------------------------------------------------------------------

/*!\brief Abstract base class for the seqan3::detail::configurator class.
 * \tparam configs_t Template parameter pack over all configs. These types must satisfy the
 *                   seqan3::detail::config_concept.
 * \ingroup core_algorithm
 *
 * \details
 *
 * This is an abstract base class adding a constructor to the seqan3::detail::configurator class with special copy semantics:
 * When a configurator is extended by a new config type it must be placed at the beginning. To copy the states of
 * the old configurator, i.e. it contains the states of the already contained config types, the constructors of the
 * old config types must be called with the previous configurator to properly submit the states to the extended
 * configurator instance.
 */

template <typename ... configs_t>
class configurator_base;

/*!\fn template <typename ... other_configs_t> configurator_base::configurator_base(configurator<other_configs_t...> const & cfg)
 * \brief Constructs a new configurator_base object while delegating to the constructors of all inherited
 *        sub-configurations.
 * \memberof seqan3::detail::configurator_base
 * \tparam remaining_configs_t A template parameter pack containing all the sub-configurations of the previous instance
 *         of the configurator.
 * \param cfg The previous configurator before extending it with a new configuration type.
 *
 * \details
 *
 * The constructor checks if there exists a common reference between `cfg` and the `remaining_configs_t` types.
 * That is this class can only be constructed from `cfg` if and only if for every type in `remaining_configs_t`
 * there is a corresponding type contained in `cfg`.
 */
//!\cond
template <detail::config_concept config_head_t, detail::config_concept ... remaining_configs_t>
class configurator_base<config_head_t, remaining_configs_t...> : public config_head_t,
                                                                 public remaining_configs_t...
{
protected:

    configurator_base()                                      = default;
    configurator_base(configurator_base const &)             = default;
    configurator_base(configurator_base &&)                  = default;
    configurator_base & operator=(configurator_base const &) = default;
    configurator_base & operator=(configurator_base &&)      = default;
    ~configurator_base()                                     = default;

    template <typename ... other_configs_t>
        requires (common_reference_concept<configurator<other_configs_t ...>, remaining_configs_t> && ...)
    constexpr explicit configurator_base(configurator<other_configs_t ...> const & cfg) :
        remaining_configs_t{cfg}...
    {}
};
//!\endcond

// ----------------------------------------------------------------------------
// configurator
// ----------------------------------------------------------------------------

/*!\brief Collection of configurations objects used to specify the run behavior of algorithms.
 * \ingroup core_algorithm
 * \implements seqan3::detail::configurator_concept
 * \tparam configs_t Template parameter pack containing all the inherited configs. Must satisfy the
 *                   seqan3::detail::config_concept
 *
 * \details
 *
 * This class provides a unified interface and additional helper functions to create and query configurations
 * for a specific algorithm. Certain bioinformatics algorithms, e.g. alignment or find interfaces, support a various set
 * of different configurations and policies that alter the execution of the algorithm. These configurations can be
 * orthogonal or might be mutual exclusive. Using this configurator the interface for the user becomes much more easier,
 * and incompatible configurations can be checked at compile time.
 *
 * \attention
 *
 * This class uses multiple inheritance to inherit from the actual configuration types. This also means, that the same
 * configuration type cannot occur more than once in the type definition of the configurator class.
 *
 * ### Usage
 *
 * The entire configuration system is designed to work completely in the background of any algorithm. The implementor
 * of an algorithm only needs to provide the configurations for the algorithm.
 * In general the type of a configuration is static (see seqan3::detail::config_base) within the context of an algorithm
 * to select the correct code branches based on a policy-driven design, albeit the stored state might still be a runtime
 * parameter.
 * However, in some cases a specific configuration might be known first at runtime but needs to be converted to a static
 * type for the algorithm in use. To enable a transparent conversion from the runtime parameters to a static type
 * configuration one can use deferred configs (seqan3::detail::deferred_config_base), which are invocable
 * configurations, that store the runtime parameter and on invocation translate this runtime parameter to a static
 * configuration. The implementor of an algorithm can achieve this by using the function
 * seqan3::detail::apply_deferred_configs, which iterates through all configurations and in case
 * it is a deferred configuration it will invoke the translation function and continue with the modified configurator,
 * which contains now the static type for the specific configuration.
 *
 * ### Pipe Notation
 *
 * To enable simple extension of configurations the configurator provides a generic pipe interface for the different
 * configurations. Thus, a config class can be easily constructed by chaining together different properties.
 * Consider the following example, where we assume that the config `bar` and the config functor `with_foo` are already
 * given:
 *
 * ```cpp
 * auto my_cfg = configurator<bar>{} | with_foo;  // my_cfg is now of type configurator<foo, bar>
 * ```
 * The `with_foo` functor captures the `configurator<bar>` object and returns a new configurator with the added
 * `foo`-config.

 * ### Accessor
 *
 * The configurator exposes a tuple interface. Thus, one can either use `seqan3::get` or `std::get` to query
 * the state of a specific config within this configurator. Considering the example from above one can get the state of
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
template <detail::config_concept ... configs_t>
//!\cond
    requires sizeof...(configs_t) >= 1
//!\endcond
class configurator : public detail::configurator_base<configs_t...>

{
    //!\brief Typedef for the configurator type containing all but the first configuration.
    // using prev_configurator_type = std::conditional_t<sizeof...(configs_t) > 1,
    //                                                  detail::transfer_template_args_onto_t<
    //                                                         meta::pop_front<type_list<configs_t...>>, configurator>,
    //                                                  decltype(std::ignore)>;
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
    configurator()                                 = default;
    configurator(configurator const &)             = default;
    configurator(configurator &&)                  = default;
    configurator & operator=(configurator const &) = default;
    configurator & operator=(configurator &&)      = default;
    ~configurator()                                = default;

    /*!\brief Copy constructs a configurator from the configurator prior to the modification with the new config.
     * \tparam other_configs_t The other configs.
     * \param cfg The configurator immediately prior to it's extension by another configuration.
     *
     * \details
     *
     * Note, that this constructor is only enabled for types that are in direct relation to this class, namely for
     * types for which the following invariant holds: `configurator<configs...>` corresponds to
     * `configurator<head, configs...>`, with `decltype(*this) = configurator<head, configs...>`.
     */
    template <typename ... other_configs_t>
    constexpr explicit configurator(configurator<other_configs_t ...> const & cfg) :
        detail::configurator_base<configs_t ...>(cfg)
    {}

    //!\}

    /*!\name Accessor
     * \{
     */

    /*!\brief Returns the state of the corresponding config at the specified position.
     * \tparam elem_no Non-type template parameter specifying the position of the config to get.
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
    template <size_t elem_no>
    constexpr auto & get() & noexcept
    {
        static_assert(elem_no < meta::size<type_list_type>::value,
                      "Not requested position is to large.");
        return static_cast<std::tuple_element_t<elem_no, configurator>&>(*this).data();
    }

    //!\copydoc configurator::get()
    template <size_t elem_no>
    constexpr auto const & get() const & noexcept
    {
        static_assert(elem_no < meta::size<type_list_type>::value,
                      "Not requested position is to large.");
        return static_cast<std::tuple_element_t<elem_no, configurator> const &>(*this).data();
    }

    //!\copydoc configurator::get()
    template <size_t elem_no>
    constexpr auto && get() && noexcept
    {
        static_assert(elem_no < meta::size<type_list_type>::value,
                      "Not requested position is to large.");
        return std::move(static_cast<std::tuple_element_t<elem_no, configurator> &&>(*this).data());
    }

    //!\copydoc configurator::get()
    template <size_t elem_no>
    constexpr auto const && get() const && noexcept
    {
        static_assert(elem_no < meta::size<type_list_type>::value,
                      "Not requested position is to large.");
        return std::move(static_cast<std::tuple_element_t<elem_no, configurator> const &&>(*this).data());
    }

    /*!\brief Returns the state of the corresponding config that matches the specified type.
     * \tparam type Template parameter specifying the type of the config to get.
     * \returns The first config, whose type matches the specified one.
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
    template <typename target_t>
    constexpr auto & get() & noexcept
    {
        constexpr size_t pos = meta::find_index<type_list_type, target_t>::value;
        static_assert(pos != meta::npos::value,
                      "The requested type does not exists in the list.");
        return static_cast<std::tuple_element_t<pos, configurator>&>(*this).data();
    }

    //!\copydoc configurator::get()
    template <typename target_t>
    constexpr auto const & get() const & noexcept
    {
        constexpr size_t pos = meta::find_index<type_list_type, target_t>::value;
        static_assert(pos != meta::npos::value,
                      "The requested type does not exists in the list.");
        return static_cast<std::tuple_element_t<pos, configurator> const &>(*this).data();
    }

    //!\copydoc configurator::get()
    template <typename target_t>
    constexpr auto && get() && noexcept
    {
        constexpr size_t pos = meta::find_index<type_list_type, target_t>::value;
        static_assert(pos != meta::npos::value,
                      "The requested type does not exists in the list.");
        return std::move(static_cast<std::tuple_element_t<pos, configurator> &&>(*this).data());
    }

    //!\copydoc configurator::get()
    template <typename target_t>
    constexpr auto const && get() const && noexcept
    {
        constexpr size_t pos = meta::find_index<type_list_type, target_t>::value;
        static_assert(pos != meta::npos::value,
                      "The requested type does not exists in the list.");
        return std::move(static_cast<std::tuple_element_t<pos, configurator> const &&>(*this).data());
    }
    //!\}
};

//!\cond
template <typename configurator_t, detail::config_concept old_config_t, detail::config_concept new_config_t>
struct replace_config_with;
//!\endcond

/*!\brief Helper traits meta-function to replace one configuration type with another by pushing it to the front of
 *        the configurations for a specific seqan3::detail::configurator.
 * \ingroup core_algorithm
 * \relates seqan3::detail::configurator
 *
 * \tparam configurator_t The configurator type to be altered. Must be a template class.
 * \tparam old_config_t   The old_config_t to be removed.
 * \tparam new_config_t   The new_config_t to be added at the front of the configurations.
 *
 * \details
 *
 * This meta-function parses the list of config types contained in the given seqan3::detail::configurator and
 * removes the `old_config_t` from the list if it exists, and adds `new_config_t` at the front of
 * the type list. Subsequently a new `configurator_t` type is defined containing the new config type list.
 *
 * \see seqan3::detail::replace_config_with_t
 */
template <template <typename ...> class configurator_t,
          detail::config_concept old_config_t,
          detail::config_concept new_config_t,
          typename ... configs_t>
struct replace_config_with<configurator_t<configs_t...>, old_config_t, new_config_t>
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

    //!\brief Typedef with the modified condigurator type.
    using type = detail::transfer_template_args_onto_t<_new_list, configurator_t>;
};

/*!\brief Shortcut typedef for seqan3::detail::replace_config_with.
 * \relates seqan3::detail::replace_config_with
 */
template <typename configurator_t, typename old_config_t, typename new_config_t>
using replace_config_with_t = typename replace_config_with<configurator_t, old_config_t, new_config_t>::type;

// ----------------------------------------------------------------------------
// configurator_fn_base
// ----------------------------------------------------------------------------

/*!\brief An abstract crtp-base class to add pipeable concept for configuration functors in combination with
          seqan3::detail::configurator.
 * \ingroup core_algorithm
 * \tparam derived_t The derived type to be extended with functionality from this base class.
 */
template <typename derived_fn_t>
class configurator_fn_base
{
    //!\brief Empty base class to enable checking for the correct proxy type.
    struct configurator_fn_proxy_base
    {};

    /*!\brief A proxy class used to defer invocation of the actual functor.
     * \tparam args_t Template parameter pack with intermediate arguments that should be applied on invocation.
     */
    template <typename ... args_t>
    struct configurator_fn_proxy : configurator_fn_proxy_base
    {
        //!\brief Constructs this proxy while forwarding arguments.
        configurator_fn_proxy(args_t && ... args) : args_cache{std::forward<args_t>(args)...}
        {}

        /*!\name Invocable interface
         * \{
         */
        template <typename configurator_t>
        auto operator()(configurator_t && cfg) &&
        {
            return explode(cfg, std::make_index_sequence<sizeof...(args_t)>{});
        }

        template <typename configurator_t>
        auto operator()(configurator_t && cfg) &
        {
            return explode(cfg, std::make_index_sequence<sizeof...(args_t)>{});
        }
        //!\}

    protected:

        /*!\name Helper functions
         * \{
         */
        template <typename configurator_t, std::size_t ... Is>
        auto explode(configurator_t && cfg, std::index_sequence<Is...> const & /*unsued*/) &&
        {
            return configurator_fn_base{}(std::forward<configurator_t>(cfg),
                                          std::forward<args_t>(std::get<Is>(args_cache))...);
        }

        // Copy the temporaries at this place.
        template <typename configurator_t, std::size_t ... Is>
        auto explode(configurator_t && cfg, std::index_sequence<Is...> const & /*unsued*/) &
        {
            // copy the prvalues instead of forwarding them.
            return configurator_fn_base{}(std::forward<configurator_t>(cfg), std::get<Is>(args_cache)...);
        }
        //!\}

        //!\brief The cached data.
        std::tuple<args_t...> args_cache;
    };

public:

    /*!\brief Invokes the configuration specific functor to extend the seqan3::detail::configurator with the associated config.
     * \tparam configurator_t The configurator to be extended with the new configuration type. Must satisfy the
     *                        seqan3::detail::configurator_concept.
     * \tparam args_t A template parameter pack with the arguments applied to the configuration functor call.
     *
     * \param configurator The configurator to be extended.
     * \param args         The argument pack to be passed to the functor.
     * \returns A new seqan3::detail::configurator with the extended configuration.
     */
    template <detail::configurator_concept configurator_t,
              typename ... args_t>
    constexpr auto operator()(configurator_t && configurator,
                              args_t && ... args) const
    {
        return static_cast<derived_fn_t const &>(*this).invoke(std::forward<configurator_t>(configurator),
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
    constexpr configurator_fn_proxy<args_t...> operator()(args_t && ... args) const
    {
        return {{std::forward<args_t>(args)...}};
    }

    /*!\name Pipe Interface
     * \{
     */

     /*!\brief Invokes the configuration functor with the given seqan3::detail::configurator.
      * \tparam configurator_t The configurator type to be extended. Must satisfy seqan3::detail::configurator_concept.
      * \tparam configurator_fn_proxy_t The proxy type which defers the functor invocation.
      *                                 Must be derived from configurator_fn_proxy_base.
      *
      * \param cfg      The configurator.
      * \param proxy_fn The proxy.
      * \returns The result of invoking the configuration functor with the passed `cfg` object and the arguments cached
      *          in `proxy_fn`.
      */
    template <detail::configurator_concept configurator_t,
              typename configurator_fn_proxy_t>
    friend constexpr auto operator|(configurator_t && cfg,
                                    configurator_fn_proxy_t && proxy_fn)
        requires std::is_base_of_v<configurator_fn_proxy_base, configurator_fn_proxy_t>
    {
        return proxy_fn(std::forward<configurator_t>(cfg));
    }

    /*!\brief Invokes the configuration functor with the given seqan3::detail::configurator.
     * \tparam configurator_t The configurator type to be extended. Must satisfy seqan3::detail::configurator_concept.
     * \tparam derived_fn_t The derived type that adds the configuration.
     *
     * \param cfg The configurator.
     * \param fn  The functor that adds the configuration to `cfg`.
     * \returns The result of invoking the configuration functor with the passed `cfg`.
     */
    template <detail::configurator_concept configurator_t>
    friend constexpr auto operator|(configurator_t && cfg,
                                    derived_fn_t && fn)
    {
        return fn(std::forward<configurator_t>(cfg));
    }
    //!\}
};

// ----------------------------------------------------------------------------
// apply_deferred_configs
// ----------------------------------------------------------------------------

//!\cond
template <std::size_t N, typename fn_t, detail::configurator_concept config_t>
constexpr auto apply_deferred_configs(fn_t & fn,
                                      config_t && config)
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

template <typename fn_t, detail::configurator_concept config_t>
constexpr auto apply_deferred_configs(fn_t & fn,
                                      config_t && config)
    requires invocable_concept<std::remove_reference_t<fn_t>, std::remove_reference_t<config_t>>
{
    using type_list_t = typename std::remove_reference_t<config_t>::type_list_type;
    return apply_deferred_configs<meta::size<type_list_t>::value>(std::forward<fn_t>(fn),
                                                                  std::forward<config_t>(config));
}
//!\endcond
} // namespace seqan3::detail
