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
#include <seqan3/core/algorithm/configuration_utility.hpp>
#include <seqan3/core/metafunction/basic.hpp>
#include <seqan3/core/metafunction/template_inspection.hpp>
#include <seqan3/core/tuple_utility.hpp>
#include <seqan3/core/type_list.hpp>
#include <seqan3/std/concepts>

namespace seqan3
{

// ----------------------------------------------------------------------------
// configuration
// ----------------------------------------------------------------------------

/*!\brief Collection of elements to configure an algorithm.
 * \ingroup algorithm
 *
 * \tparam configs_t Template parameter pack containing all configuration elements. Each must satisfy the
 *                   seqan3::detail::config_element_concept
 *
 * \details
 *
 * This class provides a unified interface to create and query such
 * configurations for a specific algorithm. It extends the standard tuple interface with some useful functions to modify
 * and query the user configurations.
 */
template <detail::config_element_concept ... configs_t>
class configuration : public std::tuple<configs_t...>
{
    //!\brief Friend declaration for other instances of the configuration.
    template <detail::config_element_concept ... _configs_t>
    friend class configuration;
public:
    //!\privatesection
    //!\brief A type alias for the base class.
    using base_type = std::tuple<configs_t...>;

    //!\publicsection
    /*!\name Constructor, destructor and assignment
     * \{
     */
    constexpr configuration()                                  = default;
    constexpr configuration(configuration const &)             = default;
    constexpr configuration(configuration &&)                  = default;
    constexpr configuration & operator=(configuration const &) = default;
    constexpr configuration & operator=(configuration &&)      = default;
    ~configuration()                                           = default;

    /*!\brief Constructs a configuration from a single configuration element.
     * \param elem The element to store.
     */
    template <typename config_t>
    constexpr configuration(config_t && elem) : base_type{std::forward<config_t>(elem)}
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

    // TODO Make to push_back
    /*!\brief Adds a new config element to the end of the configuration.
     *
     * \param[in] elem The configuration element to add.
     *
     * \returns A new seqan3::detail::configuration containing the added element.
     *
     * \details
     *
     * Creates a new seqan3::detail::configuration from `this` and appends the passed config element.
     * Note, that `this` is not modified by this operation.
     * Further the configuration checks for an invalid configuration using an algorithm specific lookup table
     * for the configuration elements and tests whether configuration elements are from the same algorithm.
     *
     * ### Example
     *
     * \snippet test/snippet/core/algorithm/configuration.cpp push_back
     *
     * ### Complexity
     *
     * Linear in the number of elements.
     *
     * ### Exception
     *
     * Strong exception guarantee.
     */
    template <detail::config_element_concept config_element_t>
    constexpr auto push_back(config_element_t elem) const &
    {
        static_assert(detail::is_configuration_valid_v<remove_cvref_t<config_element_t>,
                                                            configs_t...>,
                      "Configuration error: The passed element cannot be combined with one or more elements in the "
                      "current configuration.");

        return configuration<configs_t..., std::remove_reference_t<config_element_t>>{
            std::tuple_cat(static_cast<base_type>(*this),
            std::tuple{std::move(elem)})};
    }

    //!\copydoc push_back
    template <detail::config_element_concept config_element_t>
    constexpr auto push_back(config_element_t elem) &&
    {
        static_assert(detail::is_configuration_valid_v<remove_cvref_t<config_element_t>,
                                                            configs_t...>,
                      "Configuration error: The passed element cannot be combined with one or more elements in the "
                      "current configuration.");

        return configuration<configs_t..., std::remove_reference_t<config_element_t>>{
            std::tuple_cat(std::move(static_cast<base_type>(*this)),
            std::tuple{std::move(elem)})};
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
     * ### Complexity
     *
     * Linear in the number of elements.
     */
    template <detail::config_element_concept old_config_element_t,
              detail::config_element_concept new_config_element_t>
    constexpr auto replace_with(old_config_element_t const & SEQAN3_DOXYGEN_ONLY(old_element),
                                new_config_element_t && new_element) const &
    {
        static_assert(std::tuple_size_v<base_type> > 0, "The configuration cannot be empty.");
        static_assert(meta::find_index<detail::transfer_template_args_onto_t<base_type, type_list>,
                                       old_config_element_t>::value != meta::npos::value,
                      "The element to be replaced is not contained in the passed configuration.");

        auto && [prefix, remainder] = seqan3::tuple_split<old_config_element_t>(static_cast<base_type>(*this));

        return detail::configuration{std::tuple_cat(std::move(prefix),
                                                    std::tuple{std::forward<new_config_element_t>(new_element)},
                                                    tuple_pop_front(std::move(remainder)))};
    }

    //!\copydoc replace_with
    template <detail::config_element_concept old_config_element_t,
              detail::config_element_concept new_config_element_t>
    constexpr auto replace_with(old_config_element_t const & SEQAN3_DOXYGEN_ONLY(old_element),
                                new_config_element_t && new_element) &&
    {
        static_assert(std::tuple_size_v<base_type> > 0, "The configuration cannot be empty.");
        static_assert(meta::find_index<detail::transfer_template_args_onto_t<base_type, type_list>,
                                       old_config_element_t>::value != meta::npos::value,
                      "The element to be replaced is not contained in the passed configuration.");

        auto && [prefix, remainder] =
            seqan3::tuple_split<old_config_element_t>(std::move(static_cast<base_type>(*this)));

        return detail::configuration{std::tuple_cat(std::move(prefix),
                                                    std::tuple{std::forward<new_config_element_t>(new_element)},
                                                    tuple_pop_front(std::move(remainder)))};
    }
    //!\}

    /*!\name Observers
     * \{
     */

    /*!\brief Returns the contained value if `*this` has a value, otherwise returns `default_value`.
     * \tparam    query_t       The type to get the value from.
     * \param[in] default_value The default value if `query_t` is not contained in the configuration.
     *
     * \details
     *
     * Returns a reference to the stored configuration value by passing through the `.value` member of
     * the respective configuration element. If it does not exists than the default value is returned.
     * The existence check of a type is done at compile time.
     *
     * ### Example
     *
     * \snippet test/snippet/core/algorithm/configuration.cpp value_or
     *
     * ### Exception
     *
     * no-throw guarantee.
     *
     * ### Complexity
     *
     * Constant time.
     */
    template <typename query_t, typename default_t>
    constexpr decltype(auto) value_or(default_t && default_value) & noexcept
    {
        if constexpr (exists<query_t>())
            return get<query_t>(*this).value;
        else
            return std::forward<default_t>(default_value);
    }

    //!\copydoc value_or
    template <typename query_t, typename default_t>
    constexpr decltype(auto) value_or(default_t && default_value) const & noexcept
    {
        if constexpr (exists<query_t>())
            return get<query_t>(*this).value;
        else
            return std::forward<default_t>(default_value);
    }

    //!\copydoc value_or
    template <typename query_t, typename default_t>
    constexpr decltype(auto) value_or(default_t && default_value) && noexcept
    {
        if constexpr (exists<query_t>())
            return get<query_t>(std::move(*this)).value;
        else
            return std::forward<default_t>(default_value);
    }

    //!\copydoc value_or
    template <typename query_t, typename default_t>
    constexpr decltype(auto) value_or(default_t && default_value) const && noexcept
    {
        if constexpr (exists<query_t>())
            return get<query_t>(std::move(*this)).value;
        else
            return std::forward<default_t>(default_value);
    }

    //!\copydoc value_or
    template <template <typename...> typename query_t, typename default_t>
    constexpr decltype(auto) value_or(default_t && default_value) & noexcept
    {
        if constexpr (exists<query_t>())
            return get<query_t>(*this).value;
        else
            return std::forward<default_t>(default_value);
    }

    //!\copydoc value_or
    template <template <typename...> typename query_t, typename default_t>
    constexpr decltype(auto) value_or(default_t && default_value) const & noexcept
    {
        if constexpr (exists<query_t>())
            return get<query_t>(*this).value;
        else
            return std::forward<default_t>(default_value);
    }

    //!\copydoc value_or
    template <template <typename...> typename query_t, typename default_t>
    constexpr decltype(auto) value_or(default_t && default_value) && noexcept
    {
        if constexpr (exists<query_t>())
            return get<query_t>(std::move(*this)).value;
        else
            return std::forward<default_t>(default_value);
    }

    //!\copydoc value_or
    template <template <typename...> typename query_t, typename default_t>
    constexpr decltype(auto) value_or(default_t && default_value) const && noexcept
    {
        if constexpr (exists<query_t>())
            return get<query_t>(std::move(*this)).value;
        else
            return std::forward<default_t>(default_value);
    }

    //!\brief Checks if the given type exists in the tuple.
    template <typename query_t>
    static constexpr bool exists() noexcept
    {
        return !meta::empty<meta::find<type_list<configs_t...>, query_t>>::value;
    }
    //!\brief Checks if the given type exists in the tuple.
    template <template <typename...> typename query_t>
    static constexpr bool exists() noexcept
    {
        return !meta::empty<meta::find_if<type_list<configs_t...>, detail::is_same_configuration_f<query_t>>>::value;
    }
    //!\}

protected:

    /*!\name Internal constructor
     * \{
     */
    //!\brief Constructs from std::tuple.
    template <typename ... _configs_t>
    explicit constexpr configuration(std::tuple<_configs_t...> const & cfg) : base_type{cfg}
    {}

    //!\brief Constructs from std::tuple.
    template <typename ... _configs_t>
    explicit constexpr configuration(std::tuple<_configs_t...> && cfg) : base_type{std::move(cfg)}
    {}
    //!\}
};

/*!\name Type deduction guides
 * \relates seqan3::configuration
 * \{
 */
//!\brief Deduces the correct configuration element type from the passed element.
template <typename config_t>
configuration(config_t &&) -> configuration<remove_cvref_t<config_t>>;
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
 * \snippet test/snippet/core/algorithm/configuration.cpp get
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
    using _tail = meta::find_if<type_list<configs_t...>, detail::is_same_configuration_f<query_t>>;
    static_assert(!meta::empty<_tail>::value, "Access error: The requested type is not contained.");

    constexpr size_t pos = sizeof...(configs_t) - meta::size<_tail>::value;
    return get<pos>(config);
}

template <template <typename ...> class query_t, typename ...configs_t>
constexpr auto const & get(configuration<configs_t...> const & config) noexcept
{
    using _tail = meta::find_if<type_list<configs_t...>, detail::is_same_configuration_f<query_t>>;
    static_assert(!meta::empty<_tail>::value, "Access error: The requested type is not contained.");

    constexpr size_t pos = sizeof...(configs_t) - meta::size<_tail>::value;
    return get<pos>(config);
}

template <template <typename ...> class query_t, typename ...configs_t>
constexpr auto && get(configuration<configs_t...> && config) noexcept
{
    using _tail = meta::find_if<type_list<configs_t...>, detail::is_same_configuration_f<query_t>>;
    static_assert(!meta::empty<_tail>::value, "Access error: The requested type is not contained.");

    constexpr size_t pos = sizeof...(configs_t) - meta::size<_tail>::value;
    return get<pos>(std::move(config));
}

template <template <typename ...> class query_t, typename ...configs_t>
constexpr auto const && get(configuration<configs_t...> const && config) noexcept
{
    using _tail = meta::find_if<type_list<configs_t...>, detail::is_same_configuration_f<query_t>>;
    static_assert(!meta::empty<_tail>::value, "Access error: The requested type is not contained.");

    constexpr size_t pos = sizeof...(configs_t) - meta::size<_tail>::value;
    // TODO: change after GCC-7 bug with const && version of get in std::tuple is fixed.
    // return get<pos>(std::move(config));
    return std::move(get<pos>(config));
}
//!\}

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
    requires detail::is_algorithm_configuration_v<remove_cvref_t<config_t>>
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
    requires detail::is_algorithm_configuration_v<remove_cvref_t<config_t>> &&
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
struct tuple_size<seqan3::configuration<configs_t...>>
{
    //!\brief The number of elements.
    static constexpr size_t value = std::tuple_size_v<typename seqan3::configuration<configs_t...>::base_type>;
};

/*!\brief Returns the type of the element at the specified position within seqan3::configuration.
 * \ingroup algorithm
 */
template <size_t pos, seqan3::detail::config_element_concept ... configs_t>
struct tuple_element<pos, seqan3::configuration<configs_t...>>
{
    //!\brief The type of the config at position `pos`
    using type = std::tuple_element_t<pos, typename seqan3::configuration<configs_t...>::base_type>;
};
//!\endcond
}  //namespace std
