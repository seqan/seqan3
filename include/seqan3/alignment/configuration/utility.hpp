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
 * \brief Provides functionality to access get function by enum values.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <meta/meta.hpp>

#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/detail/reflection.hpp>

namespace seqan3::align_cfg
{

/*!\brief Specifies an id for every configuration element.
 * \ingroup configuration
 *
 * \details
 *
 * The seqan3::align_cfg::id is used to identify a specific alignment configuration element independent of
 * it's concrete type and position within the seqan3::align_cfg::alignment_configuration object.
 * Thus one can access the value of the corresponding configuration element via the special get interface.
 *
 * ### Example
 *
 * ```cpp
 * align_cfg::alignment_configuration cfg = align_cfg::gap_linear(gap_cost{-10});
 * auto cost = get<align_cfg::id::gap>(cfg);  // cost = -10;
 * ```
 */
enum struct id : uint8_t
{
    //!\brief Identifier for gap configuration.
    gap,
    //!\brief Identifier for global alignment configuration.
    global
    //!\cond
    // ATTENTION: Must always be the last item; will be used to determine the number of ids.
    ,SIZE
    //!\endcond
};

} // namespace seqan3::align_cfg

namespace seqan3::detail
{
// ----------------------------------------------------------------------------
// SEQAN3_INVALID_CONFIG
// ----------------------------------------------------------------------------

//!\cond
// A common error message for testing valid configuration combinations.
#define SEQAN3_INVALID_CONFIG(e) "Configuration error: The configuration <" #e "> is not combinable with "  \
                                 "one of the previous config elements. Please see the documentation to get "\
                                 "more information about which configurations can be combined."
//!\endcond

// ----------------------------------------------------------------------------
// on_align_config
// ----------------------------------------------------------------------------

/*!\brief Checks if a specific type corresponds to the given seqan3::align_cfg::id.
 * \ingroup configuration
 * \tparam e The id to test for.
 *
 * \details
 *
 * This class is an invocable trait meta-function which is used in combination with the meta::find_if type function.
 * It is used to find the configuration element type within the detail::configuration, which associated with the
 * seqan3::align_cfg::id given as non-type template parameter.
 * To use this, implementors need to overload this class with the specific seqan3::align_cfg::id and provide
 * a type alias called `invoke`, which resolves to a std::bool_constant (e.g. std::true_type).
 *
 * ### Example
 *
 * ```cpp
 * template <>
 * struct on_align_config<align_cfg::id::gap>
 * {
 *     template <typename t>
 *     using invoke = std::is_type_specialisation_of<t, align_config_gap>;
 * };
 * ```
 *
 * In this example `align_config_gap` is a template class representing the configuration element for the gap model.
 */
template <align_cfg::id e>
struct on_align_config
{
    /*!\brief A template alias for the invocation. Defaults to std::false_type.
     * \tparam t The type of the configuration element to test.
     */
    template <config_element_concept t>
    using invoke = std::false_type;
};

// ----------------------------------------------------------------------------
// align_config_type_to_id
// ----------------------------------------------------------------------------

/*!\brief Maps the given configuration element type to it's associated seqan3::align_cfg::id.
 * \ingroup configuration
 * \tparam config_element_t The type to get the mapped seqan3::align_cfg::id for.
 * \see seqan3::detail::align_config_type_to_id_v
 */
template <config_element_concept config_element_t>
struct align_config_type_to_id
{
    //!\brief The mapped seqan3::align_cfg::id. Defaults to seqan3::align_cfg::id::SIZE.
    static constexpr align_cfg::id value = align_cfg::id::SIZE;
};

//!\brief Helper variable template for seqan3::detail::align_config_type_to_id.
//!\ingroup configuration
template <config_element_concept config_element_t>
inline constexpr align_cfg::id align_config_type_to_id_v = align_config_type_to_id<config_element_t>::value;

// ----------------------------------------------------------------------------
// align_config_validation_matrix
// ----------------------------------------------------------------------------

/*!\brief Validation matrix to check how alignment configuration elements can be combined.
 * \ingroup configuration
 *
 * \details
 *
 * This matrix is used to check if the specified alignment configurations can be combined with each other.
 * A cell value `true`, indicates that the corresponding seqan3::align_cfg::id in the current column can be combined
 * with the associated seqan3::align_cfg::id in the current row. The size of the matrix is determined by the enum value
 * SIZE of seqan3::align_cfg::id.
 */
inline constexpr std::array<std::array<bool, static_cast<uint8_t>(align_cfg::id::SIZE)>,
                            static_cast<uint8_t>(align_cfg::id::SIZE)> align_config_validation_matrix =
{
    //  gap   global
    { { false, true },  // gap
      { true, false } } // global
};

/*!\brief Determines the first type in reverse order of the given detail::configuration that is not combinable with
 *        the given configuration identified by seqan3::align_cfg::id.
 * \ingroup configuration
 *
 * \tparam query           The seqan3::align_cfg::id to test for valid combination.
 * \tparam configuration_t The detail::configuration type.
 *
 * \returns The first type that is not combinable with the specified seqan3::align_cfg::id.
 * \see
 */
template <align_cfg::id query, typename configuration_t>
//!\cond
    requires is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
//!\endcond
struct invalid_alignment_configuration
{
protected:

    /*!\brief A type trait meta-function that checks if the given type is compatible with `query`.
     * \tparam target_t The type to check.
     */
    template <typename target_t>
    struct invoke
    {
        //!\brief Column index of `query`.
        static constexpr uint8_t pos1 = static_cast<uint8_t>(query);
        //!\brief Row index of the seqan3::align_cfg::id associated with `target_t`.
        static constexpr uint8_t pos2 = static_cast<uint8_t>(align_config_type_to_id_v<remove_cvref_t<target_t>>);

        static_assert(pos1 != static_cast<uint8_t>(align_cfg::id::SIZE),
                      "Unknown align_cfg::id! "
                      "Did you use an unknown config or did you forget to specialize align_config_type_to_id?");
        static_assert(pos2 != static_cast<uint8_t>(align_cfg::id::SIZE),
                      "Unknown align_cfg::id! "
                      "Did you use an unknown config or did you forget to specialize align_config_type_to_id?");

        //!\brief Type alias resolving to `std::false_type` if `target_t` can be combined with query, `std::true_type` otherwise.
        using type = std::conditional_t<align_config_validation_matrix[pos1][pos2], std::false_type, std::true_type>;
    };

    //!\brief seqan3::type_list alias for the types contained in `configuration_t`.
    using target_list_t = tuple_type_list_t<typename remove_cvref_t<configuration_t>::base_type>;
    //!\brief seqan3::type_list alias for tail list returned by meta::find_if.
    using tail_list_t   = meta::find_if<meta::reverse<target_list_t>, meta::quote_trait<invoke>>;
    //!\brief Defines a list of void element if the tail_list_t is empty.
    using final_list_t  = std::conditional_t<meta::empty<tail_list_t>::value, type_list<void>, tail_list_t>;

public:

    //!\brief The type alias for the type that is not compatible with query or void, if no conflicts are detected.
    using type = meta::front<final_list_t>;
};

//!\brief Helper template definition for seqan3::detail::invalid_alignment_configuration.
//!\ingroup configuration
template <align_cfg::id query, typename configuration_t>
//!\cond
    requires is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
//!\endcond
using invalid_alignment_configuration_t = typename invalid_alignment_configuration<query, configuration_t>::type;

/*!\brief Returns whether the given `query` is compatible with `configuration_t`.
 * \ingroup configuration
 *
 * \tparam query           The seqan3::align_cfg::id to check.
 * \tparam configuration_t The seqan3::detail::configuration type to test against.
 *
 * \returns std::true_type, if `query` is compatible, `std::false_type` otherwise.
 * \see seqan3::detail::is_valid_alignment_configuration_v
 */
template <align_cfg::id query, typename configuration_t>
//!\cond
    requires is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
//!\endcond
struct is_valid_alignment_configuration :
    public std::is_same<invalid_alignment_configuration_t<query, configuration_t>, void>
{};

//!\brief Helper variable template for seqan3::detail::is_valid_alignment_configuration.
//!\ingroup configuration
template <align_cfg::id query, typename configuration_t>
//!\cond
    requires is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
//!\endcond
inline constexpr bool is_valid_alignment_configuration_v =
    is_valid_alignment_configuration<query, configuration_t>::value;

} // namespace seqan3::detail

namespace seqan3
{

/*!\name Enum-based get-interface
 * \ingroup configuration
 * \brief Provides a special overload for the get interface using seqan3::align_cfg::id as identifier.
 * \{
 */

/*!\brief Access the value of the alignment configuration element identified by seqan3::align_cfg::id.
 * \relates seqan3::detail::configuration
 *
 * \tparam e                 The seqan3::align_cfg::id to get the stored configuration value for.
 * \tparam cfg_elements_t... The configuration elements stored in detail::alignment_configuration.
 *
 * \param[in] cfg            The configuration for which to get the configuration element.
 * \returns The stored configuration value associated with seqan3::align_cfg::id `e`.
 */
template <align_cfg::id e, typename ... cfg_elements_t>
constexpr auto & get(detail::configuration<cfg_elements_t...> & cfg) noexcept
{
    using type_list_t = detail::tuple_type_list_t<typename detail::configuration<cfg_elements_t...>::base_type>;

    constexpr size_t pos = meta::size<type_list_t>::value -
                           meta::size<meta::find_if<type_list_t, detail::on_align_config<e>>>::value;

    static_assert(static_cast<int8_t>(e) < static_cast<int8_t>(align_cfg::id::SIZE),
                  "Did you forget to update align_cfg::id::SIZE value?");
    static_assert(pos < meta::size<type_list_t>::value,
                  "The specified config element is not contained in the configuration.");

    return get<pos>(cfg).value;
}

/*!\brief Access the value of the alignment configuration element identified by seqan3::align_cfg::id.
 * \relates seqan3::detail::configuration
 *
 * \tparam e                 The seqan3::align_cfg::id to get the stored configuration value for.
 * \tparam cfg_elements_t... The configuration elements stored in detail::alignment_configuration.
 *
 * \param[in] cfg            The configuration for which to get the configuration element.
 * \returns The stored configuration value associated with seqan3::align_cfg::id `e`.
 */
template <align_cfg::id e, typename ... cfg_elements_t>
constexpr auto const & get(detail::configuration<cfg_elements_t...> const & cfg) noexcept
{
    using type_list_t = detail::tuple_type_list_t<typename detail::configuration<cfg_elements_t...>::base_type>;
    constexpr size_t pos = meta::size<type_list_t>::value -
                           meta::size<meta::find_if<type_list_t, detail::on_align_config<e>>>::value;

    static_assert(static_cast<int8_t>(e) < static_cast<int8_t>(align_cfg::id::SIZE),
                  "Did you forget to update align_cfg::id::SIZE value?");
    static_assert(pos < meta::size<type_list_t>::value,
                  "The specified config element is not contained in the configuration.");

    return get<pos>(cfg).value;
}

/*!\brief Access the value of the alignment configuration element identified by seqan3::align_cfg::id.
 * \relates seqan3::detail::configuration
 *
 * \tparam e                 The seqan3::align_cfg::id to get the stored configuration value for.
 * \tparam cfg_elements_t... The configuration elements stored in detail::alignment_configuration.
 *
 * \param[in] cfg            The configuration for which to get the configuration element.
 * \returns The stored configuration value associated with seqan3::align_cfg::id `e`.
 */
template <align_cfg::id e, typename ... cfg_elements_t>
constexpr auto && get(detail::configuration<cfg_elements_t...> && cfg) noexcept
{
    using type_list_t = detail::tuple_type_list_t<typename detail::configuration<cfg_elements_t...>::base_type>;
    constexpr size_t pos = meta::size<type_list_t>::value -
                           meta::size<meta::find_if<type_list_t, detail::on_align_config<e>>>::value;

    static_assert(static_cast<int8_t>(e) < static_cast<int8_t>(align_cfg::id::SIZE),
                  "Did you forget to update align_cfg::id::SIZE value?");
    static_assert(pos < meta::size<type_list_t>::value,
                  "The specified config element is not contained in the configuration.");

    return get<pos>(std::move(cfg)).value;
}

/*!\brief Access the value of the alignment configuration element identified by seqan3::align_cfg::id.
 * \relates seqan3::detail::configuration
 *
 * \tparam e                 The seqan3::align_cfg::id to get the stored configuration value for.
 * \tparam cfg_elements_t... The configuration elements stored in detail::alignment_configuration.
 *
 * \param[in] cfg            The configuration for which to get the configuration element.
 * \returns The stored configuration value associated with seqan3::align_cfg::id `e`.
 */
template <align_cfg::id e, typename ... cfg_elements_t>
constexpr auto const && get(detail::configuration<cfg_elements_t...> const && cfg) noexcept
{
    using type_list_t = detail::tuple_type_list_t<typename detail::configuration<cfg_elements_t...>::base_type>;
    constexpr size_t pos = meta::size<type_list_t>::value -
                           meta::size<meta::find_if<type_list_t, detail::on_align_config<e>>>::value;

    static_assert(static_cast<int8_t>(e) < static_cast<int8_t>(align_cfg::id::SIZE),
                  "Did you forget to update align_cfg::id::SIZE value?");
    static_assert(pos < meta::size<type_list_t>::value,
                  "The specified config element is not contained in the configuration.");

    return std::move(get<pos>(cfg).value); // TODO remove std::move when g++7 is fixed
}
//!\}
} // namespace seqan3
