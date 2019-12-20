// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides the mode configuration to define the search modes "all", "all_best", "best" and "strata".
 * \author Christopher Pockrandt <christopher.pockrandt AT fu-berlin.de>
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <variant>

#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/algorithm/dynamic_state.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>
#include <seqan3/core/detail/strong_type.hpp>
#include <seqan3/core/type_traits/basic.hpp>
#include <seqan3/search/configuration/detail.hpp>

namespace seqan3::detail
{

//!\brief Type for the "all" value for the configuration element "mode".
//!\ingroup search_configuration
struct search_mode_all {};
//!\brief Type for the "all_best" value for the configuration element "mode".
//!\ingroup search_configuration
struct search_mode_all_best {};
//!\brief Type for the "best" value for the configuration element "mode".
//!\ingroup search_configuration
struct search_mode_best {};

} // namespace seqan3::detail

namespace seqan3::search_cfg
{

//!\brief Configuration element to receive all hits within the error bounds.
//!\ingroup search_configuration
inline detail::search_mode_all constexpr all;
//!\brief Configuration element to receive all hits within the lowest number of errors.
//!\ingroup search_configuration
inline detail::search_mode_all_best constexpr all_best;
//!\brief Configuration element to receive one best hit (with the lowest number of errors).
//!\ingroup search_configuration
inline detail::search_mode_best constexpr best;

/*!\brief Configuration element to receive all hits with the best number of errors plus the strata value.
 *        A strong type of underlying type `uint8_t` that represents the number or errors for strata.
 *        All hits are found with the fewest numbererrors plus 'value'.
 * \ingroup search_configuration
 * \tparam value_t The underlying type.
 */
struct strata : detail::strong_type<uint8_t, strata, detail::strong_type_skill::convert>
{
    using detail::strong_type<uint8_t, strata, detail::strong_type_skill::convert>::strong_type;
};

/*!\brief Configuration element to determine the search mode.
 * \ingroup search_configuration
 *
 * \details
 *
 * This configuration element can be used to determine which hits are supported.
 * Currently these modes are available:
 * | Mode                         | Behaviour                                                           |
 * |------------------------------|---------------------------------------------------------------------|
 * | seqan3::search_cfg::all      | Report all hits within error bounds.                                |
 * | seqan3::search_cfg::all_best | Report all hits with the lowest number of errors within the bounds. |
 * | seqan3::search_cfg::best     | Report one best hit (hit with lowest error) within bounds.          |
 * | seqan3::search_cfg::strata   | Report all hits within best + x errors.                             |
 *
 * ### Example
 *
 * \include test/snippet/search/configuration_modes.cpp
 *
 * ### Dynamic configuration
 *
 * In most use cases the search mode is fixed throughout the entire program. Therefore, it is sufficient to
 * initialise the seqan3::search_cfg::mode as shown in the example above. However, sometimes the mode might depend
 * on a runtime parameter, for example a program option set by the caller of the program. In this case the
 * seqan3::search_cfg::mode configuration can be set dynamically by first default constructing the object and then
 * assigning the corresponding search mode as depicted in the table above.
 * The following example demonstrates the dynamic configuration:
 *
 * \include test/snippet/search/configuration_modes_dynamic.cpp
 */
template <typename search_mode_t>
//!\cond
    requires std::same_as<search_mode_t, detail::search_mode_all> ||
             std::same_as<search_mode_t, detail::search_mode_all_best> ||
             std::same_as<search_mode_t, detail::search_mode_best> ||
             std::same_as<search_mode_t, strata> ||
             std::same_as<search_mode_t, dynamic_state>
//!\endcond
class mode : public pipeable_config_element<mode<search_mode_t>, search_mode_t>
{
private:
    //!\brief The type of the base class.
    using base_t = pipeable_config_element<mode<search_mode_t>, search_mode_t>;
    //!\brief A variant type over valid search modes.
    using search_modes_t = std::variant<detail::search_mode_all,
                                        detail::search_mode_all_best,
                                        detail::search_mode_best,
                                        strata>;

public:
    /*!\name Constructor, destructor and assignment
     * \{
     */
    constexpr mode() = default; //!< Enables dynamic configuration.
    constexpr mode(mode const &) = default; //!< Defaulted.
    constexpr mode(mode &&) = default; //!< Defaulted.
    constexpr mode & operator=(mode const &) = default; //!< Defaulted.
    constexpr mode & operator=(mode &&) = default; //!< Defaulted.
    ~mode() = default; //!< Defaulted.

    /*!\brief Construction of a fixed search mode with the given mode option.
     *
     * \param[in] new_search_mode The instance of the search mode option to use.
     *
     * \details
     *
     * Note when constructing the search mode directly with an instance of the valid mode options the mode cannot
     * be reassigned another value.
     */
    template <typename new_search_mode_t>
    //!\cond
        requires !std::same_as<new_search_mode_t, mode> &&
                 (std::same_as<new_search_mode_t, detail::search_mode_all> ||
                  std::same_as<new_search_mode_t, detail::search_mode_all_best> ||
                  std::same_as<new_search_mode_t, detail::search_mode_best> ||
                  std::same_as<new_search_mode_t, strata>)
    //!\endcond
    explicit constexpr mode(new_search_mode_t new_search_mode) : base_t{std::move(new_search_mode)}
    {}

    /*!\brief Assigns a new mode option to the search mode.
     *
     * \param[in] new_search_mode The instance of the search mode option to use.
     *
     * \details
     *
     * Note this assignment operator is only available if the search mode was default constructed.
     */
    template <typename new_search_mode_t>
    //!\cond
        requires !std::same_as<new_search_mode_t, mode> &&
                 (std::same_as<new_search_mode_t, detail::search_mode_all> ||
                  std::same_as<new_search_mode_t, detail::search_mode_all_best> ||
                  std::same_as<new_search_mode_t, detail::search_mode_best> ||
                  std::same_as<new_search_mode_t, strata>)
    //!\endcond
    constexpr mode & operator=(new_search_mode_t new_search_mode)
    {
        static_assert(std::same_as<search_mode_t, dynamic_state>,
                      "Expected this configuration to be dynamic, i.e. "
                      "std::same_as<search_mode_t, seqan3::dynamic_state> evaluates to true.");
        search_mode_variant = std::move(new_search_mode);
        return *this;
    }
    //!\}

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr detail::search_config_id id{detail::search_config_id::mode};

    //!\brief Returns the variant containing the search modes.
    constexpr search_modes_t const & search_modes() const noexcept
    //!\cond
        requires std::same_as<search_mode_t, dynamic_state>
    //!\endcond
    {
        return search_mode_variant;
    }

    //!\brief Evaluates to `true` if the search mode can be configured at runtime, otherwise `false`.
    static constexpr bool has_dynamic_state = std::same_as<search_mode_t, dynamic_state>;

private:
    //!\brief Stores the dynamically selected search mode option.
    search_modes_t search_mode_variant{};
};

/*!\name Type deduction guides
 * \relates seqan3::search_cfg::mode
 * \{
 */
//!\brief Deduces the dynamic search mode type.
mode() -> mode<dynamic_state>;

//!\brief Deduces search mode option type from constructor argument.
template <typename search_mode_t>
mode(search_mode_t) -> mode<search_mode_t>;
//!\}
} // namespace seqan3::search_cfg

namespace seqan3::detail
{
/*!\brief An internal representation of the associated seqan3::search_cfg::mode.
 * \ingroup search_configuration
 *
 * \tparam search_mode_t The type of the search mode option; must be one of seqan3::detail::search_mode_all,
 *                       seqan3::detail::search_mode_best, seqan3::detail::search_mode_all_best, or
 *                       seqan3::detail::search_mode_strata.
 *
 * \details
 *
 * This configuration wraps the public seqan3::search_cfg::mode. It is used to convert the dynamically configured
 * search mode into a static mode that can be used within the actual search algorithm implementation.
 * This is a trick to add the static information to the already existing algorithm configuration containing the
 * seqan3::search_cfg::mode configuration. Since, the same configuration element cannot be added multiple times a
 * distinct type is used to forward the actual information. This means that implementors need to use this configuration
 * instead of seqan3::search_cfg::mode when implementing the algorithm.
 */
template <typename search_mode_t>
//!\cond
    requires std::same_as<search_mode_t, search_mode_all> ||
             std::same_as<search_mode_t, search_mode_all_best> ||
             std::same_as<search_mode_t, search_mode_best> ||
             std::same_as<search_mode_t, search_cfg::strata>
//!\endcond
struct internal_search_mode : public pipeable_config_element<internal_search_mode<search_mode_t>, search_mode_t>
{
    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr search_config_id id{search_config_id::internal_search_mode};
};

/*!\name Type deduction guides
 * \relates seqan3::detail::internal_search_mode
 * \{
 */
//!\brief Deduces search mode type from constructor argument.
template <typename search_mode_t>
internal_search_mode(search_mode_t) -> internal_search_mode<search_mode_t>;
//!\}
}
