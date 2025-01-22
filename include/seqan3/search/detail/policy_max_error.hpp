// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides the seqan3::detail::policy_max_error.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 * \author Lydia Buntrock <lydia.buntrock AT fu-berlin.de>
 */

#pragma once

#include <ranges>

#include <seqan3/search/detail/search_common.hpp>
#include <seqan3/search/detail/search_traits.hpp>

namespace seqan3::detail
{

//!\brief Provides the function `max_error_counts` if inherited by a search algorithm.
//!\ingroup search
struct policy_max_error
{
protected:
    //!\brief The total errors set by the user.
    search_cfg::max_error_total total{};
    //!\brief The substitution errors set by the user.
    search_cfg::max_error_substitution substitution{};
    //!\brief The insertion errors set by the user.
    search_cfg::max_error_insertion insertion{};
    //!\brief The deletion errors set by the user.
    search_cfg::max_error_deletion deletion{};

    //!\brief Flag indicating if only max error was given.
    bool only_max_error_total{false};
    //!\brief Flag indicating if max error total was given.
    bool has_max_error_total{false};

    /*!\name Constructors, destructor and assignment
     * \{
     */
    policy_max_error() = default;                                     //!< Defaulted.
    policy_max_error(policy_max_error const &) = default;             //!< Defaulted.
    policy_max_error(policy_max_error &&) = default;                  //!< Defaulted.
    policy_max_error & operator=(policy_max_error const &) = default; //!< Defaulted.
    policy_max_error & operator=(policy_max_error &&) = default;      //!< Defaulted.
    ~policy_max_error() = default;                                    //!< Defaulted.

    /*!\brief Initialises the policy with the given configuration.
     *
     * \tparam configuration_t The configuration type; must be an instance of seqan3::configuration.
     * \param[in] config The search configuration object.
     *
     * \details
     *
     * Initialises the maximal errors for the respective edit operations from the given configuration.
     */
    template <typename configuration_t>
        requires is_type_specialisation_of_v<configuration_t, seqan3::configuration>
    explicit policy_max_error(configuration_t const & config)
    {
        using search_traits_t = search_traits<configuration_t>;
        only_max_error_total = search_traits_t::only_max_error_total;
        has_max_error_total = search_traits_t::has_max_error_total;

        total = config.get_or(search_cfg::max_error_total{search_cfg::error_count{0}});
        substitution = config.get_or(search_cfg::max_error_substitution{search_cfg::error_count{0}});
        insertion = config.get_or(search_cfg::max_error_insertion{search_cfg::error_count{0}});
        deletion = config.get_or(search_cfg::max_error_deletion{search_cfg::error_count{0}});
    }
    //!\}

    /*!\brief Returns a detail::search_param object filled by the information from the configuration.
     *
     * \tparam query_t Must model std::ranges::forward_range over the index's alphabet.
     *
     * \param[in] query The current query sequence.
     *
     * \throws std::invalid_argument
     */
    template <std::ranges::forward_range query_t>
    auto max_error_counts(query_t && query)
    {
        detail::search_param errors{0, 0, 0, 0}; // total, substitution, insertion, deletion

        [[maybe_unused]] auto query_size = std::ranges::size(query);

        errors.total = to_error_count(total.error, query_size);
        errors.substitution = to_error_count(substitution.error, query_size);
        errors.insertion = to_error_count(insertion.error, query_size);
        errors.deletion = to_error_count(deletion.error, query_size);

        // If only total is set, we set all other errors to the total limit.
        if (only_max_error_total)
            errors.substitution = errors.insertion = errors.deletion = errors.total;
        // If total is not set but any other field is set than use total as the sum of all set errors.
        else if (!has_max_error_total)
            errors.total = std::min<uint32_t>(255, errors.substitution + errors.insertion + errors.deletion);

        // Validate the error configuration.
        // Checks if the given thresholds for the configured errors are valid. Otherwise throws std::invalid_argument.
        if (errors.substitution > errors.total)
            throw std::invalid_argument{"The substitution error threshold is higher than the total error threshold."};
        if (errors.insertion > errors.total)
            throw std::invalid_argument{"The insertion error threshold is higher than the total error threshold."};
        if (errors.deletion > errors.total)
            throw std::invalid_argument{"The deletion error threshold is higher than the total error threshold."};

        return errors;
    }

private:
    /*!\brief Returns a uint8_t object directly given by the search_cfg::error_count or calculated by the
     *        search_cfg::error_rate depending on the set alternative.
     *
     * \param[in] error_variant A std::variant over search_cfg::error_count and search_cfg::error_rate.
     * \param[in] query_size    The size of the query.
     *
     * \throws std::invalid_argument if error_rate is not between 0 and 1.
     *
     * \details
     *
     * If seqan3::search_cfg::error_count is set in the configuration, its value already is of type uint8_t and can be
     * returned as is.
     * If seqan3::search_cfg::error_rate is set in the configuration, the error rates are converted to error counts
     * based on the length of the query sequence.
     */
    uint8_t to_error_count(std::variant<search_cfg::error_count, search_cfg::error_rate> const & error_variant,
                           [[maybe_unused]] size_t const query_size)
    {
        return std::visit(
            [&query_size](auto error)
            {
                if constexpr (std::same_as<decltype(error), search_cfg::error_count>)
                    return error.get();
                else
                {
                    // check correct error rate values.
                    if (0.0 > error.get() || error.get() > 1.0)
                        throw std::invalid_argument{"Error rates must be between 0 and 1."};

                    // make sure that error rate can be saved as uint8_t, so it is not too big in terms of query size
                    uint8_t const calculated_error_count = std::clamp(error.get() * query_size, 0.0, 255.0);
                    return calculated_error_count;
                }
            },
            error_variant);
    }
};

} // namespace seqan3::detail
