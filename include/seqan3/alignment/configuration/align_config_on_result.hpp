// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::align_cfg::on_result.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <type_traits>

#include <seqan3/alignment/configuration/detail.hpp>
#include <seqan3/contrib/std/detail/movable_box.hpp>
#include <seqan3/core/configuration/pipeable_config_element.hpp>

namespace seqan3::align_cfg
{

/*!\brief Configuration element to provide a user defined callback function for the alignment.
 * \ingroup alignment_configuration
 *
 * \tparam callback_t The type of the callback; must model std::invocable with the generated seqan3::alignment_result
 *                    and std::move_constructible.
 *
 * \details
 *
 * Allows the user to specify a callback that should be called for every computed alignment result. The callback
 * must take exactly one argument for the alignment result and return `void`. If the user callback is
 * specified, the call to the alignment algorithm seqan3::align_pairwise will return nothing, i.e. it does not return
 * a seqan3::algorithm_result_generator_range anymore. Note that within a parallel configuration the order of the generated alignment
 * results and therefore the call to the user callback is non-deterministic. However, the continuation interface with the
 * user callback can be more efficient in a concurrent environment. If you pass an lvalue function object as callback
 * function, you need to make sure that the referenced function object outlives the call to the alignment algorithm.
 *
 * \if DEV
 * The given callback is wrapped inside a seqan::stl::detail::movable_box wrapper type. This allows to also
 * use lambdas with a capture block, which otherwise are not std::copy_assignable and therefore invalidate the
 * requirements for the configuration element (must model std::semiregular).
 * \endif
 *
 * ### Example
 *
 * The following code snippet demonstrates the basic usage:
 *
 * \include test/snippet/alignment/configuration/align_cfg_on_result.cpp
 */
template <std::move_constructible callback_t>
class on_result : private seqan3::pipeable_config_element
{
public:
    //!\brief The stored callable which will be invoked with the alignment result.
    seqan::stl::detail::movable_box_t<callback_t> callback; // Allows lambdas with capture blocks.

    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr on_result() = default;                              //!< Defaulted.
    constexpr on_result(on_result const &) = default;             //!< Defaulted.
    constexpr on_result(on_result &&) = default;                  //!< Defaulted.
    constexpr on_result & operator=(on_result const &) = default; //!< Defaulted.
    constexpr on_result & operator=(on_result &&) = default;      //!< Defaulted.
    ~on_result() = default;                                       //!< Defaulted.

    /*!\brief Constructs the configuration element with the given user callback.
     * \param[in] callback The callback to invoke for a computed seqan3::alignment_result.
     */
    constexpr explicit on_result(callback_t && callback) : callback{std::forward<callback_t>(callback)}
    {}
    //!\}

    //!\privatesection
    //!\brief Internal id to check for consistent configuration settings.
    static constexpr seqan3::detail::align_config_id id{seqan3::detail::align_config_id::on_result};
};

/*!\name Type deduction guides
 * \{
 */
//!\brief Deduces the callback type from a forwarding constructor argument.
template <std::move_constructible callback_t>
on_result(callback_t &&) -> on_result<std::decay_t<callback_t>>;
//!\}
} // namespace seqan3::align_cfg
