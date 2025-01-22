// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Provides type traits seqan3::detail::transfer_type_modifier_onto.
 */

#pragma once

#include <type_traits>

#include <seqan3/core/platform.hpp>

namespace seqan3::detail
{

/*!\brief Transfers the type modifier `&`, `&&` and `const` (and any combination) to the target type.
 * \implements seqan3::transformation_trait
 * \tparam source_t The type you wish to transfer the type modifier from.
 * \tparam target_t The type you wish to transfer the type modifier to.
 * \ingroup core
 * If the `target_t` already has a type modifier, e.g. `const`, it will keep that type modifier.
 *
 * If the resulting type would have the type modifier `&&` and `&` at the same time, it follows the rule of reference
 * collapsing, that means `&` will be preferred.
 * \sa https://en.cppreference.com/w/cpp/language/reference#Reference_collapsing
 */
template <typename source_t, typename target_t>
struct transfer_type_modifier_onto
{
private:
    //!\brief Transfers the `const` type modifier to the target type.
    using maybe_const_target_t = std::conditional_t<std::is_const_v<std::remove_reference_t<source_t>>
                                                        || std::is_const_v<std::remove_reference_t<target_t>>,
                                                    std::add_const_t<std::remove_cvref_t<target_t>>,
                                                    std::remove_cvref_t<target_t>>;

    //!\brief Transfers the `&&` type modifier to the target type.
    using maybe_rvalue_reference_t =
        std::conditional_t<std::is_rvalue_reference_v<source_t> || std::is_rvalue_reference_v<target_t>,
                           std::add_rvalue_reference_t<maybe_const_target_t>,
                           maybe_const_target_t>;

    //!\brief Transfers the `&` type modifier to the target type.
    using maybe_lvalue_reference_target_t =
        std::conditional_t<std::is_lvalue_reference_v<source_t> || std::is_lvalue_reference_v<target_t>,
                           std::add_lvalue_reference_t<maybe_rvalue_reference_t>,
                           maybe_rvalue_reference_t>;

public:
    //!\brief Transfers the type modifier `&`, `&&` and `const` (and any combination) to the target type.
    using type = maybe_lvalue_reference_target_t;
};

/*!\brief Transfers the type modifier `&`, `&&` and `const` (and any combination) to the target type
 *        (transformation_trait shortcut).
 * \tparam source_t The type you wish to transfer the type modifier from.
 * \tparam target_t The type you wish to transfer the type modifier to.
 * \see seqan3::detail::transfer_type_modifier_onto
 * \relates seqan3::detail::transfer_type_modifier_onto
 */
template <typename source_t, typename target_t>
using transfer_type_modifier_onto_t = typename transfer_type_modifier_onto<source_t, target_t>::type;

} // namespace seqan3::detail
