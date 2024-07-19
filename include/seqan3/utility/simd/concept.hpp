// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::simd::simd_concept.
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 */

#pragma once

#include <concepts>
#include <type_traits>

#include <seqan3/utility/simd/simd_traits.hpp>

namespace seqan3::detail
{
//!\cond
template <template <typename> typename rebind>
struct simd_traits_has_rebind : std::true_type
{};
// NOTE: this definition should be used for seqan3::simd, but gcc has a bug that it will not fail silently if
// simd_t is a pointer to a incomplete type. Furthermore the is_pointer_v should prevent those cases by checking the type
// beforehand, but for some reasons the short-circuit semantic of `&&` does not work in this case and gcc still evaluates
// the requires clause which in turn triggers the error.
//
// If this concept is used directly on incomplete types it will produces this compiler error:
//     error: invalid use of incomplete type ‘struct incomplete::template_type<int>’
//          requires std::same_as<decltype(a - b), simd_t>;
template <typename simd_t>
concept simd_concept = requires (simd_t a, simd_t b) {
    typename simd_traits<std::remove_reference_t<simd_t>>::scalar_type;
    typename simd_traits<std::remove_reference_t<simd_t>>::mask_type;
    typename simd_traits<std::remove_reference_t<simd_t>>::swizzle_type;
    requires simd_traits_has_rebind<simd_traits<std::remove_reference_t<simd_t>>::template rebind>::value;

    // require that static member variables are defined
    requires std::integral<decltype(simd_traits<std::remove_reference_t<simd_t>>::length)>;
    requires std::integral<decltype(simd_traits<std::remove_reference_t<simd_t>>::max_length)>;

    // assume array access that returns a scalar_type type
    { a[0] } -> std::convertible_to<typename simd_traits<std::remove_reference_t<simd_t>>::scalar_type>;

    // require comparison operators
    { a == b } -> std::same_as<typename simd_traits<std::remove_reference_t<simd_t>>::mask_type>;
    { a != b } -> std::same_as<typename simd_traits<std::remove_reference_t<simd_t>>::mask_type>;
    { a < b } -> std::same_as<typename simd_traits<std::remove_reference_t<simd_t>>::mask_type>;
    { a > b } -> std::same_as<typename simd_traits<std::remove_reference_t<simd_t>>::mask_type>;
    { a <= b } -> std::same_as<typename simd_traits<std::remove_reference_t<simd_t>>::mask_type>;
    { a >= b } -> std::same_as<typename simd_traits<std::remove_reference_t<simd_t>>::mask_type>;

    // require arithmetic operators
    { a + b } -> std::same_as<std::remove_reference_t<simd_t>>;
    { a - b } -> std::same_as<std::remove_reference_t<simd_t>>;
    { a * b } -> std::same_as<std::remove_reference_t<simd_t>>;
    { a / b } -> std::same_as<std::remove_reference_t<simd_t>>;
    { a += b } -> std::same_as<std::remove_reference_t<simd_t> &>;
    { a -= b } -> std::same_as<std::remove_reference_t<simd_t> &>;
    { a *= b } -> std::same_as<std::remove_reference_t<simd_t> &>;
    { a /= b } -> std::same_as<std::remove_reference_t<simd_t> &>;
};
//!\endcond

} // namespace seqan3::detail

namespace seqan3
{

inline namespace simd
{

/*!\interface seqan3::simd::simd_concept <>
 * \brief The generic simd concept.
 * \ingroup utility_simd
 *
 * \details
 *
 * seqan3::simd::simd_concept checks whether a given type is a simd type. One of the prerequisites is
 * that seqan3::simd::simd_traits is defined for this type.
 *
 * \if DEV
 * \todo Simplify concept to the seqan3::detail::simd_conceptonce gcc bug is fixed
 * \endif
 */
//!\cond
template <typename simd_t>
concept simd_concept = !std::is_pointer_v<std::decay_t<simd_t>> && detail::simd_concept<simd_t>;
//!\endcond

/*!\interface seqan3::simd::simd_index <>
 * \brief Refines the seqan3::simd::simd_concept requiring the underlying scalar type to model std::integral.
 * \ingroup utility_simd
 *
 * \details
 *
 * seqan3::simd::simd_index checks whether a given type is a simd type and can be used to represent an index.
 */
//!\cond
template <typename t>
concept simd_index = simd::simd_concept<t> &&
    requires ()
{
    requires std::integral<typename simd_traits<std::remove_reference_t<t>>::scalar_type>;
};
//!\endcond

} // namespace simd

} // namespace seqan3
