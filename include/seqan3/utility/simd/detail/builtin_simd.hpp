// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \author Marcel Ehrhardt <marcel.ehrhardt AT fu-berlin.de>
 * \brief Provides seqan3::detail::builtin_simd, seqan3::detail::is_builtin_simd
 *        and seqan3::simd::simd_traits<builtin_simd_t>.
 */

#pragma once

#include <bit>
#include <type_traits>

#include <seqan3/utility/detail/integer_traits.hpp>
#include <seqan3/utility/simd/detail/default_simd_length.hpp>
#include <seqan3/utility/simd/simd_traits.hpp>
#include <seqan3/utility/type_traits/detail/transformation_trait_or.hpp>

namespace seqan3::detail
{

/*!\brief A class that holds the type of a simd implementation called [vector extension]
 * (https://gcc.gnu.org/onlinedocs/gcc/Vector-Extensions.html)
 * (formerly known as "seqan simd" in seqan2).
 * \ingroup utility_simd
 * \tparam scalar_t The underlying type of a simd vector
 * \tparam length   The number of packed values in a simd vector
 *
 * \include test/snippet/utility/simd/detail/builtin_simd.cpp
 *
 * seqan3::detail::builtin_simd is basically defined as:
 *
 * ```
 * template <typename scalar_t, size_t length>
 * struct builtin_simd
 * {
 *    using type [[gnu::vector_size(sizeof(scalar_t) * length)]] = scalar_t;
 * };
 * ```
 *
 * \attention This class itself only delegates to a [vector extension]
 * (https://gcc.gnu.org/onlinedocs/gcc/Vector-Extensions.html) type,
 * which is offered by the compiler as a builtin type.
 *
 * \sa https://gcc.gnu.org/onlinedocs/gcc/Vector-Extensions.html
 */
template <typename scalar_t, size_t length>
struct builtin_simd;

//!\copydoc seqan3::detail::builtin_simd
//!\ingroup utility_simd
template <typename scalar_t, size_t length>
    requires (std::has_single_bit(length))
struct builtin_simd<scalar_t, length>
{
    //!\brief The type of the builtin simd.
#if SEQAN3_DOXYGEN_ONLY(1) 0
    using type = scalar_t __attribute__((vector_size(sizeof(scalar_t) * length))));
    // doxygen 1.8.13 does not support c++11 attributes, thus this doxygen-only definition
#elif defined(__clang__)
    using type = scalar_t __attribute__((ext_vector_type(length)));
#else
    using type [[gnu::vector_size(sizeof(scalar_t) * length)]] = scalar_t;
#endif
};

/*!\brief Helper struct for seqan3::detail::is_builtin_simd
 * \ingroup utility_simd
 * \sa seqan3::detail::is_builtin_simd
 */
template <typename builtin_simd_t>
struct builtin_simd_traits_helper : std::false_type
{};

/*!\brief Helper struct for seqan3::detail::is_builtin_simd
 * \ingroup utility_simd
 * \sa seqan3::detail::is_builtin_simd
 */
template <typename builtin_simd_t>
//!\cond
// NOTE: gcc throws a compile time error if builtin_simd_t is a pointer of an incomplete type. To tackle this we
// short-circuit the requires with is_pointer_v. See builtin_simd_test.cpp for a test case for this.
    requires (!std::is_pointer_v<std::decay_t<builtin_simd_t>>) && requires (builtin_simd_t simd) {
        { simd[0] };
    }
//!\endcond
struct builtin_simd_traits_helper<builtin_simd_t>
{
    //!\brief The scalar type of builtin_simd_t
    using scalar_type = std::remove_reference_t<decltype(std::declval<builtin_simd_t>()[0])>;
    //!\brief The length of builtin_simd_t
    static constexpr auto length = sizeof(builtin_simd_t) / sizeof(scalar_type);

    //!\brief Whether builtin_simd_t is a builtin type or not?
    //!\hideinitializer
    static constexpr bool value =
        std::has_single_bit(length)
        && std::is_same_v<builtin_simd_t, transformation_trait_or_t<builtin_simd<scalar_type, length>, void>>;
};

/*!\brief This class inherits from std::true_type, **iff**
 * seqan3::detail::builtin_simd<scalar_t, length>::type is a builtin simd type.
 * \ingroup utility_simd
 * \tparam builtin_simd_t The type to check.
 *
 * \include test/snippet/utility/simd/detail/is_builtin_simd.cpp
 * \sa https://gcc.gnu.org/onlinedocs/gcc/Vector-Extensions.html
 */
template <typename builtin_simd_t>
struct is_builtin_simd : std::bool_constant<builtin_simd_traits_helper<builtin_simd_t>::value>
{};

/*!\brief Helper variable to test whether a type is a simd builtin type.
 * \ingroup utility_simd
 * \tparam builtin_simd_t The type to check.
 * \see seqan3::detail::is_builtin_simd
 */
template <typename builtin_simd_t>
constexpr bool is_builtin_simd_v = is_builtin_simd<builtin_simd_t>::value;

/*!\brief This function specializes seqan3::detail::default_simd_max_length for
 * seqan3::detail::builtin_simd.
 * \ingroup utility_simd
 *
 * The redefinition of *default_simd_max_length* influences the default
 * *length* (i.e., seqan3::detail::default_simd_length) of seqan3::simd::simd_type for
 * seqan3::detail::builtin_simd types.
 */
template <>
inline constexpr auto default_simd_max_length<builtin_simd> = []()
{
#if defined(__AVX512F__)
    return min_viable_uint_v<64u>;
#elif defined(__AVX2__)
    return min_viable_uint_v<32u>;
#elif defined(__SSE4_1__) && defined(__SSE4_2__)
    return min_viable_uint_v<16u>;
#else
    return min_viable_uint_v<0u>;
#endif
}();

/*!\brief This class inherits from std::true_type, **iff** the builtin simd type is supported by the current
 *        architecture.
 * \ingroup utility_simd
 * \tparam builtin_simd_t The type to check.
 *
 * \details
 *
 * A builtin simd type is native if the following conditions are true:
 * * the default simd max length is not equal to `0`.
 * * the max length of the simd type is at least 16 (SSE4)
 * * the max length of the simd type is at most 64 (AVX512)
 */
template <typename builtin_simd_t>
struct is_native_builtin_simd :
    std::bool_constant<(default_simd_max_length<builtin_simd> != 0)
                       && ((builtin_simd_traits_helper<builtin_simd_t>::length
                            * sizeof(typename builtin_simd_traits_helper<builtin_simd_t>::scalar_type))
                           >= 16)
                       && ((builtin_simd_traits_helper<builtin_simd_t>::length
                            * sizeof(typename builtin_simd_traits_helper<builtin_simd_t>::scalar_type))
                           <= 64)>
{};

/*!\brief Helper variable to test whether a type is a native simd builtin type.
 * \ingroup utility_simd
 * \tparam builtin_simd_t The type to check.
 * \see seqan3::detail::is_native_builtin_simd_v
 */
template <typename builtin_simd_t>
constexpr bool is_native_builtin_simd_v = is_native_builtin_simd<builtin_simd_t>::value;

} // namespace seqan3::detail

namespace seqan3
{

inline namespace simd
{

/*!\brief This class specializes seqan3::simd::simd_traits for seqan3::detail::builtin_simd types
 * \tparam builtin_simd_t A simd type that satisfies seqan3::detail::is_builtin_simd_v<builtin_simd_t>.
 * \ingroup utility_simd
 * \sa seqan3::simd::simd_traits for more information
 */
template <typename builtin_simd_t>
// \cond
    requires detail::is_builtin_simd<builtin_simd_t>::value
// \endcond
struct simd_traits<builtin_simd_t>
{
    //!\copydoc seqan3::simd::simd_traits::scalar_type
    using scalar_type = typename detail::builtin_simd_traits_helper<builtin_simd_t>::scalar_type;
    //!\copydoc seqan3::simd::simd_traits::length
    static constexpr auto length = detail::builtin_simd_traits_helper<builtin_simd_t>::length;
    //!\copydoc seqan3::simd::simd_traits::max_length
    static constexpr auto max_length = length == 1u ? length : sizeof(scalar_type) * length;

    static_assert(std::is_integral_v<scalar_type>, "For now we assume that builtin simd can only be integers");
    //!\copydoc seqan3::simd::simd_traits::mask_type
    using mask_type = decltype(std::declval<builtin_simd_t>() == std::declval<builtin_simd_t>());
    //!\copydoc seqan3::simd::simd_traits::swizzle_type
    using swizzle_type = typename detail::builtin_simd<uint8_t, max_length>::type;

    //!\copydoc seqan3::simd::simd_traits::rebind
    template <typename new_scalar_type>
    // \cond
        requires (sizeof(scalar_type) == sizeof(new_scalar_type))
    // \endcond
    using rebind = typename detail::builtin_simd<new_scalar_type, length>::type;
};

} // namespace simd

} // namespace seqan3
