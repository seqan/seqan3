// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <seqan3/utility/type_traits/detail/transformation_trait_or.hpp>

template <typename T>
struct A;

template <>
struct A<int>
{
    using type = int;
};

// A<unsigned>::type is not defined, thus falling back to `void`
static_assert(std::is_same_v<void, seqan3::detail::transformation_trait_or_t<A<unsigned>, void>>);

// A<int>::type is defined, use A<int>::type
static_assert(std::is_same_v<int, seqan3::detail::transformation_trait_or_t<A<int>, void>>);
