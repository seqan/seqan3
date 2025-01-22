// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <type_traits>

#include <seqan3/alignment/aligned_sequence/aligned_sequence_concept.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/core/detail/is_class_template_declarable.hpp>

template <typename t>
    requires std::is_integral_v<t>
struct foo
{
    t value;
};

// foo is declarable with int, i.e. foo<int> is a valid expression
static_assert(seqan3::detail::is_class_template_declarable_with_v<foo, int>);
// foo is not declarable with double, because it does not fulfil the requires clause of foo.
static_assert(!seqan3::detail::is_class_template_declarable_with_v<foo, double>);

// This also works with std::enable_if and producing a substitution failure.
template <typename t, typename = std::enable_if_t<std::is_same_v<t, int>>>
struct bar
{
    t value;
};

// bar is declarable with int, i.e. bar<int> is a valid expression
static_assert(seqan3::detail::is_class_template_declarable_with_v<bar, int>);
// bar is not declarable with double, because it produces an substitution failure (SFINAE).
static_assert(!seqan3::detail::is_class_template_declarable_with_v<bar, double>);

// is_class_template_declarable_with_v works well with lazy_conditional_t
template <typename t>
using maybe_foo_t = seqan3::detail::
    lazy_conditional_t<seqan3::detail::is_class_template_declarable_with_v<foo, t>, seqan3::detail::lazy<foo, t>, t>;

int main()
{
    foo<int> a = maybe_foo_t<int>{10};                // foo is instantiable with int, thus use foo<int>
    seqan3::debug_stream << "a: " << a.value << '\n'; // prints 10
    float b = maybe_foo_t<float>{0.4f};               // foo is not instantiable with float, thus use float directly
    seqan3::debug_stream << "b: " << b << '\n';       // prints 0.4
    return 0;
}
