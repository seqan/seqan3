// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <string>
#include <type_traits>

#include <seqan3/core/detail/deferred_crtp_base.hpp>

// Defines a crtp_base class with an additional value type.
template <typename derived_t, int value>
class base1
{
public:
    int func1() const
    {
        return value;
    }
};

// Defines a crtp_base class with an additional value type and a parameter type.
template <typename derived_t, typename value_t, typename parameter_t>
class base2
{
public:
    value_t func2(parameter_t const p) const
    {
        return static_cast<value_t>(p);
    }
};

// The derived class that inherits from a variadic crtp pattern, which are augmented with additional trait types.
// These types must be wrapped in a deferred layer, otherwise the compilation fails as incomplete types are not allowed.
// But during the definition of the base classes, the derived class cannot be known.
// In addition the deferred type must be invoked with the derived class using the `invoke_deferred_crtp_base` helper
// template to instantiate the correct crtp base type.
// Note that it is possible to define base classes with type template parameters (see base2) or
// non-type template parameters (see base1), but non-type and type template parameters cannot be mixed in one
// base class definition.
template <typename... deferred_bases_t>
class derived : public seqan3::detail::invoke_deferred_crtp_base<deferred_bases_t, derived<deferred_bases_t...>>...
{};

int main()
{
    // Define deferred base with non-type template parameter
    using deferred_base1 = seqan3::detail::deferred_crtp_base_vargs<base1, 10>;
    // Define deferred base with type template parameter.
    using deferred_base2 = seqan3::detail::deferred_crtp_base<base2, uint8_t, uint32_t>;

    // Instantiate the derived class with the deferred crtp base classes.
    derived<deferred_base1, deferred_base2> d{};
    // Check the inherited interfaces.
    static_assert(std::is_same_v<decltype(d.func1()), int>, "Return type must be int");
    static_assert(std::is_same_v<decltype(d.func2(10u)), uint8_t>, "Return type must be uint8_t");
}
