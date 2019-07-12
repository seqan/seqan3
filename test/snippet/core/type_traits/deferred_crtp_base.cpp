#include <string>
#include <type_traits>

#include <seqan3/core/type_traits/deferred_crtp_base.hpp>

// Defines a crtp_base class with an additional value type.
template <typename derived_t, typename value_t>
class base1
{
public:

    value_t func1() const
    {
        return {"instance of base1"};
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
template <typename ...deferred_bases_t>
class derived : public seqan3::detail::invoke_deferred_crtp_base<deferred_bases_t, derived<deferred_bases_t...>>...
{};

int main()
{
    // Wraps the actual classes into deferred crtp base classes.
    // This step is necessary, since during their declaration the type of the derived class cannot be known.
    using deferred_base1 = seqan3::detail::deferred_crtp_base<base1, std::string>;
    using deferred_base2 = seqan3::detail::deferred_crtp_base<base2, uint8_t, uint32_t>;

    // Instantiate the derived class with the deferred crtp base classes.
    derived<deferred_base1, deferred_base2> d{};
    // Check the inherited interfaces.
    static_assert(std::is_same_v<decltype(d.func1()), std::string>, "Return type must be std::string");
    static_assert(std::is_same_v<decltype(d.func2(10u)), uint8_t>, "Return type must be uint8_t");
}
