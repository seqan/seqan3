#include <seqan3/core/algorithm/configuration.hpp>

using namespace seqan3;

struct bar
{
    float value = 0;
};

struct foo
{
    int value = 0;
};

struct with_foo_adaptor : public detail::configuration_fn_base<with_foo_adaptor>
{
    template <typename configuration_t>
    constexpr auto invoke(configuration_t && cfg,
                          int new_v) const
        requires detail::is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
    {
        return std::forward<configuration_t>(cfg).push_front(foo{new_v});
    }

    template <typename configuration_t>
    constexpr auto invoke(configuration_t && cfg) const
        requires detail::is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
    {
        return std::forward<configuration_t>(cfg).push_front(foo{2});
    }
};

struct with_bar_adaptor : public detail::configuration_fn_base<with_bar_adaptor>
{
    template <typename configuration_t>
    constexpr auto invoke(configuration_t && cfg,
                          float new_v) const
        requires detail::is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
    {
        return std::forward<configuration_t>(cfg).push_front(bar{new_v});
    }

    template <typename configuration_t>
    constexpr auto invoke(configuration_t && cfg) const
        requires detail::is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
    {
        return std::forward<configuration_t>(cfg).push_front(bar{2});
    }
};

inline constexpr with_foo_adaptor with_foo;
inline constexpr with_bar_adaptor with_bar;

int main()
{

//! [combine]
auto my_cfg = detail::configuration<bar>{} | with_foo(1);  // my_cfg is now of type configuration<foo, bar>
//! [combine]

[[maybe_unused]]
//! [access]
auto bar_value = std::get<1>(my_cfg).value;
//! [access]

{
//! [constructor]
detail::configuration cfg = with_foo(1);
detail::configuration cfg1 = with_foo;
//! [constructor]
}

{
//! [combine_2]
// case 1: adaptor with adaptor
auto cfg1 = with_foo | with_bar;

// case 2: adaptor with proxy
auto cfg2 = with_foo | with_bar(2);

// case 3: proxy with adaptor
auto cfg3 = with_foo(1) | with_bar;

// case 4: proxy with proxy
auto cfg4 = with_foo(1) | with_bar(2);
//! [combine_2]
}
}
