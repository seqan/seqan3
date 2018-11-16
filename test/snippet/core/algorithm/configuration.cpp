#include <iostream>

//! [configuration_setup]
#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>

using namespace seqan3;

enum struct my_id : int
{
    bar_id,
    foo_id
};

struct bar : public pipeable_config_element
{
    static constexpr my_id id{my_id::bar_id};

    bar()                        = default;
    bar(bar const &)             = default;
    bar(bar &&)                  = default;
    bar & operator=(bar const &) = default;
    bar & operator=(bar &&)      = default;
    ~bar()                       = default;

    bar(float v) : value(v)
    {}

    float value;
};

template <typename t>
struct foo : public pipeable_config_element
{
    static constexpr my_id id{my_id::foo_id};

    foo()                        = default;
    foo(foo const &)             = default;
    foo(foo &&)                  = default;
    foo & operator=(foo const &) = default;
    foo & operator=(foo &&)      = default;
    ~foo()                       = default;

    foo(t v) : value(v)
    {}

    t value;
};
//! [configuration_setup]

//! [compatibility]
namespace seqan3::detail
{
template <>
inline constexpr std::array<std::array<int, 2>, 2> compatibility_table<my_id>
{
    {
        {0, 1},
        {1, 0}
    }
};
} // namespace seqan3::detail
//! [compatibility]
int main()
{
{
//! [combine]
configuration my_cfg = bar{1.3} | foo<int>{{4}};  // my_cfg is now of type configuration<bar, foo<int>>
//! [combine]

//! [get]
std::cout << get<1>(my_cfg).value   << '\n';  // prints 4
std::cout << get<bar>(my_cfg).value << '\n';  // prints 1.3
std::cout << get<foo>(my_cfg).value << '\n';  // prints 4
//! [get]
}

{
//! [push_back]
configuration my_cfg = configuration{foo<int>{4}}.push_back(bar{});
//! [push_back]
}

{
//! [value_or]
configuration my_cfg{bar{1.3}};
std::cout << my_cfg.value_or<bar>("not there") << '\n';  // prints: 1.3
std::cout << my_cfg.value_or<foo>("not there") << '\n';  // prints: not there
//! [value_or]
}
{
//! [constructor]
configuration cfg{bar{1.2}};
//! [constructor]
}
}
