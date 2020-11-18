#include <seqan3/core/algorithm/pipeable_config_element.hpp>
#include <seqan3/core/configuration/configuration.hpp>

enum struct my_id : int
{
    bar_id,
    foo_id
};

class bar : public seqan3::pipeable_config_element<bar>
{
public:

    bar() = default;
    bar(bar const &) = default;
    bar(bar &&) = default;
    bar & operator=(bar const &) = default;
    bar & operator=(bar &&) = default;
    ~bar() = default;

    static constexpr my_id id{my_id::bar_id};
};

template <typename t>
struct foo : public seqan3::pipeable_config_element<foo<t>>
{
public:

    foo() = default;
    foo(foo const &) = default;
    foo(foo &&) = default;
    foo & operator=(foo const &) = default;
    foo & operator=(foo &&) = default;
    ~foo() = default;

    static constexpr my_id id{my_id::foo_id};
};

template <typename t>
foo(t) -> foo<t>;
