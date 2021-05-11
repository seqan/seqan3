#include <seqan3/core/configuration/configuration.hpp>
#include <seqan3/core/configuration/pipeable_config_element.hpp>

enum struct my_id : int
{
    bar_id,
    foo_id
};

class bar : private seqan3::pipeable_config_element
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
struct foo : private seqan3::pipeable_config_element
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
