#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>

enum struct my_id : int
{
    bar_id,
    foo_id
};

struct bar : public seqan3::pipeable_config_element<bar, float>
{
    static constexpr my_id id{my_id::bar_id};
};

template <typename t>
struct foo : public seqan3::pipeable_config_element<foo<t>, t>
{
    static constexpr my_id id{my_id::foo_id};
};

template <typename t>
foo(t) -> foo<t>;
