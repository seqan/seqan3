#include <seqan3/core/algorithm/configuration.hpp>
#include <seqan3/core/algorithm/pipeable_config_element.hpp>
#include <seqan3/core/debug_stream.hpp>

// Initial setup used in the actual example:
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

int main()
{
    seqan3::configuration my_cfg{bar{1.3}};
    seqan3::debug_stream << my_cfg.value_or<bar>("not there") << '\n';  // prints: 1.3
    seqan3::debug_stream << my_cfg.value_or<foo>("not there") << '\n';  // prints: not there
}
