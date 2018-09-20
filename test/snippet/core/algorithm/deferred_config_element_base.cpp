#include <seqan3/core/algorithm/deferred_config_element_base.hpp>

using namespace seqan3;

int main() {}
//! [usage]
template <size_t I>
struct my_config
{
    size_t value{I};  // Has to be named `value`.
};

struct my_deferred_config
{
    template <typename fn_t, typename configuration_t>
    constexpr auto invoke(fn_t && fn, configuration_t && config) const
    requires detail::is_algorithm_configuration_v<remove_cvref_t<configuration_t>>
    {
        if (value == 0)
            return fn(std::forward<configuration_t>(config).replace_with(my_deferred_config{}, my_config<0>{}));
        else
            return fn(std::forward<configuration_t>(config).replace_with(my_deferred_config{}, my_config<1>{}));
    }

    int value{0};  // Has to be named `value`.
};
//! [usage]
