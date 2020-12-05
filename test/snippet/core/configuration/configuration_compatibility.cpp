#include <seqan3/core/configuration/configuration.hpp>

enum struct my_id : int
{
    bar_id,
    foo_id
};

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
