#include <type_traits>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/core/detail/pack_algorithm.hpp>

// With c++20 you could also write it like this
// auto fn = []<typename value_t>(value_t && value)
// {
// ...
// };
auto fn = [](auto value)
{
    // id is the original type not wrapped in std::type_identity.
    using value_t = decltype(value);

    if constexpr(std::is_same_v<value_t, bool>)
        return value == false;
    else if constexpr(std::is_same_v<value_t, int>)
        return value == 3;
    else if constexpr(std::is_same_v<value_t, double>)
        return std::abs(value - 1.2) < 0.00001;
    else
        return false;
};

static_assert(seqan3::detail::all_of(fn, 3, 1.2, false)); // evalates to true
static_assert(!seqan3::detail::all_of(fn, 3, 1.2, false, "something else")); // evalates to false
