#include <utility>                          // for std::pair
#include <seqan3/std/concepts>

template <typename t>
struct square_root_type;

template <std::integral t>
struct square_root_type<t>
{
    using type = std::pair<float, float>;   // real and imaginary part
};

template <std::unsigned_integral t>
struct square_root_type<t>
{
    using type = float;                     // doesn't need imaginary part
};

// `int` models std::integral but not std::unsigned_integral:
static_assert(std::same_as<typename square_root_type<int>::type,           std::pair<float, float>>);

// `unsigned` models std::integral and std::unsigned_integral, but the latter is more refined:
static_assert(std::same_as<typename square_root_type<unsigned>::type,      float>);
