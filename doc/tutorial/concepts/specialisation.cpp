#include <utility>                          // for std::pair
#include <seqan3/std/concepts>

template <typename t>
struct square_root_type;

template <std::Integral t>
struct square_root_type<t>
{
    using type = std::pair<float, float>;   // real and imaginary part
};

template <std::UnsignedIntegral t>
struct square_root_type<t>
{
    using type = float;                     // doesn't need imaginary part
};

// `int` models std::Integral but not std::UnsignedIntegral:
static_assert(std::Same<typename square_root_type<int>::type,           std::pair<float, float>>);

// `unsigned` models std::Integral and std::UnsignedIntegral, but the latter is more refined:
static_assert(std::Same<typename square_root_type<unsigned>::type,      float>);
