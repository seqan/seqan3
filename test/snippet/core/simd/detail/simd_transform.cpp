#include <seqan3/core/simd/simd.hpp>
#include <seqan3/core/simd/simd_algorithm.hpp>

int main()
{
{
//! [definition]
using simd_type = seqan3::simd::simd_type_t<uint16_t, 8>;
using scalar_t = typename seqan3::simd::simd_traits<simd_type>::scalar_type; // same as uint16_t

auto fn = [] (size_t, scalar_t ai, scalar_t bi, /*scalar_t ci, ...,*/ scalar_t zi) -> scalar_t
{
    return ai + bi + /*ci + ... +*/ zi;
};

simd_type a{}, b{}, /*c{}, ...,*/ z{};
simd_type result = seqan3::detail::simd_transform<simd_type>(fn, a, b, /*c{}, ...,*/ z);

// same as
result = simd_type
{
    fn(0, a[0], b[0], /*c[0], ...,*/ z[0]),
    fn(1, a[1], b[1], /*c[1], ...,*/ z[1]),
    fn(2, a[2], b[2], /*c[2], ...,*/ z[2]),
    // ...
    fn(7, a[7], b[7], /*c[7], ...,*/ z[7])
};
result = simd_type{};

for (int i = 0; i < 7; ++i)
{
    result[i] = fn(i, a[i], b[i], /*c[i], ...,*/ z[i]);
}
//! [definition]
    static_cast<void>(result);
}
{
//! [generator]
using simd_type = seqan3::simd::simd_type_t<uint16_t, 8>;
using scalar_t = typename seqan3::simd::simd_traits<simd_type>::scalar_type; // same as uint16_t

auto iota = [i = scalar_t{0}] (size_t) mutable
{
    return i++;
};

simd_type result = seqan3::detail::simd_transform<simd_type>(iota);

// same as
result = simd_type
{
    iota(0),
    iota(1),
    iota(2),
    // ...
    iota(7)
};
//! [generator]
    static_cast<void>(result);
}
{
//! [binary_max]
using simd_type = seqan3::simd::simd_type_t<uint16_t, 8>;

auto max = [] (size_t, auto ai, auto bi)
{
    return ai > bi ? ai : bi;
};

simd_type a{}, b{};
simd_type result = seqan3::detail::simd_transform<simd_type>(max, a, b);

// same as
result = simd_type
{
    max(0, a[0], b[0]),
    max(1, a[1], b[1]),
    max(2, a[2], b[2]),
    // ...
    max(7, a[7], b[7])
};
//! [binary_max]
    static_cast<void>(result);
}
    return 0;
}
