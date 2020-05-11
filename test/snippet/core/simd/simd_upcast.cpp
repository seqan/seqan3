#include <seqan3/core/debug_stream.hpp>
#include <seqan3/core/simd/all.hpp>

using int16x8_t = seqan3::simd::simd_type_t<int16_t, 8>;
using int32x4_t = seqan3::simd::simd_type_t<int32_t, 4>;

int main()
{
    int16x8_t a{0, -1, 2, -3, 4, -5, 6, -7};
    int32x4_t b = seqan3::simd::upcast<int32x4_t>(a);

    seqan3::debug_stream << b << '\n'; // [0,-1,2,-3]
    return 0;
}
