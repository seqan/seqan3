#include <seqan3/core/simd/all.hpp>

using namespace seqan3;

using int16x8_t = simd_type_t<int16_t, 8>;
using int32x4_t = simd_type_t<int32_t, 4>;

int main()
{
    int16x8_t a{0, -1, 2, -3, 4, -5, 6, -7};
    int32x4_t b = upcast<int32x4_t>(a);

    debug_stream << b << "\n"; // [0,-1,2,-3]
    return 0;
}
