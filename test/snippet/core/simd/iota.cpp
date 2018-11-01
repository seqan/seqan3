#include <seqan3/core/simd/all.hpp>

using namespace seqan3;

using uint8x16_t = simd_type_t<uint16_t, 8>;

int main()
{
    uint8x16_t a = iota<uint8x16_t>(1);
    debug_stream << a << "\n"; // [1,2,3,4,5,6,7,8]

    // or:

    uint8x16_t b = simd::iota<uint8x16_t>(1);
    debug_stream << b << "\n"; // [1,2,3,4,5,6,7,8]
    return 0;
}
