#include <seqan3/core/simd/all.hpp>

using namespace seqan3;

using uint8x16_t = simd_type_t<uint16_t, 8>;

int main()
{
    uint8x16_t a = fill<uint8x16_t>(4);
    debug_stream << a << "\n"; // [4,4,4,4,4,4,4,4]

    //or:

    uint8x16_t b = simd::fill<uint8x16_t>(4);
    debug_stream << b << "\n"; // [4,4,4,4,4,4,4,4]
    return 0;
}
