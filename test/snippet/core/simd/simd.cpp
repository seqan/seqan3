#include <seqan3/core/simd/all.hpp>

using namespace seqan3;

using uint8x16_t = simd_type_t<uint16_t, 8>;

int main()
{
    uint8x16_t a = fill<uint8x16_t>(4);
    uint8x16_t b = simd::fill<uint8x16_t>(5); // you can also explicitly use simd::
    uint8x16_t c = a + b;

    debug_stream << c << "\n"; // [9,9,9,9,9,9,9,9]
    return 0;
}
