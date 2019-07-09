#include <seqan3/core/debug_stream.hpp>
#include <seqan3/core/simd/all.hpp>

using namespace seqan3;

using uint16x8_t = simd_type_t<uint16_t, 8>;

int main()
{
    uint16x8_t a = fill<uint16x8_t>(4);
    uint16x8_t b = simd::fill<uint16x8_t>(5); // you can also explicitly use simd::
    uint16x8_t c = a + b;

    debug_stream << c << "\n"; // [9,9,9,9,9,9,9,9]
    return 0;
}
