#include <seqan3/core/debug_stream.hpp>
#include <seqan3/core/simd/all.hpp>

using namespace seqan3;

using uint16x8_t = simd_type_t<uint16_t, 8>;

int main()
{
    uint16x8_t a = fill<uint16x8_t>(4);
    debug_stream << a << "\n"; // [4,4,4,4,4,4,4,4]

    // or:

    uint16x8_t b = simd::fill<uint16x8_t>(4);
    debug_stream << b << "\n"; // [4,4,4,4,4,4,4,4]
    return 0;
}
