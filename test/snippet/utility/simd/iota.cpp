#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/simd/all.hpp>

using uint16x8_t = seqan3::simd::simd_type_t<uint16_t, 8>;

int main()
{
    uint16x8_t a = seqan3::iota<uint16x8_t>(1);
    seqan3::debug_stream << a << "\n"; // [1,2,3,4,5,6,7,8]

    // or:

    uint16x8_t b = seqan3::simd::iota<uint16x8_t>(1);
    seqan3::debug_stream << b << "\n"; // [1,2,3,4,5,6,7,8]
    return 0;
}
