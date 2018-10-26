#include <seqan3/core/simd/all.hpp>

using namespace seqan3;
using namespace seqan3::simd;

int main()
{
    constexpr auto max_length = detail::default_simd_max_length<detail::default_simd_backend>;
    using uint8_simd_t = simd_type_t<uint8_t, max_length>;

    uint8_simd_t a = fill<uint8_simd_t>(4);
    uint8_simd_t b = fill<uint8_simd_t>(5);
    uint8_simd_t c = a + b;

    debug_stream << c << "\n"; // (9,9,9,9,9,9,9,9)
    return 0;
}
