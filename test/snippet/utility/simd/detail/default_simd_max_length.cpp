#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/simd/all.hpp>

int main()
{
    constexpr auto max_length = seqan3::detail::default_simd_max_length<seqan3::detail::default_simd_backend>;
    using uint8_simd_t = seqan3::simd::simd_type_t<uint8_t, max_length == 0 ? 1 : max_length>;

    uint8_simd_t a = seqan3::fill<uint8_simd_t>(4);
    uint8_simd_t b = seqan3::fill<uint8_simd_t>(5);
    uint8_simd_t c = a + b;

    seqan3::debug_stream << c << "\n"; // [9,9,9,9,9,9,9,9]
    return 0;
}
