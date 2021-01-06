#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/simd/all.hpp>

// The simd type with 8 unsigned shorts.
using uint16x8_t = seqan3::simd_type_t<uint16_t, 8>;

int main()
{
    // Generate ascending simd index.
    for (auto && simd_id : seqan3::views::iota_simd<uint16x8_t>(0, 10))
        seqan3::debug_stream << simd_id << '\n'; // [0, 0, ..., 0], [1, 1, ..., 1], ... [9, 9, ..., 9]

    return 0;
}
