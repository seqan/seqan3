#include <seqan3/utility/simd/detail/builtin_simd.hpp>

// 8x 16bit integers = 128bit
using int16x8_t = seqan3::detail::builtin_simd<int16_t, 8u>::type; // sse4
