#include <seqan3/simd/detail/ume_simd.hpp>

#if __has_include(<umesimd/UMESimd.h>)

// 8x 16bit integers = 128bit
using int16x8_t = seqan3::detail::ume_simd<int16_t, 8u>::type; // sse4

#endif // __has_include(<umesimd/UMESimd.h>)
