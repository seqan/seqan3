#include <seqan3/simd/detail/ume_simd.hpp>
using namespace seqan3::detail;

#if __has_include(<umesimd/UMESimd.h>)

using int8x16_t = ume_simd<int8_t, 16>::type;

static_assert(is_ume_simd<int8x16_t>::value);
static_assert(!is_ume_simd<int8_t>::value);

#endif // __has_include(<umesimd/UMESimd.h>)
