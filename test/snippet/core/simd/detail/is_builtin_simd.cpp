#include <seqan3/core/simd/detail/builtin_simd.hpp>
using namespace seqan3::detail;

using int8x16_t = builtin_simd<int8_t, 16>::type;

static_assert(is_builtin_simd<int8x16_t>::value);
static_assert(!is_builtin_simd<int8_t>::value);
