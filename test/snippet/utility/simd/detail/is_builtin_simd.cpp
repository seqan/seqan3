#include <seqan3/utility/simd/detail/builtin_simd.hpp>

using int8x16_t = seqan3::detail::builtin_simd<int8_t, 16>::type;

static_assert(seqan3::detail::is_builtin_simd<int8x16_t>::value);
static_assert(!seqan3::detail::is_builtin_simd<int8_t>::value);
