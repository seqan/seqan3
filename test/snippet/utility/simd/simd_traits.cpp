#include <seqan3/utility/simd/all.hpp>

using uint16x8_t = seqan3::simd::simd_type_t<uint16_t, 8>;

static_assert(std::is_same_v<seqan3::simd::simd_traits<uint16x8_t>::scalar_type, uint16_t>);
static_assert(seqan3::simd::simd_traits<uint16x8_t>::length == 8);
static_assert(seqan3::simd::simd_traits<uint16x8_t>::max_length == 16);
