#include <seqan3/core/simd/detail/builtin_simd_macros.hpp>

template<typename scalar_t, size_t length>
struct builtin_simd;

#define SEQAN3_BUILTIN_SIMD_CLASS(scalar_t, length, max_length, simd_t) \
template <> \
struct builtin_simd<scalar_t, length> \
{ using type = simd_t; };

SEQAN3_BUILTIN_SIMD_PARTIAL_TEMPLATE_SPECIALIZATION(SEQAN3_BUILTIN_SIMD_CLASS)
