#include <gtest/gtest.h>

#include <seqan3/core/simd/all.hpp>
#include <seqan3/test/simd_utility.hpp>

TEST(simd, simd_eq)
{
    using int16x4_t = seqan3::simd::simd_type_t<int16_t, 8u>;
    int16x4_t a{4, 3, 2, 1, 0, -1, -2, -3};

    SIMD_EQ(a, int16x4_t{4, 3, 2, 1, 0, -1, -2, -3});
}
