#include <gtest/gtest.h>

#include <seqan3/utility/simd/all.hpp>

#include <seqan3/test/simd_utility.hpp>

TEST(simd, simd_eq)
{
    using int16x8_t = seqan3::simd::simd_type_t<int16_t, 8u>;
    int16x8_t a{4, 3, 2, 1, 0, -1, -2, -3};

    SIMD_EQ(a, int16x8_t{4, 3, 2, 1, 0, -1, -2, -3});
}
