#include <seqan3/core/debug_stream.hpp>
#include <seqan3/core/simd/all.hpp>

using int16x8_t = seqan3::simd::simd_type_t<int16_t, 8>;

int main()
{
    int16x8_t a{0, -1, 2, -3, 4, -5, 6, -7};

    int16x8_t b = seqan3::detail::extract_halve<1>(a);
    int16x8_t c = seqan3::detail::extract_quarter<1>(a);
    int16x8_t d = seqan3::detail::extract_eighth<7>(a);

    seqan3::debug_stream << b << '\n'; // [4,-5,6,-7,0,0,0,0]
    seqan3::debug_stream << c << '\n'; // [2,-3,0,0,0,0,0,0]
    seqan3::debug_stream << d << '\n'; // [-7,0,0,0,0,0,0,0]
    return 0;
}
