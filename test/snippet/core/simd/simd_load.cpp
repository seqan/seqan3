#include <seqan3/core/debug_stream.hpp>
#include <seqan3/core/simd/all.hpp>

using uint16x8_t = seqan3::simd::simd_type_t<uint16_t, 8>;

int main()
{
    std::vector<uint16_t> memory{0, 1, 2, 3, 4, 5, 6, 7};
    uint16x8_t a = seqan3::simd::load<uint16x8_t>(memory.data());
    seqan3::debug_stream << a << '\n'; // [0,1,2,3,4,5,6,7]
    return 0;
}
