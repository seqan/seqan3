#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/simd/all.hpp>

using uint16x8_t = seqan3::simd::simd_type_t<uint16_t, 8>;

int main()
{
    uint16x8_t a = seqan3::simd::iota<uint16x8_t>(0);
    std::vector<uint16_t> memory(seqan3::simd::simd_traits<uint16x8_t>::length);
    seqan3::simd::store(memory.data(), a);
    seqan3::debug_stream << memory << '\n'; // [0,1,2,3,4,5,6,7]
    return 0;
}
