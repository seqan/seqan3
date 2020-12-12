#include <seqan3/std/bit>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    uint8_t  t0 = 0b1100'1001;
    uint16_t t1 = 0b0100'0001'1110'1001;
    uint32_t t2 = 0b0000'0000'0000'0000'0000'0000'0000'0000;
    uint64_t t3 = 0b1000'0000'1111'0000'0000'0000'0000'0110'0000'0000'0000'0000'1110'0000'0000'0001;

    seqan3::debug_stream << std::popcount(t0) << '\n'; // 4
    seqan3::debug_stream << std::popcount(t1) << '\n'; // 7
    seqan3::debug_stream << std::popcount(t2) << '\n'; // 0
    seqan3::debug_stream << std::popcount(t3) << '\n'; // 11

    return 0;
}
