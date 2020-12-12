#include <seqan3/std/bit>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    uint8_t  t0 = 0b0000'1001;
    uint16_t t1 = 0b0100'0001'0000'1001;
    uint32_t t2 = 0b0000'0000'0000'0000'0000'0000'0000'0001;
    uint64_t t3 = 0b0100'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000;

    seqan3::debug_stream << std::bit_width(t0) << '\n'; // 4
    seqan3::debug_stream << std::bit_width(t1) << '\n'; // 15
    seqan3::debug_stream << std::bit_width(t2) << '\n'; // 1
    seqan3::debug_stream << std::bit_width(t3) << '\n'; // 63

    return 0;
}
