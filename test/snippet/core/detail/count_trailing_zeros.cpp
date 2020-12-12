#include <seqan3/std/bit>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    uint8_t  t0 = 0b0000'1001;
    uint16_t t1 = 0b0100'0001'0000'1000;
    uint32_t t2 = 0b0000'0000'0000'0000'1000'0000'0000'0000;
    uint64_t t3 = 0b0000'0000'0000'0000'0010'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000;

    seqan3::debug_stream << std::countr_zero(t0) << '\n'; // 0
    seqan3::debug_stream << std::countr_zero(t1) << '\n'; // 3
    seqan3::debug_stream << std::countr_zero(t2) << '\n'; // 15
    seqan3::debug_stream << std::countr_zero(t3) << '\n'; // 45

    return 0;
}
