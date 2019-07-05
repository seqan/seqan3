#include <seqan3/core/bit_manipulation.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    uint8_t  t0 = 0b0000'1001;
    uint16_t t1 = 0b0100'0001'0000'1001;
    uint32_t t2 = 0b0000'0000'0000'0000'0000'0000'0000'0001;
    uint64_t t3 = 0b0100'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000'0000;

    seqan3::debug_stream << seqan3::detail::most_significant_bit_set(t0) << '\n'; // 3
    seqan3::debug_stream << seqan3::detail::most_significant_bit_set(t1) << '\n'; // 14
    seqan3::debug_stream << seqan3::detail::most_significant_bit_set(t2) << '\n'; // 0
    seqan3::debug_stream << seqan3::detail::most_significant_bit_set(t3) << '\n'; // 62

    return 0;
}
