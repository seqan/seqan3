#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/container/dynamic_bitset.hpp>

int main()
{
    seqan3::dynamic_bitset t1{"1100"};

    t1.resize(8);        // Enlarge to 8.
    seqan3::debug_stream << t1 << '\n'; // 0000'1100
    t1.resize(5);        // Shrink to 5, last three bits (5, 6, 7) are set to false.
    seqan3::debug_stream << t1 << '\n'; // 0110'0
    t1.resize(10, true); // Enlarge to 10 and set new bits to true.
    seqan3::debug_stream << t1 << '\n'; // 1111'1011'00
}
