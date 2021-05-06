#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/container/dynamic_bitset.hpp>

int main()
{
    seqan3::dynamic_bitset t1{"10001100"};
    seqan3::dynamic_bitset const t2{0b1011'1000};

    t1 |= t2;
    seqan3::debug_stream << t1 << '\n'; // 1011'1100
}
