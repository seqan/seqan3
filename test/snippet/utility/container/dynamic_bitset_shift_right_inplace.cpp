#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/container/dynamic_bitset.hpp>

int main()
{
    seqan3::dynamic_bitset t1{"10001100"};

    t1 >>= 3;
    seqan3::debug_stream << t1 << '\n'; // 0001'0001
}
