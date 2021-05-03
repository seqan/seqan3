#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/container/dynamic_bitset.hpp>

int main()
{
    seqan3::dynamic_bitset t1{0b1011'1000'1111};

    // begin() refers to the rightmost position.
    for (auto it = t1.begin(); it != t1.end(); ++it)
        seqan3::debug_stream << *it; // 1111'0001'1101

    seqan3::debug_stream << '\n';
}
