#include <iostream>

#include <seqan3/alphabet/gap/gap.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    seqan3::gap my_gap = seqan3::gap{};
    seqan3::gap another_gap{};
    another_gap.assign_char('A'); // this does not change anything

    seqan3::debug_stream << my_gap.to_char(); // outputs '-'
    if (my_gap.to_char() == another_gap.to_char())
        seqan3::debug_stream << "Both gaps are the same!";
}
