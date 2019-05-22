#include <iostream>
#include <seqan3/alphabet/gap/gap.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

int main()
{
//! [general]
gap my_gap = gap{};
gap another_gap{};
another_gap.assign_char('A'); // setting this does not change anything

debug_stream << my_gap.to_char(); // outputs '-'
if (my_gap.to_char() == another_gap.to_char())
    debug_stream << "Both gaps are the same!";
//! [general]
}
