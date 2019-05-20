#include <iostream>
#include <seqan3/alphabet/mask/mask.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

int main()
{
//! [general]
mask my_mask = mask::MASKED;
mask another_mask{};

my_mask.assign_rank(false);  // will assign my_mask the value mask::UNMASKED
another_mask.assign_rank(0); // will also assign another_mask the value mask::UNMASKED

if (my_mask.to_rank() == another_mask.to_rank())
    debug_stream << "Both are UNMASKED!\n";
//! [general]
}
