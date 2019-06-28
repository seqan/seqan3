#include <seqan3/alphabet/mask/mask.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    seqan3::mask my_mask = seqan3::mask::MASKED;
    seqan3::mask another_mask{};

    my_mask.assign_rank(false);  // will assign my_mask the value mask::UNMASKED
    another_mask.assign_rank(0); // will also assign another_mask the value mask::UNMASKED

    if (my_mask.to_rank() == another_mask.to_rank())
        seqan3::debug_stream << "Both are UNMASKED!\n";
}
