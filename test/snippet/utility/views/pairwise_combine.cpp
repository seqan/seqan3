#include <vector>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/utility/views/pairwise_combine.hpp>

int main()
{
    std::vector vec{'a', 'b', 'c', 'a'};
    for (auto res : vec | seqan3::views::pairwise_combine)
    {
        seqan3::debug_stream << res << '\n';
    }

    // Possible Output:
    // (a,b)
    // (a,c)
    // (a,a)
    // (b,c)
    // (b,a)
    // (c,a)
}
