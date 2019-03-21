#include <vector>
#include <seqan3/io/stream/debug_stream.hpp>
#include <seqan3/range/view/pairwise_combine.hpp>

using namespace seqan3;

int main()
{
    std::vector vec{'a', 'b', 'c', 'd'};
    for (auto res : vec | view::pairwise_combine)
    {
        debug_stream << res << "\n";
    }
}
