//![start]
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/std/ranges>

using namespace seqan3;

//![start]
auto my_reverse_complement = std::view::reverse | std::view::transform([] (auto const d)
{
    return complement(d);
});

//![end]
int main()
{
    std::vector<dna5> vec{"ACCAGATTA"_dna5};
    debug_stream << vec << '\n';                    // will print "ACCAGATTA"

    auto v = vec | my_reverse_complement;

    debug_stream << v << '\n';                      // prints "TAATCTGGT"
}
//![end]
