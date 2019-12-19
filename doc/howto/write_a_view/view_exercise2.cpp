//![start]
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/std/ranges>

using seqan3::operator""_dna5;

//![start]
auto my_reverse_complement = std::views::reverse | std::views::transform([] (auto const d)
{
    return seqan3::complement(d);
});

//![end]
int main()
{
    std::vector<seqan3::dna5> vec{"ACCAGATTA"_dna5};
    seqan3::debug_stream << vec << '\n';            // will print "ACCAGATTA"

    auto v = vec | my_reverse_complement;

    seqan3::debug_stream << v << '\n';              // prints "TAATCTGGT"
}
//![end]
