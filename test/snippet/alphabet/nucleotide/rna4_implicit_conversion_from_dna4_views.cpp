// generated from test/snippet/alphabet/nucleotide/@target_alphabet@_implicit_conversion_from_@source_alphabet@_views.cpp.in

//![main]
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/utility/views/convert.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::rna4_vector vector = "ACG"_rna4;

    auto dna4_view = vector | seqan3::views::convert<seqan3::dna4>;

    for (auto && chr : dna4_view) // converts lazily on-the-fly
    {
        static_assert(std::same_as<decltype(chr), seqan3::dna4 &&>);
    }
}
//![main]
