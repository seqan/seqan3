// generated from test/snippet/alphabet/nucleotide/@target_alphabet@_implicit_conversion_from_@source_alphabet@_views.cpp.in

//![main]
#include <vector>

#include <seqan3/alphabet/nucleotide/dna15.hpp>
#include <seqan3/alphabet/nucleotide/rna15.hpp>
#include <seqan3/utility/views/convert.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::rna15_vector vector = "ACG"_rna15;

    auto dna15_view = vector | seqan3::views::convert<seqan3::dna15>;

    for (auto && chr : dna15_view) // converts lazily on-the-fly
    {
        static_assert(std::same_as<decltype(chr), seqan3::dna15 &&>);
    }
}
//![main]
