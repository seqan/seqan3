#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/view/to_char.hpp>
#include <seqan3/std/ranges>

int main()
{
    using seqan3::operator""_dna4;
    seqan3::dna4_vector s{"ACTTTGATAA"_dna4};
    using iterator = seqan3::dna4_vector::iterator;
    auto v1 = std::ranges::subrange<iterator, iterator>{std::ranges::begin(s) + 2, std::ranges::end(s)}
            | seqan3::view::to_char; // == "TTTGATAA"

    seqan3::debug_stream << v1 << '\n';
}
