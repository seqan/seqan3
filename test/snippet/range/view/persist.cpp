#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/view/persist.hpp>
#include <seqan3/range/view/to_char.hpp>

int main()
{
    using seqan3::operator""_dna4;

    // explicitly create an l-value of our dna vector:
    auto vec = "ACGT"_dna4;
    auto v = vec | seqan3::view::to_char;

    // using seqan3::view::persist you can bind the temporary directly:
    auto v2 = "ACGT"_dna4 | seqan3::view::persist | seqan3::view::to_char;

    // note that seqan3::view::persist must follow immediately after the temporary,
    // thus the function notation might be more intuitive:
    auto v3 = seqan3::view::persist("ACGT"_dna4) | seqan3::view::to_char;
}
