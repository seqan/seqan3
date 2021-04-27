#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/views/to_char.hpp>
#include <seqan3/core/detail/persist_view.hpp>

int main()
{
    using namespace seqan3::literals;

    // explicitly create an l-value of our dna vector:
    auto vec = "ACGT"_dna4;
    auto v = vec | seqan3::views::to_char;

    // using seqan3::detail::persist you can bind the temporary directly:
    auto v2 = "ACGT"_dna4 | seqan3::detail::persist | seqan3::views::to_char;

    // note that seqan3::detail::persist must follow immediately after the temporary,
    // thus the function notation might be more intuitive:
    auto v3 = seqan3::detail::persist("ACGT"_dna4) | seqan3::views::to_char;
}
