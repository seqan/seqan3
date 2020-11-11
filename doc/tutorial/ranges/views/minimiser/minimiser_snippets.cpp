#include <iostream>
#include <vector>
#include <seqan3/std/ranges> // include all of the standard library's views
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/range/views/kmer_hash.hpp>
#include <seqan3/range/views/minimiser.hpp>

using seqan3::operator""_dna4;

int main()
{
std::vector<seqan3::dna4> text{"CCACGTCGACGGTT"_dna4};
{
//![minimiser_win8]
//![def]
auto minimisers = text | seqan3::views::kmer_hash(seqan3::ungapped{4}) | seqan3::views::minimiser(5);
//![minimiser_win8]
}

    {
//![minimiser_win4]
//![def]
        auto minimisers = text | seqan3::views::kmer_hash(seqan3::ungapped{4}) | seqan3::views::minimiser(1);
//![minimiser_win4]
    }

{
//![minimiser_seed]
//![def]
uint64_t const seed = 0x8F3F73B5CF1C9ADE;
auto minimisers = text | seqan3::views::kmer_hash(seqan3::ungapped{4})
                       | std::views::transform([seed] (uint64_t i)
                                               {return i ^ seed;})
                       | seqan3::views::minimiser(5);
//![minimiser_seed]
}
}
