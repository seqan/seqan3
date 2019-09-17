#include <string>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/decorator/gap_decorator.hpp>
#include <seqan3/range/views/pseudo_random_access.hpp>

int main()
{
    using seqan3::operator""_dna4;

    // A gap decorator is a pseudo random access range using logarithmic time complexity internally.
    auto seq = "ACGTACGACT"_dna4;
    seqan3::gap_decorator dec{seq};

    // The default interface models at most std::bidirectional_range.
    auto it = std::ranges::begin(dec);          // Returned iterator models std::bidirectional_iterator.
    std::ranges::advance(it, 3);                // Advancing the iterator takes linear time.
    seqan3::debug_stream << *it << '\n';        // "T"

    // After adapting it with pseudo_random_access, the returned range models std::random_access_range.
    auto dec_ra = dec | seqan3::views::pseudo_random_access;
    auto it_ra = std::ranges::begin(dec_ra);    // Returned iterator models std::random_access_iterator.
    std::ranges::advance(it_ra, 3);             // Advancing the iterator takes now logarithmic time.
    seqan3::debug_stream << *it_ra << '\n';     // "T"

}
