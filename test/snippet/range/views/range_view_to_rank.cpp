#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/aliases.hpp>
#include <seqan3/alphabet/views/to_rank.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using seqan3::operator""_dna4;

    seqan3::dna4_vector vec = "ACTTTGATA"_dna4;
    auto v = vec | seqan3::views::to_rank;
    seqan3::debug_stream << v << '\n'; // [0,1,3,3,3,2,0,3,0]

    std::vector<seqan3::phred42> qvec{{0}, {7}, {5}, {3}, {7}, {4}, {30}, {16}, {23}};
    auto v3 = qvec | seqan3::views::to_rank;
    seqan3::debug_stream << v3 << '\n'; // [0,7,5,3,7,4,30,16,23]

    std::vector<seqan3::dna4q> qcvec{{'C'_dna4, seqan3::phred42{0}}, {'A'_dna4, seqan3::phred42{7}},
                                     {'G'_dna4, seqan3::phred42{5}}, {'T'_dna4, seqan3::phred42{3}},
                                     {'G'_dna4, seqan3::phred42{7}}, {'A'_dna4, seqan3::phred42{4}},
                                     {'C'_dna4, seqan3::phred42{30}}, {'T'_dna4, seqan3::phred42{16}},
                                     {'A'_dna4, seqan3::phred42{23}}};
    auto v4 = qcvec | seqan3::views::to_rank;
    seqan3::debug_stream << v4 << '\n'; // [42,7,89,129,91,4,72,142,23]
}
