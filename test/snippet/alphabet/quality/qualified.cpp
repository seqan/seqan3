#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using seqan3::operator""_dna4;
    using seqan3::get;

    seqan3::qualified<seqan3::dna4, seqan3::phred42> letter{'A'_dna4, seqan3::phred42{7}};
    seqan3::debug_stream << int(seqan3::to_rank(letter)) << ' '
                         << int(seqan3::to_rank(get<0>(letter))) << ' '
                         << int(seqan3::to_rank(get<1>(letter))) << '\n';
    // 28 0 7

    seqan3::debug_stream << seqan3::to_char(letter) << ' '
                         << seqan3::to_char(get<0>(letter)) << ' '
                         << seqan3::to_char(get<1>(letter)) << '\n';
    // A A (

    seqan3::debug_stream << seqan3::to_phred(letter) << ' '
                         << seqan3::to_phred(get<1>(letter)) << '\n';
    // 7 7

    // modify via structured bindings and references:
    auto & [ seq_l, qual_l ] = letter;
    seq_l = 'G'_dna4;
    seqan3::debug_stream << seqan3::to_char(letter) << '\n';
    // G
}
