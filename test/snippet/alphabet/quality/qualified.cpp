#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;
    using seqan3::get;

    seqan3::qualified<seqan3::dna4, seqan3::phred42> letter{'A'_dna4, '('_phred42};

    seqan3::debug_stream << seqan3::to_rank(letter) << ' '
                         << seqan3::to_rank(get<0>(letter)) << ' '
                         << seqan3::to_rank(get<1>(letter)) << '\n';
    // prints "7 0 7"

    seqan3::debug_stream << seqan3::to_char(letter) << ' '
                         << seqan3::to_char(get<0>(letter)) << ' '
                         << seqan3::to_char(get<1>(letter)) << '\n';
    // prints "A A ("

    seqan3::debug_stream << seqan3::to_phred(letter) << ' '
                         << seqan3::to_phred(get<1>(letter)) << '\n';
    // prints "7 7"

    // Modify via structured bindings and references:
    auto & [ sequence_letter, quality_letter ] = letter;
    sequence_letter = 'G'_dna4;
    seqan3::debug_stream << seqan3::to_char(letter) << '\n';
    // prints "G"
}
