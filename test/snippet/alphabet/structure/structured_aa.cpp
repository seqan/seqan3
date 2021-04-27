#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/alphabet/structure/dssp9.hpp>
#include <seqan3/alphabet/structure/structured_aa.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;
    using seqan3::get;

    seqan3::structured_aa<seqan3::aa27, seqan3::dssp9> letter{'W'_aa27, 'B'_dssp9};

    seqan3::debug_stream << seqan3::to_rank(letter) << ' '
                         << seqan3::to_rank(get<0>(letter)) << ' '
                         << seqan3::to_rank(get<1>(letter)) << '\n';
    // prints "199 22 1"

    seqan3::debug_stream << seqan3::to_char(letter) << ' '
                         << seqan3::to_char(get<0>(letter)) << ' '
                         << seqan3::to_char(get<1>(letter)) << '\n';
    // prints "W W B"

    // Modify via structured bindings and references:
    auto & [ sequence_letter, structure_letter ] = letter;
    sequence_letter = 'V'_aa27;
    seqan3::debug_stream << seqan3::to_char(letter) << '\n';
    // prints "V"
}
