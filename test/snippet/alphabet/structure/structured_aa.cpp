#include <seqan3/alphabet/structure/structured_aa.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using seqan3::operator""_aa27;
    using seqan3::operator""_dssp9;
    using seqan3::get;

    seqan3::structured_aa<seqan3::aa27, seqan3::dssp9> letter{'W'_aa27, 'B'_dssp9};

    seqan3::debug_stream << seqan3::to_rank(letter) << ' '
                         << seqan3::to_rank(get<0>(letter)) << ' '
                         << seqan3::to_rank(get<1>(letter)) << '\n';
    // 49 22 1

    seqan3::debug_stream << seqan3::to_char(letter) << ' '
                         << seqan3::to_char(get<0>(letter)) << ' '
                         << seqan3::to_char(get<1>(letter)) << '\n';
    // W W B

    // modify via structured bindings and references:
    auto & [ seq_l, structure_l ] = letter;
    seq_l = 'V'_aa27;
    seqan3::debug_stream << seqan3::to_char(letter) << '\n';
    // V
}
