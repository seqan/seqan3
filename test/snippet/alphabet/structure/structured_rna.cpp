#include <seqan3/alphabet/structure/structured_rna.hpp>
#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/alphabet/structure/dot_bracket3.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using seqan3::operator""_rna4;
    using seqan3::operator""_db3;
    using seqan3::get;

    seqan3::structured_rna<seqan3::rna4, seqan3::dot_bracket3> letter{'G'_rna4, '('_db3};
    seqan3::debug_stream << seqan3::to_rank(letter) << ' '
                         << seqan3::to_rank(get<0>(letter)) << ' '
                         << seqan3::to_rank(get<1>(letter)) << '\n';
    // 6 2 1

    seqan3::debug_stream << seqan3::to_char(letter) << ' '
                         << seqan3::to_char(get<0>(letter)) << ' '
                         << seqan3::to_char(get<1>(letter)) << '\n';
    // G G (

    // modify via structured bindings and references:
    auto & [ seq_l, structure_l ] = letter;
    seq_l = 'U'_rna4;
    seqan3::debug_stream << seqan3::to_char(letter) << '\n';
    // U
}
