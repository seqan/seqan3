#include <range/v3/view/zip.hpp>

#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/nucleotide/rna15.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

int main()
{

{
//! [two letters]
nucleotide_scoring_scheme scheme; // hamming is default
debug_stream << "Score between DNA5 A and G: " << (int) scheme.score('A'_dna5, 'G'_dna5) << "\n"; // == -1
debug_stream << "Score between DNA5 A and A: " << (int) scheme.score('A'_dna5, 'A'_dna5) << "\n"; // == 0

scheme.set_simple_scheme(match_score{3}, mismatch_score{-2});
debug_stream << "Score between DNA5 A and RNA15 G: " << (int) scheme.score('A'_dna5, 'G'_rna15) << "\n"; // == -2
debug_stream << "Score between DNA5 A and RNA15 A: " << (int) scheme.score('A'_dna5, 'A'_rna15) << "\n"; // == 3
// you can score differenct nucleotides  ^
//! [two letters]
}

{
//! [edit matrix]
nucleotide_scoring_scheme scheme; // hamming distance is default
debug_stream << "Score between DNA A and G before edit: "
          << (int) scheme.score('A'_dna15, 'G'_dna15) << "\n"; // == -1
scheme.score('A'_dna15, 'G'_dna15) = 3;
debug_stream << "Score after editing: " << (int) scheme.score('A'_dna15, 'G'_dna15) << "\n"; // == 3
//! [edit matrix]
}

{
//! [score sequences]
std::vector<dna15> one = "AGAATA"_dna15;
std::vector<dna15> two = "ATACTA"_dna15;
nucleotide_scoring_scheme scheme; // hamming distance is default

int score = 0;
for (auto pair : std::view::zip(one, two))
    score += scheme.score(std::get<0>(pair), std::get<1>(pair));
debug_stream << "Score: " << score << "\n"; // == 0 - 1 + 0 - 1 + 0 + 0 = -2
//! [score sequences]
}

}
