#include <range/v3/view/zip.hpp>

#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/nucleotide/rna15.hpp>
#include <seqan3/io/stream/debug_stream.hpp>

using namespace seqan3;

int main()
{

{
//! [two letters]
nucleotide_scoring_scheme scheme; // hamming is default
debug_stream << "Score between DNA5 A and G: " << (int) scheme.score(dna5::A, dna5::G) << "\n"; // == -1
debug_stream << "Score between DNA5 A and A: " << (int) scheme.score(dna5::A, dna5::A) << "\n"; // == 0

scheme.set_simple_scheme(match_score{3}, mismatch_score{-2});
debug_stream << "Score between DNA5 A and RNA15 G: " << (int) scheme.score(dna5::A, rna15::G) << "\n"; // == -2
debug_stream << "Score between DNA5 A and RNA15 A: " << (int) scheme.score(dna5::A, rna15::A) << "\n"; // == 3
// you can score differenct nucleotides  ^
//! [two letters]
}

{
//! [edit matrix]
nucleotide_scoring_scheme scheme; // hamming distance is default
debug_stream << "Score between DNA A and G before edit: "
          << (int) scheme.score(dna15::A, dna15::G) << "\n"; // == -1
scheme.score(dna15::A, dna15::G) = 3;
debug_stream << "Score after editing: " << (int) scheme.score(dna15::A, dna15::G) << "\n"; // == 3
//! [edit matrix]
}

{
//! [score sequences]
using namespace seqan3::literal;
std::vector<dna15> one = "AGAATA"_dna15;
std::vector<dna15> two = "ATACTA"_dna15;
nucleotide_scoring_scheme scheme; // hamming distance is default

int score = 0;
for (auto pair : ranges::view::zip(one, two))
    score += scheme.score(std::get<0>(pair), std::get<1>(pair));
debug_stream << "Score: " << score << "\n"; // == 0 - 1 + 0 - 1 + 0 + 0 = -2
//! [score sequences]
}

}
