#include <range/v3/view/zip.hpp>

#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alphabet/aminoacid/all.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

int main()
{

{
//! [two letters]
aminoacid_scoring_scheme scheme{aminoacid_similarity_matrix::BLOSUM62};
debug_stream << "BLOSUM62 score for T and S: " << (int) scheme.score('T'_aa27, 'S'_aa27) << "\n"; // == 1

scheme.set_similarity_matrix(aminoacid_similarity_matrix::BLOSUM80);
debug_stream << "BLOSUM80 score for 'T'_aa27 and 'S'_aa20: " << (int) scheme.score('T'_aa27, 'S'_aa20) << "\n"; // == 2
// you can score aa20 against aa27

scheme.set_hamming_distance();
debug_stream << "Hamming distance between T and S: " << (int) scheme.score('T'_aa27, 'S'_aa20) << "\n"; // == -1
debug_stream << "Hamming distance between T and T: " << (int) scheme.score('T'_aa27, 'T'_aa20) << "\n"; // == 0
//! [two letters]
}

{
//! [edit matrix]
aminoacid_scoring_scheme scheme{aminoacid_similarity_matrix::BLOSUM80};
debug_stream << "BLOSUM80 score between T and S: " << (int) scheme.score('T'_aa27, 'S'_aa27) << "\n"; // == 2
auto & cell = scheme.score('T'_aa27, 'S'_aa27);
cell = 3;
debug_stream << "New score after editing entry: " << (int) scheme.score('T'_aa27, 'S'_aa27) << "\n"; // == 3
//! [edit matrix]
}

{
//! [score sequences]
std::vector<aa27> one = "ALIGATOR"_aa27;
std::vector<aa27> two = "ANIMATOR"_aa27;

aminoacid_scoring_scheme scheme{aminoacid_similarity_matrix::BLOSUM62};

int score = 0;
for (auto pair : std::view::zip(one, two))
    score += scheme.score(std::get<0>(pair), std::get<1>(pair));
debug_stream << "Score: " << score << "\n"; // 4 + -3 + 4 + -3 + 4 + 5 + -1 + 5 = 15
//! [score sequences]
}

}
