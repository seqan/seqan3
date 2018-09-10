#include <range/v3/view/zip.hpp>

#include <seqan3/alignment/scoring/aminoacid_scoring_scheme.hpp>
#include <seqan3/alphabet/aminoacid/all.hpp>
#include <seqan3/alphabet/aminoacid/aa27.hpp>

using namespace seqan3;

int main()
{

{
//! [two letters]
aminoacid_scoring_scheme scheme{aminoacid_similarity_matrix::BLOSUM62};
std::cout << "BLOSUM62 score for T and S: " << (int) scheme.score(aa27::T, aa27::S) << "\n"; // == 1

scheme.set_similarity_matrix(aminoacid_similarity_matrix::BLOSUM80);
std::cout << "BLOSUM80 score for aa27::T and aa20::S: " << (int) scheme.score(aa27::T, aa20::S) << "\n"; // == 2
// you can score aa20 against aa27

scheme.set_hamming_distance();
std::cout << "Hamming distance between T and S: " << (int) scheme.score(aa27::T, aa20::S) << "\n"; // == -1
std::cout << "Hamming distance between T and T: " << (int) scheme.score(aa27::T, aa20::T) << "\n"; // == 0
//! [two letters]
}

{
//! [edit matrix]
aminoacid_scoring_scheme scheme;
scheme.set_similarity_matrix(aminoacid_similarity_matrix::BLOSUM80);
std::cout << "BLOSUM80 score between T and S: " << (int) scheme.score(aa27::T, aa27::S) << "\n"; // == 2
auto & cell = scheme.score(aa27::T, aa27::S);
cell = 3;
std::cout << "New score after editing entry: " << (int) scheme.score(aa27::T, aa27::S) << "\n"; // == 3
//! [edit matrix]
}

{
//! [score sequences]
using namespace seqan3::literal;
std::vector<aa27> one = "ALIGATOR"_aa27;
std::vector<aa27> two = "ANIMATOR"_aa27;

aminoacid_scoring_scheme scheme{aminoacid_similarity_matrix::BLOSUM62};

int score = 0;
for (auto pair : ranges::view::zip(one, two))
    score += scheme.score(std::get<0>(pair), std::get<1>(pair));
std::cout << "Score: " << score << "\n"; // 4 + -3 + 4 + -3 + 4 + 5 + -1 + 5 = 15
//! [score sequences]
}

}
