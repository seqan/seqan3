//! [get_cigar_vector]
#include <seqan3/alphabet/cigar/all.hpp>
#include <seqan3/alphabet/gap/gapped.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>

#include <range/v3/algorithm/equal.hpp>

using namespace seqan3;

int main()
{
    // The aligned sequences:
    std::vector<gapped<dna4>> ref{dna4::T, dna4::A, gap::GAP, dna4::C, dna4::T, dna4::T};  // TA-CTT
    std::vector<gapped<dna4>> read{dna4::A, dna4::A, dna4::G, dna4::C, gap::GAP, dna4::T}; // AAGC-T

    // get the corresponding cigar string
    cigar_vector basic_cigar = get_cigar_vector(std::make_pair(ref, read));
    // This is the same as the following (default values)
    cigar_vector basic_cigar2 = get_cigar_vector(std::make_pair(ref, read), 0, 0, false);

    std::cout << basic_cigar << std::endl; // 2M1I1M1D1M

    if (ranges::equal(basic_cigar, basic_cigar2))
        std::cout << "yeah" << std::endl; // yeah

    // You can add soft clipping to your cigar
    cigar_vector soft_clipped_cigar = get_cigar_vector(std::make_pair(ref, read), 5, 10);

    std::cout << soft_clipped_cigar << std::endl; // 5S2M1I1M1D1M10S

    // You can also get the cigar in the extended cigar alphabet
    // where we differentiate between a match (=) and a mismatch (X)
    cigar_vector extended_cigar = get_cigar_vector(std::make_pair(ref, read), 5, 10, true);

    std::cout << soft_clipped_cigar << std::endl; // 5S1X1=1I1=1D1=10S
}
//! [get_cigar_vector]
