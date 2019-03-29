#include <seqan3/alphabet/cigar/cigar_op.hpp>

using namespace seqan3;

int main()
{
//! [general]
// Initialze an seqan3::cigar_op:
cigar_op match{'M'_cigar_op};

// Usually you want to use the cigar_op directly in a cigar element
// cigar my_cigar_elem1{match, 20u};

// you print cigar_op and cigar values:
std::cout << match.to_char() << std::endl;          // M
// std::cout << my_cigar_elem1 << std::endl; // 20M

// You can also make your own cigar_vector and print it
// cigar my_cigar_elem2{cigar_op::D, 5u};
// cigar my_cigar_elem3{cigar_op::M, 10u};
//
// cigar_vector my_cigar_vector{my_cigar_elem1, my_cigar_elem2, my_cigar_elem3};
//
// std::cout << my_cigar_vector << std::endl; // 20M5D10M
// }
//! [general]
}
