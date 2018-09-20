#include <seqan3/alphabet/quality/phred42.hpp>

using namespace seqan3;

int main()
{
//! [general]
// doesn't work:
// phred = 2;
// works:
phred42 phred;
phred.assign_rank(2); // wrapper for assign_phred(2)
std::cout << (int) phred.to_phred() << "\n"; // 2
std::cout << phred.to_char() << "\n";        // '#'
std::cout << (int) phred.to_rank() << "\n";  // 2

// doesn't work:
// phred42{4};
// works:
phred42 another_phred;
another_phred.assign_rank(49); // converted down to 41
std::cout << (int) another_phred.to_phred() << "\n"; // 41
// we need to cast to(int)for human readable console output
//! [general]
}
