#include <seqan3/alphabet/quality/phred68legacy.hpp>

using namespace seqan3;
int main()
{
//! [general]
phred68legacy phred;
phred.assign_phred(-2);
std::cout << (int) phred.to_phred() << "\n";  // -2
std::cout << phred.to_char() << "\n";  // '>'
std::cout << (int) phred.to_rank() << "\n";  // 3
// this doesn't work:
// phred68legacy{4};
// phred = 75; // <- throws assertion
//! [general]
}
