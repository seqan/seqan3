#include <seqan3/alphabet/quality/phred68legacy.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;
int main()
{
//! [general]
phred68legacy phred;
phred.assign_phred(-2);
debug_stream << (int) phred.to_phred() << "\n"; // -2
debug_stream << phred.to_char() << "\n";        // '>'
debug_stream << (int) phred.to_rank() << "\n";  // 3
// this doesn't work:
// phred68legacy{4};
// phred = 75; // <- throws assertion
//! [general]
}
