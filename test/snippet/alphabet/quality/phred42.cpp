//! [general]
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

int main()
{
// doesn't work:
// phred = 2;
// works:
phred42 phred;
phred.assign_rank(2); // wrapper for assign_phred(2)
debug_stream << (int) phred.to_phred() << "\n"; // 2
debug_stream << phred.to_char() << "\n";        // '#'
debug_stream << (int) phred.to_rank() << "\n";  // 2

// doesn't work:
// phred42{4};
// another_phred.assign_rank(49);
// works:
phred42 another_phred;
another_phred.assign_phred(49); // converted down to 41
debug_stream << (int) another_phred.to_phred() << "\n"; // 41
// we need to cast to(int)for human readable console output
}
//! [general]
