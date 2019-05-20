#include <seqan3/alphabet/quality/phred63.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

int main()
{
//! [general]
// doesn't work:
// phred = 2;
// works:
phred63 phred;
phred.assign_rank(2); // wrapper for assign_phred(2)
debug_stream << (int) phred.to_phred() << "\n"; // 2
debug_stream << phred.to_char() << "\n";        // '#'
debug_stream << (int) phred.to_rank() << "\n";  // 2

// phred = 75; // <- throws assertion
// doesn't work:
// phred63{4};
// works:
phred63 another_phred{49};
debug_stream << (int) another_phred.to_phred() << "\n"; // 49
// we need to cast to(int)for human readable console output
//! [general]
}
