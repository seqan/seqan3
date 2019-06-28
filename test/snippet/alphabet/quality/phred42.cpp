#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    seqan3::phred42 phred;
    phred.assign_rank(2); // wrapper for assign_phred(2)
    seqan3::debug_stream << (int) phred.to_phred() << "\n"; // 2
    seqan3::debug_stream << phred.to_char() << "\n";        // '#'
    seqan3::debug_stream << (int) phred.to_rank() << "\n";  // 2

    seqan3::phred42 another_phred;
    another_phred.assign_phred(49); // converted down to 41
    seqan3::debug_stream << another_phred.to_phred() << "\n"; // 41
    // we need to cast to(int)for human readable console output
}
