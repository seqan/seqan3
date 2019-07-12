#include <seqan3/alphabet/quality/phred63.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    seqan3::phred63 phred;
    phred.assign_rank(2); // wrapper for assign_phred(2)
    seqan3::debug_stream << static_cast<int>(phred.to_phred()) << "\n"; // 2
    seqan3::debug_stream << phred.to_char() << "\n";        // '#'
    seqan3::debug_stream << static_cast<int>(phred.to_rank()) << "\n";  // 2

    seqan3::phred63 another_phred{49};
    seqan3::debug_stream << another_phred.to_phred() << "\n"; // 49
    // we need to cast to (int) for human readable console output
}
