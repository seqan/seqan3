#include <seqan3/alphabet/quality/phred68legacy.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    seqan3::phred68legacy phred;
    phred.assign_phred(-2);
    seqan3::debug_stream << (int) phred.to_phred() << "\n"; // -2
    seqan3::debug_stream << phred.to_char() << "\n";        // '>'
    seqan3::debug_stream << (int) phred.to_rank() << "\n";  // 3
}
