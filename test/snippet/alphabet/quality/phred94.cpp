#include <seqan3/alphabet/quality/phred94.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    seqan3::phred94 phred;
    phred.assign_phred(2);
    seqan3::debug_stream << phred.to_phred() << "\n"; // 2
    seqan3::debug_stream << phred.to_char() << "\n";  // '#'
    seqan3::debug_stream << phred.to_rank() << "\n";  // 2

    seqan3::phred94 another_phred{75};
    seqan3::debug_stream << another_phred.to_phred() << "\n"; // 75

    seqan3::phred94 a_third_phred;
    another_phred.assign_phred(105); // converted down to 93
    seqan3::debug_stream << another_phred.to_phred() << "\n"; // 93
}
