#include <seqan3/alphabet/quality/phred94.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::phred94 letter{'@'_phred94};

    letter.assign_char('!');
    seqan3::debug_stream << letter.to_phred() << '\n'; // prints "0"
    seqan3::debug_stream << letter.to_char() << '\n';  // prints "!"

    letter.assign_phred(99); // Values exceeding the maximum are implicitly limited to the maximum phred value.
    seqan3::debug_stream << letter.to_phred() << '\n'; // prints "93"
}
