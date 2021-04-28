#include <seqan3/alphabet/quality/phred68solexa.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::phred68solexa letter{'@'_phred68solexa};

    letter.assign_char(';');
    seqan3::debug_stream << letter.to_phred() << '\n'; // prints "-5"
    seqan3::debug_stream << letter.to_char() << '\n';  // prints ";"

    letter.assign_phred(72); // Values exceeding the maximum are implicitly limited to the maximum phred value.
    seqan3::debug_stream << letter.to_phred() << '\n'; // prints "62"
}
