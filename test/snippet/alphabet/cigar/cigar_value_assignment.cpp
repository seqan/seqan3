#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::cigar letter{10, 'I'_cigar_operation};
    seqan3::debug_stream << "letter: " << letter << '\n'; // 10I

    letter = 'D'_cigar_operation;
    seqan3::debug_stream << "letter: " << letter << '\n'; // 10D

    letter = 20;
    seqan3::debug_stream << "letter: " << letter << '\n'; // 20D
}
