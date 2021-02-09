#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/alphabet/cigar/cigar_op.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using seqan3::operator""_cigar_op;

    seqan3::cigar letter{10, 'I'_cigar_op};
    seqan3::debug_stream << "letter: " << letter << '\n'; // 10I

    letter = 'D'_cigar_op;
    seqan3::debug_stream << "letter: " << letter << '\n'; // 10D

    letter = 20;
    seqan3::debug_stream << "letter: " << letter << '\n'; // 20D
}
