#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/alphabet/cigar/cigar_op.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using seqan3::operator""_cigar_op;

    seqan3::cigar letter1{10, 'I'_cigar_op};

    letter1 = 'D'_cigar_op;  // yields 10D
    letter1 = 20; // yields 20D

    if (letter1 == seqan3::cigar{20, 'D'_cigar_op})
        seqan3::debug_stream << "yeah\n";
}
