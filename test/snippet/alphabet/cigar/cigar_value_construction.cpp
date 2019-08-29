#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/alphabet/cigar/cigar_op.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using seqan3::operator""_cigar_op;

    seqan3::cigar letter1{0};
    // creates 0M, as the cigar_op field is not provided.
    seqan3::cigar letter2{'M'_cigar_op};
    // creates 0M, as the integer field is not provided.

    if (letter1 == letter2)
        seqan3::debug_stream << "yeah\n"; // yeah

    seqan3::cigar letter3{10, 'I'_cigar_op};
    // creates 10I, as both fields are explicitly given.
}
