#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/core/debug_stream.hpp>

using seqan3::operator""_cigar_operation;

int main()
{
    // Initialise a seqan3::cigar::operation:
    seqan3::cigar::operation match{'M'_cigar_operation};

    // you can print cigar_op values:
    seqan3::debug_stream << match << '\n'; // M
}
