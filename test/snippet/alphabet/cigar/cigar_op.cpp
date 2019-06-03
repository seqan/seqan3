#include <seqan3/alphabet/cigar/cigar_op.hpp>
#include <seqan3/core/debug_stream.hpp>

using seqan3::operator""_cigar_op;

int main()
{
// Initialze an seqan3::cigar_op:
seqan3::cigar_op match{'M'_cigar_op};

// you can print cigar_op values:
seqan3::debug_stream << match << std::endl;          // M
}
