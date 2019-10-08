#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/core/debug_stream.hpp>

using seqan3::operator""_cigar_op;

int main()
{
    seqan3::cigar c{13, 'M'_cigar_op};
    seqan3::debug_stream << c << '\n'; // "13M"

    c.assign_string("14X");
    seqan3::debug_stream << c << '\n'; // "14X"
}
