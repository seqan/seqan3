#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/alphabet/cigar/cigar_op.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using seqan3::get;
    using seqan3::operator""_cigar_op;

    seqan3::cigar letter{10, 'M'_cigar_op};

    uint32_t size{get<uint32_t>(letter)};                       // Note this is equivalent to get<0>(letter)
    seqan3::cigar_op cigar_char{get<seqan3::cigar_op>(letter)}; // Note this is equivalent to get<1>(letter)

    seqan3::debug_stream << "Size is "       << size       << '\n';
    seqan3::debug_stream << "Cigar char is " << cigar_char << '\n'; // seqan3::debug_stream converts to char on the fly.

}
