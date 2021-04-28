#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using seqan3::get;
    using namespace seqan3::literals;

    seqan3::cigar letter{10, 'M'_cigar_operation};

    // Note that this is equivalent to get<uint32_t>(letter)
    uint32_t size{get<0>(letter)};

    // Note that this is equivalent to get<seqan3::cigar::operation>(letter)
    seqan3::cigar::operation cigar_char{get<1>(letter)};

    seqan3::debug_stream << "Size is "       << size       << '\n';
    seqan3::debug_stream << "Cigar char is " << cigar_char << '\n'; // seqan3::debug_stream converts to char on the fly.

}
