#include <sstream>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using seqan3::operator""_dna5;

    std::ostringstream o;
    seqan3::debug_stream.set_underlying_stream(o);

    seqan3::debug_stream << "ACGT"_dna5;

    o.flush();
    seqan3::debug_stream << o.str(); // prints the string stream's buffer: "ACGT"
}
