#include <seqan3/alphabet/quality/phred63.hpp>
#include <seqan3/core/debug_stream.hpp>

int main()
{
    using seqan3::operator""_phred63;

    // directly assign to a std::vector<phred63> using a string literal
    std::vector<seqan3::phred63> qual_vec = "###!"_phred63;

    // This is the same as a sequence of char literals:
    std::vector<seqan3::phred63> qual_vec2 = {'#'_phred63, '#'_phred63, '#'_phred63, '!'_phred63};

    seqan3::debug_stream << ranges::equal(qual_vec, qual_vec2) << std::endl; // prints 1 (true)
}
