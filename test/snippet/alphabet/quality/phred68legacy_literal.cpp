#include <seqan3/alphabet/quality/phred68legacy.hpp>
#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

int main()
{
    // directly assign to a std::vector<phred68legacy> using a string literal
    std::vector<phred68legacy> qual_vec = "###!"_phred68legacy;

    // This is the same as a sequence of char literals:
    std::vector<phred68legacy> qual_vec2 = {'#'_phred68legacy, '#'_phred68legacy, '#'_phred68legacy, '!'_phred68legacy};

    debug_stream << ranges::equal(qual_vec, qual_vec2) << std::endl; // prints 1 (true)
}
