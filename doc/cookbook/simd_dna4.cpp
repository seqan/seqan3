// Checks that defining simd_dna4 works without putting it into the seqan3 namespace.

#define SEQAN3_USE_NAMESPACE 0
#include <seqan3/test/performance/simd_dna4.hpp>

int main()
{
    [[maybe_unused]] simd_dna4 letter{};
}
