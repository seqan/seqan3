#include <seqan3/core/metafunction/range.hpp>

using namespace seqan3;

int main()
{
//! [usage]
// these evaluate to true:
static_assert(seqan3::compatible_concept<std::string,              std::vector<char>>);
static_assert(seqan3::compatible_concept<std::vector<std::string>, std::vector<std::vector<char>>>);
//! [usage]
}
