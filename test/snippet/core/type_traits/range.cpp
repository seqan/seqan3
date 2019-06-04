#include <seqan3/core/type_traits/range.hpp>

using namespace seqan3;

int main()
{
//! [usage]
// these evaluate to true:
static_assert(seqan3::Compatible<std::string,              std::vector<char>>);
static_assert(seqan3::Compatible<std::vector<std::string>, std::vector<std::vector<char>>>);
//! [usage]
}
