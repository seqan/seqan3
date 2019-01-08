#include <seqan3/alphabet/aminoacid/concept.hpp>

namespace your_namespace // optional
{

// your own aminoacid definition
struct your_aa
{
    //...
};

} // namespace your_namespace end

// You can copy and paste this code and adjust `your_namespace::your_aa`
// if you want your amino acid to be recognized as such by seqan3.
namespace seqan3
{

template <>
struct is_aminoacid<your_namespace::your_aa> : std::true_type {};

} // namespace seqan3 end
