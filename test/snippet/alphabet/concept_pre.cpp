#include <seqan3/alphabet/concept_pre.hpp>
#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/alphabet/structure/dot_bracket3.hpp>
#include <seqan3/alphabet/structure/structured_rna.hpp>

using namespace seqan3;

int main()
{

{
//! [include order]
#include <seqan3/alphabet/concept_pre.hpp>

// your custom alphabet

#include <seqan3/alphabet/concept.hpp>
//! [include order]
}

{
using alphabet_type = rna4;
//! [value retrieval]
auto i = seqan3::alphabet_size<alphabet_type>::value;
// or
auto j = seqan3::alphabet_size_v<alphabet_type>;
//! [value retrieval]
(void) i;
(void) j;
}

{
using alphabet_type = structured_rna<rna4, dot_bracket3>;
//! [pseudoknot value retrieval]
uint8_t pk_support = seqan3::pseudoknot_support<alphabet_type>::value;
// or
uint8_t pk_support_2 = seqan3::pseudoknot_support_v<alphabet_type>;
//! [pseudoknot value retrieval]
(void) pk_support;
(void) pk_support_2;
}

}
