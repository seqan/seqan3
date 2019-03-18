#include <seqan3/io/detail/record.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>

using namespace seqan3;

int main()
{
//! [usage]
using types         = type_list<std::string, dna4_vector, std::vector<phred42>>;
using types_as_ids  = fields<field::ID,      field::SEQ,  field::QUAL>;
using selected_ids  = fields<field::QUAL, field::ID>;

using selected_types = detail::select_types_with_ids_t<types, types_as_ids, selected_ids>;
// resolves to type_list<std::vector<phred42>, std::string>>
//! [usage]
[[maybe_unused]] selected_types t;
}
