#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/core/type_list/type_list.hpp>
#include <seqan3/io/detail/record.hpp>

int main()
{
    using types         = seqan3::type_list<std::string,       seqan3::dna4_vector, std::vector<seqan3::phred42>>;
    using types_as_ids  = seqan3::fields<   seqan3::field::id, seqan3::field::seq,  seqan3::field::qual>;
    using selected_ids  = seqan3::fields<seqan3::field::qual, seqan3::field::id>;

    using selected_types = seqan3::detail::select_types_with_ids_t<types, types_as_ids, selected_ids>;
    // resolves to type_list<std::vector<phred42>, std::string>>
}
