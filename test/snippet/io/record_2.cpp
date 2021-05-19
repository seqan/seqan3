#include <string>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/io/record.hpp>
#include <seqan3/utility/type_list/type_list.hpp>

int main()
{
    using namespace seqan3::literals;

    // The order of the types below represent a mapping between the type and the key.
    using types        = seqan3::type_list<seqan3::dna4_vector, std::string, std::vector<seqan3::phred42>>;
    using types_as_ids = seqan3::fields<seqan3::field::seq, seqan3::field::id, seqan3::field::qual>;
    using record_type  = seqan3::record<types, types_as_ids>;
    // record_type now mimics std::tuple<std::string, dna4_vector, std::vector<phred42>>,
    // the order also depends on selected_ids

    record_type my_record{};
    std::get<1>(my_record) = "the most important sequence in the database";            // access via index
    std::get<std::string>(my_record) = "the least important sequence in the database"; // access via type

#ifdef SEQAN3_DEPRECATED_310
    using seqan3::get;
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
    get<seqan3::field::seq>(my_record) = "ACGT"_dna4;                             // access via seqan3::field
#pragma GCC diagnostic pop
#endif //SEQAN3_DEPRECATED_310
}
