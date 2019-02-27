#include <seqan3/io/record.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>

using namespace seqan3;

int main()
{

auto tmp_dir = std::filesystem::temp_directory_path();

// Create a /tmp/my.fasta for the snippet below
{
    sequence_file_output fout{tmp_dir/"my.fasta"};
    fout.emplace_back("ACGT"_dna4, "TEST1");
    fout.emplace_back("AGGCTGA"_dna4, "Test2");
    fout.emplace_back("ACTGA"_dna4, "Test2");
}
{
//! [usage_1]
// specify custom field combination/order to file:
sequence_file_input fin{tmp_dir/"my.fasta", fields<field::ID, field::SEQ>{}};

auto record = fin.front(); // get current record, in this case the first

// record is tuple-like type, but allows access via field identifiers:
auto & id = get<field::ID>(record);
auto & seq = get<field::SEQ>(record);
//! [usage_1]
(void) id;
(void) seq;
}

{
//! [usage_2]
using types        = type_list<dna4_vector, std::string, std::vector<phred42>>;
using types_as_ids =    fields<field::SEQ,  field::ID,   field::QUAL>;

using record_type  = record<types, types_as_ids>;
// record_type now mimics std::tuple<std::string, dna4_vector, std::vector<phred42>>, the order also depends on selected_ids

record_type my_record;
get<1>(my_record) = "the most important sequence in the database";   // access via index
get<field::SEQ>(my_record) = "ACGT"_dna4;                            // access via seqan3::field
get<std::string>(my_record) = "the least important sequence in the database";   // access via type
//! [usage_2]
}
}
