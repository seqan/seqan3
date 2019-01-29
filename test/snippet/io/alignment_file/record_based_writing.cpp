//! [all]
#include <seqan3/io/alignment_file/all.hpp>

using namespace seqan3;

int main()
{
    alignment_file_output fout{filesystem::temp_directory_path()/"my.sam"};

    std::string ref_id;
    std::string read_id;

    std::vector<dna5> read;

    // ... e.g. compute and alignment

    using alignment_type = std::pair<std::vector<gapped<dna5>>, std::vector<gapped<dna5>>>;

    alignment_type dummy_alignment{}; // an empty dummy alignment

    using types        = type_list<std::vector<dna5>, std::string, alignment_type>;
    using types_as_ids = fields<field::SEQ, field::ID, field::ALIGNMENT>;
    // the record type specifies the fields we want to write
    using record_type  = record<types, types_as_ids>;

    // initialize record
    record_type rec{read, ref_id, dummy_alignment};

    // Write the record
    fout.push_back(rec);

    // same as
    fout.push_back(record_type{read, ref_id, dummy_alignment});

    // as all our fields are empty so this would print an
}
//! [all]
