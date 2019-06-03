#include <iostream>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/range/view/persist.hpp>
#include <seqan3/std/ranges>

using namespace seqan3;

// For snippet col_based_writing
//![data_storage]
struct data_storage_t
{
    concatenated_sequences<dna4_vector> sequences{"ACGT"_dna4, "AAA"_dna4};
    concatenated_sequences<std::string> ids{std::string{"ID1"}, std::string{"ID2"}};
};

data_storage_t data_storage; // a global or globally used variable in your program
//![data_storage]

auto tmp_dir = std::filesystem::temp_directory_path();

int main()
{
{
// First make a /tmp/input.fastq file.
sequence_file_output fout{tmp_dir/"input.fastq"};
fout.emplace_back("ACGT"_dna4, "TEST1", "##!#"_phred42);
fout.emplace_back("AGGCTGA"_dna4, "Test2", "##!#!!!"_phred42);
fout.emplace_back("GGAGTATAATATATATATATATAT"_dna4, "Test3", "##!###!###!###!###!###!#"_phred42);
}

{
//! [template_deduction]
sequence_file_output fout{tmp_dir/"my.fasta"}; // FastA format detected, std::ofstream opened for file
//! [template_deduction]
}

{
//! [cout_write]
sequence_file_output fout{std::cout, format_fasta{}};
//              ^ no need to specify the template arguments

fout.emplace_back("ACGTN"_dna5, "example_id"); // default order for fasta: SEQ, ID
//! [cout_write]
}

{
//! [record_wise_iteration]
sequence_file_output fout{tmp_dir/"my.fasta"};

for (int i = 0; i < 5; i++) // ...
{
    std::string id{"test_id"};
    dna5_vector seq{"ACGT"_dna5};

    // ...

    fout.emplace_back(seq, id);          // as individual variables
    // or:
    fout.push_back(std::tie(seq, id));   // as a tuple
}
//! [record_wise_iteration]
}

{
//! [fields_trait_1]
sequence_file_output fout{tmp_dir/"output.fastq", fields<field::ID, field::SEQ_QUAL>{}};

for (int i = 0; i < 5; i++)
{
    std::string id{"test_id"};
    // vector of combined data structure:
    std::vector<qualified<dna5, phred42>> seq_qual{qualified<dna5, phred42>{'N'_dna5, phred42{7}}};

    // ...

    fout.emplace_back(id, seq_qual);       // note also that the order the argumets is now different, because
    // or:                                    you specified that ID should be first in the fields template argument
    fout.push_back(std::tie(id, seq_qual));
}
//! [fields_trait_1]
}

{
//! [fields_trait_2]
sequence_file_input  fin{tmp_dir/"input.fastq", fields<field::SEQ, field::ID, field::QUAL>{}};
sequence_file_output fout{tmp_dir/"output.fastq"}; // doesn't have to match the configuration

for (auto & r : fin)
{
    if (true) // r fulfills some criterium
    {
        fout.push_back(r);
    }
}
//! [fields_trait_2]
}

{
//! [batch_write]
sequence_file_output fout{tmp_dir/"my.fasta"};

std::vector<std::tuple<dna5_vector, std::string>> range
{
    { "ACGT"_dna5, "First" },
    { "NATA"_dna5, "2nd" },
    { "GATA"_dna5, "Third" }
}; // a range of "records"

fout = range; // will iterate over the records and write them
//! [batch_write]
}

{
//! [direct_writing]
// file format conversion in one line:
sequence_file_output{tmp_dir/"output.fasta"} = sequence_file_input{tmp_dir/"input.fastq"};

// with sequence_file_output as a variable:
sequence_file_output fout{tmp_dir/"output.fasta"};
sequence_file_input fin{tmp_dir/"input.fastq"};
fout = fin;

// or in pipe notation:
sequence_file_input{tmp_dir/"input.fastq"} | sequence_file_output{tmp_dir/"output.fasta"};
//! [direct_writing]
}

{
//! [view_pipeline]

// minimum_average_quality_filter and minimum_sequence_length_filter need to be implemented first
auto minimum_sequence_length_filter = std::view::filter([] (auto rec)
{
    return size(get<field::SEQ>(rec)) >= 50;
});

auto minimum_average_quality_filter = std::view::filter([] (auto const & rec)
{
    double qual_sum{0}; // summation of the qualities
    for (auto chr : get<field::QUAL>(rec))
        qual_sum += chr.to_phred();

    return qual_sum / (size(get<field::QUAL>(rec))) >= 20; // check if average quality is greater than 20.
});

sequence_file_input{tmp_dir/"input.fastq"} | view::persist
                                           | minimum_average_quality_filter
                                           | minimum_sequence_length_filter
                                           | std::view::take(5)
                                           | sequence_file_output{tmp_dir/"output.fasta"};
//! [view_pipeline]
}

{
// see data_storage snippet above main function for completeness
//! [col_based_writing]
// ... in your file writing function:

sequence_file_output fout{tmp_dir/"my.fasta"};

fout = std::tie(data_storage.sequences, data_storage.ids);
//! [col_based_writing]
}

{
//! [range_interface]
sequence_file_output fout{tmp_dir/"my.fasta"};

auto it = fout.begin();

for(int i = 0; i < 5; i++) // some criteria
{
    std::string id{"test_id"};
    dna5_vector seq{"ACGT"_dna5};

    // ...

    // assign to iterator
    *it = std::tie(seq, id);
    // is the same as:
    fout.push_back(std::tie(seq, id));
}
//! [range_interface]
}

{
//! [push_back_record]
sequence_file_output fout{tmp_dir/"my.fasta"};
for(int i = 0; i < 5; i++) // some criteria
{
    record<type_list<dna5_vector, std::string>, fields<field::SEQ, field::ID>> r{"ACGT"_dna5, "ID1"};

    // ...

    fout.push_back(r);
}
//! [push_back_record]
}

{
//! [push_back_tuple]
sequence_file_output fout{tmp_dir/"my.fasta"};

for(int i = 0; i < 5; i++) // some criteria
{
    std::string id{"test_id"};
    dna5_vector seq{"ACGT"_dna5};

    // ...

    fout.push_back(std::tie(seq, id));
}
//! [push_back_tuple]
}

{
//! [emplace_back]
sequence_file_output fout{tmp_dir/"my.fasta"};

for(int i = 0; i < 5; i++) // some criteria
{
    std::string id{"test_id"};
    dna5_vector seq{"ACGT"_dna5};

    // ...

    fout.emplace_back(seq, id);
}
//! [emplace_back]
}

{
//! [batch_write_2]
sequence_file_output fout{tmp_dir/"my.fasta"};

std::vector<std::tuple<dna5_vector, std::string>> range
{
    { "ACGT"_dna5, "First" },
    { "NATA"_dna5, "2nd" },
    { "GATA"_dna5, "Third" }
}; // a range of "records"

range | fout;
// the same as:
fout = range;
//! [batch_write_2]
}
}
