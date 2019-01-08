#include <iostream>

#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/io/alignment_file/all.hpp>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/range/view/persist.hpp>

using namespace seqan3;

int main()
{

auto tmp_dir = std::filesystem::temp_directory_path();

// The output snippets will create a /tmp/my.sam
{
//! [filename_construction]
alignment_file_output fout{tmp_dir/"my.sam"}; // SAM format detected, std::ofstream opened for file
//! [filename_construction]
}

{
//! [format_construction]
// no need to specify the template arguments <...> for format specialization:
alignment_file_output fout{tmp_dir/"my.sam", fields<field::MAPQ>{}};
//! [format_construction]
}

{
//! [write_range]
alignment_file_output fout{tmp_dir/"my.sam"};

std::vector<std::tuple<dna5_vector, std::string>> range
{
    { "ACGT"_dna5, "First" },
    { "NATA"_dna5, "2nd" },
    { "GATA"_dna5, "Third" }
}; // a range of "records"

fout = range; // will iterate over the records and write them
//! [write_range]
}

{
//! [set_header]
alignment_file_output fout{tmp_dir/"my.sam"};

// add information to the header of the file.
fout.header().comments.push_back("This is a comment");
//! [set_header]
}

{
//! [custom_fields]
//! \todo Uncomment once seqan3::alignment_file_input is implemented.
// alignment_file_input  fin{"input.sam", fields<field::REF_OFFSET, field::FLAG>{}};
// alignment_file_output fout{"output.sam"}; // doesn't have to match the configuration

// for (auto & r : fin)
// {
//     fout.push_back(r); // copy all the records.
// }
//! [custom_fields]
}

// Create a /tmp/input.fasta for the snippet below
{
    sequence_file_output fout{tmp_dir/"input.fasta"};
    fout.emplace_back("ACGT"_dna4, "TEST1");
    fout.emplace_back("AGGCTGA"_dna4, "Test2");
    fout.emplace_back("ACTGA"_dna4, "Test2");
}
// Create a /tmp/input.sam for the snippet below
{
    alignment_file_output fout{tmp_dir/"input.sam"};

    std::vector<std::tuple<dna5_vector, std::string>> range
    {
        { "ACGT"_dna5, "First" },
        { "NATA"_dna5, "2nd" },
        { "GATA"_dna5, "Third" }
    }; // a range of "records"

    fout = range; // will iterate over the records and write them
}

{
//! [input_range]
//! \todo Uncomment once seqan3::alignment_file_input is implemented.
// // copying a file in one line:
// alignment_file_output{tmp_dir/"output.sam"} = alignment_file_input{tmp_dir/"input.sam"};

// with alignment_file_output as a variable:
// alignment_file_output fout{tmp_dir/"output.sam"};
// alignment_file_input fin{tmp_dir/"input.sam"};
// fout = fin;

// // or in pipe notation:
// alignment_file_input{tmp_dir/"input.sam"} | alignment_file_output{tmp_dir/"output.sam"};
//! [input_range]
}

{
//! [io_pipeline]
//! \todo Uncomment once seqan3::alignment_file_input is implemented.
// alignment_file_input{tmp_dir/"input.sam"} | view::persist
//                                           | ranges::view::take(5) // take only the first 5 records
//                                           | alignment_file_output{tmp_dir/"output.sam"};
//! [io_pipeline]
}

{
//! [range]
alignment_file_output fout{tmp_dir/"my.sam"};

std::vector<std::tuple<dna5_vector, std::string>> range
{
    { "ACGT"_dna5, "First" },
    { "NATA"_dna5, "2nd" },
    { "GATA"_dna5, "Third" }
}; // a range of "records"

fout = range; // will iterate over the records and write them
//! [range]
}
}
