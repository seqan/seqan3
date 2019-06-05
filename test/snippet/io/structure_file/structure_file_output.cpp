#include <set>

#include <seqan3/io/structure_file/output.hpp>
#include <seqan3/io/structure_file/input.hpp>
#include <seqan3/alphabet/nucleotide/all.hpp>
#include <seqan3/alphabet/aminoacid/all.hpp>
#include <seqan3/alphabet/structure/wuss.hpp>
#include <seqan3/alphabet/structure/dssp9.hpp>

using namespace seqan3;

// Code for col_based snippet
//![data_storage]
struct data_storage_t
{
    concatenated_sequences<rna5_vector>         sequences{"AACGUU"_rna5};
    concatenated_sequences<std::string>         ids{std::string{"seq1"}};
    concatenated_sequences<std::vector<wuss51>> structures{".(())."_wuss51};
};

data_storage_t data_storage; // a global or globally used variable in your program
//![data_storage]

auto tmp_dir = std::filesystem::temp_directory_path();

int main()
{

{ // Create a file /tmp/input.dbn for reading.
structure_file_output fout{tmp_dir/"input.dbn"};
fout.emplace_back("GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA"_rna4,
                  "S.cerevisiae_tRNA-PHE M10740/1-73",
                  "(((((((..((((........)))).((((.........)))).....(((((.......))))))))))))."_wuss51);
fout.emplace_back("UUGGAGUACACAACCUGUACACUCUUUC"_rna4,
                  "example",
                  "..(((((..(((...)))..)))))..."_wuss51);

// Create a file /tmp/input_aa.dbn for reading.
structure_file_output fout2{tmp_dir/"input_aa.dbn"};
fout2.emplace_back("ACEWACEW"_aa20,
                  "S.cerevisiae_tRNA-PHE M10740/1-73",
                  "HGEBHHHH"_dssp9);
fout2.emplace_back("ACEWACEWACEWACEW"_aa20,
                  "example",
                  "HGEBHHHHHGEBHHHH"_dssp9);
}

{
//! [temp_param_deduc]
structure_file_output fout{tmp_dir/"output.dbn"}; // Vienna format detected, std::ofstream opened for file
//! [temp_param_deduc]
}

{
//! [write_std_out]
structure_file_output fout{std::cout, format_vienna{}};
//               ^ no need to specify the template arguments

fout.emplace_back("AACGUU"_rna4, "example_id", ".(())."_wuss51); // default order for vienna: SEQ, ID, STRUCTURE
//! [write_std_out]
}

{
//! [iter_by_rec]
structure_file_output fout{tmp_dir/"my.dbn"};

for (int i = 0; i < 10; i++) // ...
{
    std::string id{"test_id"};
    rna5_vector seq{"ACGU"_rna5};
    std::vector<wuss51> structure{".()."_wuss51};

    // ...

    fout.emplace_back(seq, id, structure);        // as individual variables
    // or:
    fout.push_back(std::tie(seq, id, structure)); // as a tuple
}
//! [iter_by_rec]
}

{
//! [write_fields]
structured_rna<rna5, wuss51> sr{'G'_rna5, '.'_wuss51};

structure_file_output fout{tmp_dir/"my.dbn", fields<field::ID, field::STRUCTURED_SEQ>{}};

for (int i = 0; i < 10; i++)// ...
{
    std::string id{"test_id"};
    std::vector<structured_rna<rna5, wuss51>> structured_seq{sr, sr, sr, sr}; // vector of combined data structure

    // ...

    fout.emplace_back(id, structured_seq); // note also that the order the arguments is now different, because
    // or:                                    you specified that ID should be first in the fields template argument
    fout.push_back(std::tie(id, structured_seq));
}
//! [write_fields]
}

{
bool criteria = true;
//! [pass_rec]
structure_file_input fin{tmp_dir/"input.dbn", fields<field::ID, field::SEQ, field::STRUCTURE>{}};
structure_file_output   fout{tmp_dir/"my_wrong.dbn"}; // doesn't have to match the configuration

for (auto & r : fin)
{
    if (criteria) // r fulfills some filter criterium
    {
        fout.push_back(r);
    }
}
//! [pass_rec]
}

{
//! [mult_rec]
structure_file_output fout{tmp_dir/"my.dbn"};

std::vector<std::tuple<rna5_vector, std::string, std::vector<wuss51>>> range
{
    { "ACGT"_rna5, "First", "...."_wuss51 },
    { "NATA"_rna5, "2nd",   "...."_wuss51 },
    { "GATA"_rna5, "Third", "...."_wuss51 }
}; // a range of "records"

fout = range; // will iterate over the records and write them
//! [mult_rec]
}

{
//! [file_conv]
// file format conversion in one line:
structure_file_output{tmp_dir/"output.dbn"} = structure_file_input{tmp_dir/"input.dbn"};

// with structure_file_output as a variable:
structure_file_output fout{tmp_dir/"output.dbn"};
fout = structure_file_input{tmp_dir/"input.dbn"};

// or in pipe notation:
structure_file_input{tmp_dir/"input.dbn"} | structure_file_output{tmp_dir/"output.dbn"};
//! [file_conv]
}

{
//! [pipeline]
structure_file_input my_in{tmp_dir/"input.dbn"};
my_in | view::take(5) | structure_file_output{"output.dbn"};
//! [pipeline]
}

{
// See data_storage snippet above main function for completeness
//! [col_based]
// ... in your file writing function:

structure_file_output fout{tmp_dir/"my.dbn"};

fout = std::tie(data_storage.sequences, data_storage.ids, data_storage.structures);
//! [col_based]
}

{
//! [push_back]
structure_file_output fout{tmp_dir/"my.dbn"};

auto it = fout.begin();

for (int i = 0; i < 10; i++) // ...
{
    std::string id{"test_id"};
    rna5_vector seq{"AGGGUU"_rna5};
    std::vector<wuss51> structure{"..().."_wuss51};

    // ...

    // assign to iterator
    *it = std::tie(seq, id, structure);
    // is the same as:
    fout.push_back(std::tie(seq, id, structure));
}
//! [push_back]
}

{
//! [push_back_2]
structure_file_output fout{tmp_dir/"my.dbn"};

auto it = fout.begin();

for (int i = 0; i < 10; i++) // ...
{
    std::string id{"test_id"};
    rna5_vector seq{"AGGGUU"_rna5};
    std::vector<wuss51> structure{"..().."_wuss51};

    // ...

    fout.push_back(std::tie(seq, id, structure));
}
//! [push_back_2]
}

{
//! [emplace_back]
structure_file_output fout{tmp_dir/"my.dbn"};

auto it = fout.begin();

for (int i = 0; i < 10; i++) // ...
{
    std::string id{"test_id"};
    rna5_vector seq{"AGGGUU"_rna5};
    std::vector<wuss51> structure{"..().."_wuss51};

    // ...

    fout.emplace_back(seq, id, structure);
}
//! [emplace_back]
}

{
//! [equal]
structure_file_output fout{tmp_dir/"my.dbn"};

std::vector<std::tuple<rna5_vector, std::string, std::vector<wuss51>>> range
{
    { "ACGT"_rna5, "First", "...."_wuss51 },
    { "NATA"_rna5, "2nd",   "...."_wuss51 },
    { "GATA"_rna5, "Third", "...."_wuss51 }
}; // a range of "records"

fout = range; // will iterate over the records and write them
//! [equal]
}

{
//! [pipe_func]
structure_file_output fout{tmp_dir/"my.dbn"};

std::vector<std::tuple<rna5_vector, std::string, std::vector<wuss51>>> range
{
    { "ACGT"_rna5, "First", "...."_wuss51 },
    { "NATA"_rna5, "2nd",   "...."_wuss51 },
    { "GATA"_rna5, "Third", "...."_wuss51 }
}; // a range of "records"

range | fout;
// the same as:
fout = range;
//! [pipe_func]
}

}
