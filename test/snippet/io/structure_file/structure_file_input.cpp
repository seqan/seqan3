#include <iostream>

#include <seqan3/alphabet/nucleotide/rna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/structure_file/input.hpp>
#include <seqan3/io/structure_file/output.hpp>
#include <seqan3/range/view/get.hpp>
#include <seqan3/std/ranges>

using namespace seqan3;

// For col_read snippet:
//![data_storage]
struct data_storage_t
{
    concatenated_sequences<rna5_vector>         sequences;
    concatenated_sequences<std::string>         ids;
    concatenated_sequences<std::vector<wuss51>> structures;
};

data_storage_t data_storage; // a global or globally used variable in your program
//![data_storage]

auto tmp_dir = std::filesystem::temp_directory_path();

int main()
{
// First, make /tmp/input.dbn and /tmp/input_aa.dbn
{
structure_file_output fout{tmp_dir/"input.dbn"};
fout.emplace_back("GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA"_rna4,
                  "S.cerevisiae_tRNA-PHE M10740/1-73",
                  "(((((((..((((........)))).((((.........)))).....(((((.......))))))))))))."_wuss51);
fout.emplace_back("UUGGAGUACACAACCUGUACACUCUUUC"_rna4,
                  "example",
                  "..(((((..(((...)))..)))))..."_wuss51);

structure_file_output fout2{tmp_dir/"input_aa.dbn"};
fout2.emplace_back("ACEWACEW"_aa20,
                  "S.cerevisiae_tRNA-PHE M10740/1-73",
                  "HGEBHHHH"_dssp9);
fout2.emplace_back("ACEWACEWACEWACEW"_aa20,
                  "example",
                  "HGEBHHHHHGEBHHHH"_dssp9);
}
{
//! [auto_temp_deduc]
structure_file_input sf{tmp_dir/"input.dbn"}; // Vienna with RNA sequences assumed, use regular std::ifstream as stream
//! [auto_temp_deduc]
}

{
//! [stringstream_read]
std::string const input
{
    ">S.cerevisiae_tRNA-PHE M10740/1-73\n"
    "GCGGAUUUAGCUCAGUUGGGAGAGCGCCAGACUGAAGAUUUGGAGGUCCUGUGUUCGAUCCACAGAAUUCGCA\n"
    "(((((((..((((........)))).((((.........)))).....(((((.......)))))))))))). (-17.50)\n"
    "> example\n"
    "UUGGAGUACACAACCUGUACACUCUUUC\n"
    "..(((((..(((...)))..)))))... (-3.71)\n"
};

std::istringstream iss(input);

structure_file_input fin{iss, format_vienna{}};
//                 ^ no need to specify the template arguments
//! [stringstream_read]
}

{
//! [arg_spec]
structure_file_input<structure_file_input_default_traits_rna> fin{tmp_dir/"input.dbn"};
//! [arg_spec]
}

{
std::string const input
{
    ">S.cerevisiae_tRNA-PHE M10740/1-73\n"
    "ACEW\n"
    "HBEG\n"
    "> example\n"
    "ACEWACEW\n"
    "HGEBHHHH\n"
};

//! [trait_def]
// ... input had amino acid sequences
std::istringstream iss(input);

structure_file_input<structure_file_input_default_traits_aa,
                     fields<field::SEQ, field::ID, field::STRUCTURE>,
                     type_list<format_vienna>,
                     char> fin{iss, format_vienna{}};
//! [trait_def]
}

{
//! [record_iter]
structure_file_input<structure_file_input_default_traits_aa,
                     fields<field::SEQ, field::ID, field::STRUCTURE>,
                     type_list<format_vienna>> fin{tmp_dir/"input_aa.dbn"};

for (auto & rec : fin)
{
    debug_stream << "ID: " << get<field::ID>(rec) << '\n';
    debug_stream << "SEQ: " << (get<field::SEQ>(rec) | view::to_char) << '\n'; // sequence is converted to char on-the-fly
    debug_stream << "STRUCTURE: " << (get<field::STRUCTURE>(rec) | view::to_char) << '\n';
}
//! [record_iter]
}

{
//! [data_out]
structure_file_input fin{tmp_dir/"input.dbn"};

using record_type = typename decltype(fin)::record_type;
std::vector<record_type> records;

for (auto & rec : fin)
    records.push_back(std::move(rec));
//! [data_out]
}

{
//! [structured_bindings]
structure_file_input<structure_file_input_default_traits_aa,
                     fields<field::SEQ, field::ID, field::STRUCTURE>,
                     type_list<format_vienna>> fin{tmp_dir/"input_aa.dbn"};

for (auto & [ seq, id, structure ] : fin)
{
    debug_stream << "ID: " << id << '\n';
    debug_stream << "SEQ: " << (seq | view::to_char) << '\n'; // sequence is converted to char on-the-fly
    debug_stream << "STRUCTURE: " << (structure | view::to_char) << '\n';
}
//! [structured_bindings]
}

{
//! [skip_fields]
structure_file_input fin{tmp_dir/"input.dbn", fields<field::ID, field::STRUCTURED_SEQ>{}};

// note that the order is now different, "id" comes first, because it was specified first
for (auto & [ id, structured_seq ] : fin)
{
    debug_stream << "ID: " << id << '\n';
    // sequence and structure are part of the same vector, of type std::vector<structured_rna<rna5, wuss51>>
    debug_stream << "SEQ: "  << (structured_seq | view::get<0> | view::to_char) << '\n'; // sequence string is extracted and converted to char on-the-fly
    debug_stream << "STRUCTURE: " << (structured_seq | view::get<1> | view::to_char) << '\n'; // structure string is extracted and converted to char on-the-fly
}
//! [skip_fields]
}

{
//! [filter_criteria]
structure_file_input fin{tmp_dir/"input.dbn"};

auto minimum_length5_filter = std::view::filter([] (auto const & rec)
{
    return size(get<field::SEQ>(rec)) >= 5;
});

for (auto & rec : fin | minimum_length5_filter) // only record with sequence length >= 5 will "appear"
{
    debug_stream << (get<0>(rec) | view::to_char) << '\n';
}
//! [filter_criteria]
}

// For col_read snippet
structure_file_input fin{tmp_dir/"input.dbn"};

data_storage.sequences  = std::move(get<field::SEQ>(fin));       // we move the buffer directly into our storage
data_storage.ids        = std::move(get<field::ID>(fin));        // we move the buffer directly into our storage
data_storage.structures = std::move(get<field::STRUCTURE>(fin)); // we move the buffer directly into our storage

{
//! [ref_return]
structure_file_input fin{tmp_dir/"input.dbn"};
auto it = begin(fin);

// the following are equivalent:
auto & rec0 = *it;
auto & rec1 = fin.front();

// both become invalid after incrementing "it"!
//! [ref_return]
(void) rec0;
(void) rec1;
}

{
//! [move]
structure_file_input fin{tmp_dir/"input.dbn"};

auto rec0 = std::move(fin.front());
//! [move]
(void) rec0;
}

{
// See data_storage snippet above main function for completeness
//! [col_read]
// ... in your file reading function:

structure_file_input fin{tmp_dir/"input.dbn"};

data_storage.sequences  = std::move(get<field::SEQ>(fin));       // we move the buffer directly into our storage
data_storage.ids        = std::move(get<field::ID>(fin));        // we move the buffer directly into our storage
data_storage.structures = std::move(get<field::STRUCTURE>(fin)); // we move the buffer directly into our storage
//! [col_read]
}

}
