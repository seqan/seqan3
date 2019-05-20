#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/all.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/range/container/bitcompressed_vector.hpp>
#include <seqan3/range/view/get.hpp>
#include <seqan3/std/ranges>

using namespace seqan3;

auto tmp_dir = std::filesystem::temp_directory_path();

//! [data_storage]
struct data_storage_t
{
    concatenated_sequences<dna5_vector>  sequences;
    concatenated_sequences<std::string>  ids;
};

data_storage_t data_storage; // a global or globally used variable in your program
//! [data_storage]

// Todo: change std::vector to bitcompressed_vector once #764 is fixed
//! [trait_overwrite]
struct my_traits : sequence_file_input_default_traits_dna
{
    using sequence_alphabet = dna4;               // instead of dna5

    template <typename alph>
    using sequence_container = std::vector<alph>; // must be defined as a template!
};

// [...] main

    sequence_file_input<my_traits> fin{tmp_dir/"my.fasta"};

// [...]
//! [trait_overwrite]

int main()
{
{
// Create a /tmp/my.fasta file.
sequence_file_output fout{tmp_dir/"my.fasta"};
fout.emplace_back("ACGT"_dna4, "TEST1");
fout.emplace_back("AGGCTGA"_dna4, "Test2");
fout.emplace_back("GGAGTATAATATATATATATATAT"_dna4, "Test3");
}

{
// Create a /tmp/my_seq_qual.fasta file.
sequence_file_output fout{tmp_dir/"my.fastq", fields<field::ID, field::SEQ_QUAL>{}};
std::vector<qualified<dna4, phred42>> my_qualified{{'A'_dna4, '@'_phred42},
                                                   {'C'_dna4, '@'_phred42},
                                                   {'G'_dna4, '@'_phred42},
                                                   {'T'_dna4, '@'_phred42}};
fout.emplace_back("TEST1", my_qualified);
fout.emplace_back("Test2", my_qualified);
fout.emplace_back("Test3", my_qualified);
}

{
//! [template_deduction]
sequence_file_input fin{tmp_dir/"my.fasta"}; // FastA with DNA sequences assumed, regular std::ifstream taken as stream
//! [template_deduction]
(void) fin;
}

{
//! [istringstream]
std::string input
{
    "> TEST1\n"
    "ACGT\n"
    "> Test2\n"
    "AGGCTGN\n"
    "> Test3\n"
    "GGAGTATAATATATATATATATAT\n"
};

std::istringstream iss(input);

sequence_file_input fin{std::move(iss), sequence_file_format_fasta{}};
//              ^ no need to specify the template arguments
//! [istringstream]
(void) fin;
}

{
//! [aminoacid]
sequence_file_input<sequence_file_input_default_traits_aa> fin{tmp_dir/"my.fasta"};
//! [aminoacid]
(void) fin;
}

{
//! [template_specification]
// ... input had amino acid sequences
std::string input
{
    "> TEST1\n"
    "FQTWE\n"
    "> Test2\n"
    "KYRTW\n"
    "> Test3\n"
    "EEYQTWEEFARAAEKLYLTDPMKV\n"
};
std::istringstream iss(input);

sequence_file_input<sequence_file_input_default_traits_aa /*Use amino acid traits here*/,
                    fields<field::SEQ, field::ID, field::QUAL>,
                    type_list<sequence_file_format_fasta>, char> fin{iss, sequence_file_format_fasta{}};
//! [template_specification]
(void) fin;
}

{
//! [record_iter]
sequence_file_input fin{tmp_dir/"my.fasta"};

for (auto & rec : fin)
{
    debug_stream << "ID:  " << get<field::ID>(rec) << '\n';
    debug_stream << "SEQ: " << get<field::SEQ>(rec) << '\n';
    // a quality field also exists, but is not printed, because we know it's empty for FastA files.
}
//! [record_iter]
}

{
//! [auto_ref]
sequence_file_input fin{tmp_dir/"my.fasta"};

using record_type = typename decltype(fin)::record_type;
std::vector<record_type> records;

for (auto & rec : fin)
    records.push_back(std::move(rec));
//! [auto_ref]
}

{
//! [decomposed]
sequence_file_input fin{tmp_dir/"my.fasta"};

for (auto & [ seq, id, qual ] : fin)
{
    debug_stream << "ID:  " << id << '\n';
    debug_stream << "SEQ: " << seq << '\n';
    debug_stream << "EMPTY QUAL." << qual << '\n'; // qual is empty for FastA files
}
//! [decomposed]
}

{
//! [custom_fields]
sequence_file_input fin{tmp_dir/"my.fastq", fields<field::ID, field::SEQ_QUAL>{}};

for (auto & [ id, seq_qual ] : fin) // the order is now different, "id" comes first, because it was specified first
{
    debug_stream << "ID:  " << id << '\n';
    // sequence and qualities are part of the same vector, of type std::vector<dna5q>
    debug_stream << "SEQ: "  << (seq_qual | view::get<0>) << '\n'; // sequence string is extracted
    debug_stream << "QUAL: " << (seq_qual | view::get<1>) << '\n'; // quality string is extracted
}
//! [custom_fields]
}

{
//! [file_view]
sequence_file_input fin{tmp_dir/"my.fasta"};

auto minimum_length5_filter = std::view::filter([] (auto const & rec)
{
    return std::ranges::size(get<field::SEQ>(rec)) >= 5;
});

for (auto & rec : fin | minimum_length5_filter) // only record with sequence length >= 5 will "appear"
{
    debug_stream << "IDs of seq_length >= 5: " << get<field::ID>(rec) << '\n';
    // ...
}
//! [file_view]
}

{
// see snippet data_storage above main function for completeness
//! [col_read]
// ... in your file reading function:

sequence_file_input fin{tmp_dir/"my.fasta"};

data_storage.sequences = std::move(get<field::SEQ>(fin)); // we move the buffer directly into our storage
data_storage.ids       = std::move(get<field::ID>(fin));  // we move the buffer directly into our storage
//! [col_read]
}

{
//! [return_record]
sequence_file_input fin{tmp_dir/"my.fasta"};
auto it = begin(fin);

// the following are equivalent:
auto & rec0 = *it;
auto & rec1 = fin.front();

// Note: rec0 and rec1 are references and become invalid after incrementing "it"!
//! [return_record]
(void) rec0;
(void) rec1;
}

{
//! [record_move]
sequence_file_input fin{tmp_dir/"my.fasta"};

auto rec0 = std::move(fin.front());
//! [record_move]
(void) rec0;
}
}
