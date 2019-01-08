#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/alphabet/quality/all.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/io/stream/debug_stream.hpp>
#include <seqan3/range/container/bitcompressed_vector.hpp>
#include <seqan3/range/view/get.hpp>
#include <seqan3/std/view/filter.hpp>

using namespace seqan3;

struct my_traits : sequence_file_input_default_traits_dna
{
    using sequence_alphabet = dna4;                        // instead of dna5

    template <typename alph>
    using sequence_container = bitcompressed_vector<alph>; // must be defined as a template!
};

// Code for col_read
struct data_storage_t
{
    concatenated_sequences<dna5_vector>  sequences;
    concatenated_sequences<std::string>  ids;
};

data_storage_t data_storage; // a global or globally used variable in your program

int main()
{

{
// First make a /tmp/my.fasta file.
sequence_file_output fout{"/tmp/my.fasta"};
fout.emplace_back("ACGT"_dna4, "TEST1");
fout.emplace_back("AGGCTGA"_dna4, "Test2");
fout.emplace_back("GGAGTATAATATATATATATATAT"_dna4, "Test3");
}

{
// Also make a /tmp/my_seq_qual.fasta file.
// TODO: Uncomment once @smehringer fixes bug.
// sequence_file_output fout{"/tmp/my.fastq", fields<field::ID, field::SEQ_QUAL>{}};
// std::vector<qualified<dna4, phred42>> my_qualified{qualified<dna4, phred42>{'A'_dna4, phred42{7}},
//                                                    qualified<dna4, phred42>{'C'_dna4, phred42{6}},
//                                                    qualified<dna4, phred42>{'G'_dna4, phred42{5}},
//                                                    qualified<dna4, phred42>{'T'_dna4, phred42{4}}};
// fout.emplace_back("TEST1", my_qualified);
// fout.emplace_back("Test2", my_qualified);
// fout.emplace_back("Test3", my_qualified);
}

{
#if 0
//! [trait_overwrite]
struct my_traits : sequence_file_input_default_traits_dna
{
    using sequence_alphabet = dna4;                        // instead of dna5

    template <typename alph>
    using sequence_container = bitcompressed_vector<alph>; // must be defined as a template!
};

sequence_file_input<my_traits> fin{"/tmp/my.fasta"};

//...
//! [trait_overwrite]
#endif

// TODO: Uncomment once sequence_container_container concepts are fixed by @h-2.
// [[maybe_unused]] sequence_file_input<my_traits> fin{"/tmp/my.fasta"};
}

{
//! [template_deduction]
sequence_file_input fin{"/tmp/my.fasta"}; // FastA with DNA sequences assumed, regular std::ifstream taken as stream
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
sequence_file_input<sequence_file_input_default_traits_aa> fin{"/tmp/my.fasta"};
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
sequence_file_input fin{"/tmp/my.fasta"};

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
sequence_file_input fin{"/tmp/my.fasta"};

using record_type = typename decltype(fin)::record_type;
std::vector<record_type> records;

for (auto & rec : fin)
    records.push_back(std::move(rec));
//! [auto_ref]
}

{
//! [decomposed]
sequence_file_input fin{"/tmp/my.fasta"};

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
// TODO: Uncomment once @smehringer fixes bug.
// sequence_file_input fin{"/tmp/my.fastq", fields<field::ID, field::SEQ_QUAL>{}};
//
// for (auto & [ id, seq_qual ] : fin) // note that the order is now different, "id" comes first, because it was specified first
// {
//     debug_stream << "ID:  " << id << '\n';
//     // sequence and qualities are part of the same vector, of type std::vector<dna5q>
//     debug_stream << "SEQ: "  << (seq_qual | view::get<0>) << '\n'; // sequence string is extracted
//     debug_stream << "QUAL: " << (seq_qual | view::get<1>) << '\n'; // quality string is extracted
// }
//! [custom_fields]
}

{
//! [file_view]
sequence_file_input fin{"/tmp/my.fasta"};

auto minimum_length5_filter = view::filter([] (auto const & rec)
{
    return size(get<field::SEQ>(rec)) >= 5;
});

for (auto & rec : fin | minimum_length5_filter) // only record with sequence length >= 5 will "appear"
{
    debug_stream << "IDs of seq_length >= 5: " << get<field::ID>(rec) << '\n';
    // ...
}
//! [file_view]
}

{
#if 0
//! [col_read]
struct data_storage_t
{
    concatenated_sequences<dna5_vector>  sequences;
    concatenated_sequences<std::string>  ids;
};

data_storage_t data_storage; // a global or globally used variable in your program

// ... in your file reading function:

sequence_file_input fin{"/tmp/my.fasta"};

data_storage.sequences = std::move(get<field::SEQ>(fin)); // we move the buffer directly into our storage
data_storage.ids       = std::move(get<field::ID>(fin)); // we move the buffer directly into our storage
//! [col_read]
#endif
data_storage_t data_storage; // a global or globally used variable in your program

// ... in your file reading function:

sequence_file_input fin{"/tmp/my.fasta"};

data_storage.sequences = std::move(get<field::SEQ>(fin)); // we move the buffer directly into our storage
data_storage.ids       = std::move(get<field::ID>(fin)); // we move the buffer directly into our storage
}

{
//! [return_record]
sequence_file_input fin{"/tmp/my.fasta"};
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
sequence_file_input fin{"/tmp/my.fasta"};

auto rec0 = std::move(fin.front());
//! [record_move]
(void) rec0;
}
}
