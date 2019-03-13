#include <iostream>
#include <fstream>

//![include_ranges_chunk]
#include <range/v3/view/chunk.hpp>
//![include_ranges_chunk]
#include <range/v3/view/filter.hpp>
//![include_ranges_take]
#include <range/v3/view/take.hpp>
//![include_ranges_take]

#include <seqan3/io/stream/parse_condition.hpp>
//![include]
#include <seqan3/io/sequence_file/all.hpp>
//![include]
//![include_debug_stream]
#include <seqan3/io/stream/debug_stream.hpp>
//![include_debug_stream]
#include <seqan3/range/detail/misc.hpp>
#include <seqan3/range/view/persist.hpp>
#include <seqan3/range/view/single_pass_input.hpp>
//![include_filter]
#include <seqan3/std/view/filter.hpp>
//![include_filter]

struct write_file_dummy_struct
{
    write_file_dummy_struct()
    {

auto file_raw = R"//![fastq_file](
@seq1
AGCTAGCAGCGATCG
+
IIIIIHIIIIIIIII
@seq2
CGATCGATC
+
IIIIIIIII
@seq3
AGCGATCGAGGAATATAT
+
IIIIHHGIIIIHHGIIIH
)//![fastq_file]";

        std::ofstream file{"/tmp/my.fastq"};
        std::string str{file_raw};
        file << str.substr(1); // skip first newline

        std::ofstream file2{"/tmp/my.qq"};
        file2 << str.substr(1); // skip first newline

        std::ofstream file3{"/tmp/my.fasta"};
        file3 << ">seq1\nAVAV\n>seq2\nAVAVA\n";
    }
};

write_file_dummy_struct go{}; // write file

using namespace seqan3;

int main()
{

{
//![file_extensions]
debug_stream << sequence_file_format_fastq::file_extensions << std::endl; // prints [fastq,fq]
//![file_extensions]

//![modify_file_extensions]
sequence_file_format_fastq::file_extensions.push_back("qq");
sequence_file_input fin{"/tmp/my.qq"}; // detects FASTQ format
//![modify_file_extensions]
}

{
/*
//![construct_from_cin]
sequence_file_input fin{std::cin, sequence_file_format_fasta{}};
//![construct_from_cin]
*/
}

{
//![amino_acid_type_trait]
sequence_file_input<sequence_file_input_default_traits_aa> fin{"/tmp/my.fasta"};
//![amino_acid_type_trait]
}

{
//![custom_fields]
sequence_file_input fin{"/tmp/my.fastq", fields<field::SEQ>{}};

for (auto & rec : fin)
{
    debug_stream << "SEQ:  "  << get<field::SEQ>(rec) << '\n';
    // won't work: get<field::ID>(rec)
    // won't work: get<field::QUAL>(rec)
}

// or with structured bindings:

for (auto & [ seq ] : fin)
{
    debug_stream << "SEQ:  "  << seq << '\n';
}
//![custom_fields]
}

{
//![read_in_batches]
sequence_file_input fin{"/tmp/my.fastq"};

for (auto && records : fin | ranges::view::chunk(10)) // && is important! because view::chunk returns temporaries
{
    // `records` contains 10 elements (or less at the end)
    debug_stream << "Taking the next 10 sequences:\n";
    debug_stream << "ID:  " << get<field::ID>(*records.begin()) << '\n'; // prints first ID in batch
}
//![read_in_batches]
}

{
//![writing_custom_fields]
sequence_file_output fout{"/tmp/output.fasta", fields<field::ID, field::SEQ>{}}; // now the order is changed

std::string id{"test_id"};
dna5_vector seq{"ACGT"_dna5};

fout.emplace_back(id, seq);          // first ID then SEQ !
//![writing_custom_fields]
}

{
//![piping_in_out]
sequence_file_input fin{"/tmp/my.fastq"};
sequence_file_output fout{"/tmp/output.fastq"};

// the following are equivalent:
fin | fout;

fout = fin;

sequence_file_output{"/tmp/output.fastq"} = sequence_file_input{"/tmp/my.fastq"};
//![piping_in_out]
}

{
//![file_conversion]
sequence_file_output{"/tmp/output.fasta"} = sequence_file_input{"/tmp/my.fastq"};
//![file_conversion]
}

}
