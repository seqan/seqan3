#include <iostream>
#include <fstream>

#include <seqan3/core/char_operations/predicate.hpp>
#include <seqan3/io/stream/iterator.hpp>
#include <seqan3/range/detail/misc.hpp>
#include <seqan3/range/view/persist.hpp>
#include <seqan3/range/view/single_pass_input.hpp>
#include <seqan3/range/view/take_until.hpp>
#include <seqan3/std/algorithm>
#include <seqan3/std/ranges>

struct write_file_dummy_struct
{
    write_file_dummy_struct()
    {

auto sam_file_raw = R"//![sam_file](
@HD VN:1.6 SO:coordinate
@SQ SN:ref LN:45
r001   99 ref  7 30 8M2I4M1D3M = 37  39 TTAGATAAAGGATACTG *
r003    0 ref  9 30 5S6M       *  0   0 GCCTAAGCTAA       * SA:Z:ref,29,-,6H5M,17,0;
r003 2064 ref 29 17 6H5M       *  0   0 TAGGC             * SA:Z:ref,9,+,5S6M,30,1;
r001  147 ref 37 30 9M         =  7 -39 CAGCGGCAT         * NM:i:1
)//![sam_file]";

        std::ofstream file_stream{"/tmp/my.sam"};
        seqan3::ostreambuf_iterator stream_it{file_stream};

        std::string sam_file{sam_file_raw};
        auto sam_view = sam_file | seqan3::view::single_pass_input;

        ++std::ranges::begin(sam_view); // skip first new line from raw string
        std::ranges::copy(sam_view | seqan3::view::take_until(seqan3::is_char<' '>), stream_it);

        while (begin(sam_view) != end(sam_view))
        {
            stream_it = '\t';
            seqan3::detail::consume(sam_view | seqan3::view::take_until(!seqan3::is_char<' '>));
            std::ranges::copy(sam_view | seqan3::view::take_until(seqan3::is_char<' '>), stream_it);
        }
    }
};

write_file_dummy_struct go{}; // just to have a temporary file accessible such that the snippet does not fail

//! [my_traits]
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/alignment_file/all.hpp>
#include <seqan3/range/container/bitcompressed_vector.hpp>

using namespace seqan3;

struct my_traits : alignment_file_input_default_traits<>
{
    using sequence_alphabet = dna4;                        // instead of dna5

    template <typename alph>
    using sequence_container = bitcompressed_vector<alph>; // must be defined as a template!
};

// ... within main you can then use:

alignment_file_input<my_traits> fin{"/tmp/my.sam"};

// ...
//! [my_traits]

int main()
{

{
//![get_header]
alignment_file_input fin{"/tmp/my.sam"};

// access the header information
debug_stream << fin.header().format_version << std::endl; // 1.6
debug_stream << fin.header().ref_dict << std::endl;       // [(ref,(45,))] (this only works with our debug_stream!)
//![get_header]
}

{
//![construction_from_filename]
alignment_file_input fin{"/tmp/my.sam"}; // SAM format assumed, regular std::ifstream taken as stream
//![construction_from_filename]
}

{
//![construction_from_stream]
std::string input
{
    "@HD\tVN:1.6\tSO:coordinate\n"
    "r001\t99\tref\t7\t30\t8M2I4M1D3M\t=\t37\t39\tTTAGATAAAGGATACTG\t*\n"
};

std::istringstream iss(input);

alignment_file_input fin{iss, format_sam{}};
//              ^ no need to specify the template arguments
//![construction_from_stream]
}

{
std::string input
{
    "@HD\tVN:1.6\tSO:coordinate\n"
    "r001\t99\tref\t7\t30\t8M2I4M1D3M\t=\t37\t39\tTTAGATAAAGGATACTG\t*\n"
};
std::istringstream iss(input);
//![construction_without_automatic_type_deduction]
alignment_file_input<alignment_file_input_default_traits<>,  // The expected format.
                     fields<field::SEQ,                      // The default types; you can adjust this list
                            field::ID,                       // if you don't want to read all this data.
                            field::OFFSET,
                            field::REF_SEQ,
                            field::REF_ID,
                            field::REF_OFFSET,
                            field::ALIGNMENT,
                            field::MAPQ,
                            field::QUAL,
                            field::FLAG,
                            field::MATE,
                            field::TAGS,
                            field::EVALUE,
                            field::BIT_SCORE,
                            field::HEADER_PTR>,
                     type_list<format_sam>, // Which formats are allowed.
                     char>                                 // Value type of the stream.
                     fin{iss, format_sam{}};
//![construction_without_automatic_type_deduction]
}

{
//![reading_range_based_for_loop]
alignment_file_input fin{"/tmp/my.sam"};

for (auto & rec : fin)
{
    debug_stream << "id:  "              << get<field::ID>(rec)         << '\n';
    debug_stream << "read sequence: "    << get<field::SEQ>(rec)        << '\n';
    debug_stream << "mapping position: " << get<field::REF_OFFSET>(rec) << '\n';
    debug_stream << "mapping quality: "  << get<field::MAPQ>(rec)       << '\n';

    // there are more fields read on default
}
//![reading_range_based_for_loop]
}

{
//![reading_move_record]
alignment_file_input fin{"/tmp/my.sam"};

using record_type = typename decltype(fin)::record_type;
std::vector<record_type> records; // store all my records in a vector

for (auto & rec : fin)
    records.push_back(std::move(rec));
//![reading_move_record]
}

{
//![reading_custom_fields]
alignment_file_input fin{"/tmp/my.sam", fields<field::FLAG, field::MAPQ>{}};

for (auto & rec : fin)
{
    debug_stream << "flag:  "            << get<field::FLAG>(rec) << '\n';
    debug_stream << "mapping quality:  " << get<field::MAPQ>(rec) << '\n';
    // get<field::SEQ>(rec) would fail as it was not read
}
//![reading_custom_fields]
}

{
//![reading_structured_bindings]
alignment_file_input fin{"/tmp/my.sam", fields<field::FLAG, field::MAPQ>{}};

for (auto & [flag, mapq] : fin) // the order is the same as specified in fields!
{
    debug_stream << "flag:  "            << flag << '\n';
    debug_stream << "mapping quality:  " << mapq << '\n';
}
//![reading_structured_bindings]
}

{
//![reading_filter]
alignment_file_input fin{"/tmp/my.sam"};

auto minimum_length5_filter = std::view::filter([] (auto const & rec)
{
    return size(get<field::SEQ>(rec)) >= 10;
});

for (auto & rec : fin | minimum_length5_filter) // only records with sequence length >= 10 will "appear"
{
    debug_stream << get<field::ID>(rec) << '\n';
}
//![reading_filter]
}

{
//![begin_and_front]
alignment_file_input fin{"/tmp/my.sam"};
auto it = fin.begin();

// the following are equivalent:
auto & rec0 = *it;
auto & rec1 = fin.front();
// Note: both become invalid after incrementing "it"!
//![begin_and_front]
(void) rec0;
(void) rec1;
}

{
//![front]
alignment_file_input fin{"/tmp/my.sam"};

auto rec = std::move(fin.front()); // rec now stores the data permanently
//![front]
(void) rec;
}

}
