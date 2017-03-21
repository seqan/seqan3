#pragma once

#include <vector>
#include <string>

namespace seqan3
{

class seq_file_in_format_fasta
{
public:
    seq_file_in_format_fasta() :
    {}

    // rule of six constructors deleted
    seq_file_in_format_fasta() = delete;
    seq_file_in_format_fasta(seq_file_in const &) = delete;
    seq_file_in_format_fasta & operator=(seq_file_in const &) = delete;
    seq_file_in_format_fasta(seq_file_in &&) = delete;
    seq_file_in_format_fasta & operator=(seq_file_in &&) = delete;

    static std::vector<std::string> file_extensions
    {
        { "fasta" },
        { "fa" }
        // ...
    };

    //TODO make the requirements stricter
    template <typename seq_type,
              typename id_type,
              typename qual_type,
              typename stream_type,
              typename options_type>
        requires container_concept<seq_type> &&
                 container_concept<id_type> &&
                 container_concept<qual_type>
    void read(seq_type && seq,
              id_type && id,
              qual_type && qual,
              stream_type & stream,
              options_type const & options);

    template <typename seqs_type,
              typename ids_type,
              typename quals_type,
              typename stream_type,
              typename options_type
              size_t max_records = 0>
        requires container_concept<typename seqs_type::value> &&
                 container_concept<typename ids_type::value> &&
                 container_concept<typename quals_type::value>
    void read(seqs_type && seqs,
              ids_type && ids,
              quals_type && quals,
              stream_type & stream,
              options_type const & options);
};

static_assert(align_file_out_format_concept<align_file_out_format_sam>,
              "align_file_out_format_sam does not satisfy align_file_out_format_concept");

/** implementations **/


} // namespace seqan3