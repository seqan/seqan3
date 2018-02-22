#pragma once

#include <vector>
#include <string>

#if 0
//TODO(rrahn): this is a prototype and needs more refinement, disabling for now
//!\cond

namespace seqan3
{

class sequence_file_format_fasta
{
public:
    sequence_file_format_fasta() :
    {}

    // rule of six constructors deleted
    sequence_file_format_fasta() = delete;
    sequence_file_format_fasta(sequence_file_in const &) = delete;
    sequence_file_format_fasta & operator=(sequence_file_in const &) = delete;
    sequence_file_format_fasta(sequence_file_in &&) = delete;
    sequence_file_format_fasta & operator=(sequence_file_in &&) = delete;

    static std::vector<std::string> file_extensions
    {
        { "fasta" },
        { "fa" }
        // ...
    };

    //TODO make the requirements stricter
    template <typename sequence_type,
              typename id_type,
              typename qual_type,
              typename stream_type,
              typename options_type>
        requires container_concept<sequence_type> &&
                 container_concept<id_type> &&
                 container_concept<qual_type>
    void read(sequence_type && seq,
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

//!\endcond
#endif
