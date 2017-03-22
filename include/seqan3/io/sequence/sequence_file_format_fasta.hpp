#pragma once

#include <limits>
#include <vector>
#include <string>
#include <tuple>
#include "for_test.hpp"

namespace seqan3
{

struct sequence_file_format_fasta
{
public:
    // sequence_file_format_fasta(std::string file_name)
    // {
    //     stream.open(_file_name, std::ios::in);
    // }

    // rule of six constructors deleted
    // sequence_file_format_fasta() = default;
    // sequence_file_format_fasta(sequence_file_format_fasta const & ) = default;
    // sequence_file_format_fasta & operator=(sequence_file_format_fasta const &) = default;
    // sequence_file_format_fasta(sequence_file_format_fasta &&) = default;
    // sequence_file_format_fasta & operator=(sequence_file_format_fasta &&) = default;

    static inline std::vector<std::string> file_extensions {{"fasta"},{"fa"}};

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
              typename options_type,
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

/** implementations **/
template <typename sequence_type,
          typename id_type,
          typename qual_type,
          typename stream_type,
          typename options_type>
requires container_concept<sequence_type> &&
         container_concept<id_type> &&
         container_concept<qual_type>
void sequence_file_format_fasta::read(sequence_type && seq,
                                      id_type && id,
                                      qual_type && qual,
                                      stream_type & stream,
                                      options_type const & options)
{
    skip_until(stream, '>');    // forward to the next '>'
    skip_one(stream, '>');    // assert and skip '>'
    // read_line(id, options.id_filter, stream);                        // read Fasta id
    // read_until(seq, options.sequence_filter, stream, '>'); // read Fasta sequence
}
template <typename seqs_type,
          typename ids_type,
          typename quals_type,
          typename stream_type,
          typename options_type,
          size_t max_records = 0>
requires container_concept<typename seqs_type::value> &&
         container_concept<typename ids_type::value> &&
         container_concept<typename quals_type::value>
void sequence_file_format_fasta::read(seqs_type && seqs,
                                      ids_type && ids,
                                      quals_type && quals,
                                      stream_type & stream,
                                      options_type const & options)
{
    typedef typename seqs_type::value_type seq_type;
    typedef typename ids_type::value_type id_type;
    typedef typename quals_type::value_type qual_type;
    uint32_t    records_to_read = max_records;
    if (max_records == 0)
    {
        std::numeric_limits<uint32_t>::max();
    }
    for (; !stream.at_end() && records_to_read > 0; --records_to_read)
    {
        // read(seqs[last], temp_id, temp_qual, stream, options);
        // seqs.push_back(new_seq);
        // ids.push_back(new_id);
        // quals.push_back(new_qual);
    }
}
} // namespace seqan3