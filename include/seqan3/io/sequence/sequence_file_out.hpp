#pragma once

#include <string>
#include <variant>
#include <vector>

namespace seqan3
{

// ==================================================================
// sequence_file_out_traits
// ==================================================================

  
template <typename t>
concept bool sequence_file_out_traits_concept = requires (t v)
{
    t::stream_type;
    t::valid_formats;
    requires detail::meets_sequence_file_format_concept<0, typename t::valid_formats>();
    
    t::valid_compression_formats;
};


struct sequence_file_out_default_traits
{
    using stream_type = std::ofstream;
    using valid_formats = std::variant<sequence_file_format_fasta<stream_type>,
                                       sequence_file_format_fastq<stream_type>,
                                       sequence_file_format_embl<stream_type>,
                                       sequence_file_format_genbank<stream_type>,
                                       sequence_file_format_raw<stream_type>>;
    static constexpr std::vector<std::pair<std::string, void>> valid_compression_formats{};
};



// ==================================================================
// sequence_file_out
// ==================================================================

template <typename sequence_file_out_traits = sequence_file_out_default_traits>
    requires sequence_file_traits_concept<sequence_file_out_traits>
class sequence_file_out : protected sequence_file_base<sequence_file_out_traits>
{
public:
    /* types */
    
    /* constructors */
    
    // constructor with arg
    sequence_file_out(std::string const & _file_name) :
      sequence_file_base{_file_name}
    {} 

    // copy construction and assignment are deleted
    // implicitly because we don't want multiple access to file
    sequence_file_out = delete;
    sequence_file_out(sequence_file_out const &) = delete;
    sequence_file_out & operator=(sequence_file_out const &) = delete;

    // move construction and assignment are defaulted
    sequence_file_out(sequence_file_out &&) = default;
    sequence_file_out & operator=(sequence_file_out &&) = default;

    //TODO make the requirements stricter
    template <typename seq_type, typename id_type, typename qual_type>
        requires container_concept<seq_type> &&
                 container_concept<id_type> &&
                 container_concept<qual_type>
    void write(seq_type && seq, id_type && id, qual_type && qual);

    template <typename seqs_type, typename ids_type, typename quals_type, size_t max_records = 0>
        requires container_concept<typename seqs_type::value> &&
                 container_concept<typename ids_type::value> &&
                 container_concept<typename quals_type::value>
    void write(seqs_type && seqs, ids_type && ids, quals_type && quals);


protected:
    ~sequence_file_out() = default;


};



template <typename seq_type, typename id_type, typename qual_type>
    requires container_concept<seq_type> &&
             container_concept<id_type> &&
             container_concept<qual_type>
inline void sequence_file_out::write(seq_type && seq, id_type && id = std::string{}, qual_type && qual = std::string{})
{
    assert(!format.valueless_by_exception);
    std::visit([&] (sequence_file_out_format_concept & f) { f->write(seq, id, qual, stream, options); }, format);
}

template <typename seqs_type, typename ids_type, typename quals_type, size_t num_records = 0>
    requires container_concept<typename seqs_type::value> &&
             container_concept<typename ids_type::value> &&
             container_concept<typename quals_type::value>
inline void sequence_file_out::write(seqs_type  && seqs,
                              ids_type   && ids   = std::vector<std::string>{},
                              quals_type && quals = std::vector<std::string>{})
{
    assert(!format.valueless_by_exception);
    std::visit([&] (sequence_file_out_format_concept & f)
    {
        f->write<num_records>(seqs, ids, quals, stream, options);
    }, format);
}



} // namespace seqan