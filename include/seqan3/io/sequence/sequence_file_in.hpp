#pragma once

#include <string>
#include <variant>
#include <vector>

namespace seqan3
{

// ==================================================================
// seq_file_in_traits
// ==================================================================

template <typename t>
concept bool seq_file_in_traits_concept = requires (t v)
{
    t::stream_type;
    t::valid_formats;
    requires detail::meets_concept_sequence_file_format<0, typename t::valid_formats>();

    t::valid_compression_formats;
};

struct seq_file_in_default_traits
{
    using stream_type = std::ofstream;
    using valid_formats = std::variant<seq_file_in_format_fasta<stream_type>,
                                       seq_file_in_format_fastq<stream_type>,
                                       seq_file_in_format_embl<stream_type>,
                                       seq_file_in_format_genbank<stream_type>,
                                       seq_file_in_format_raw<stream_type>>;
    static constexpr std::vector<std::pair<std::string, void>> valid_compression_formats{};
};

// ==================================================================
// seq_file_in
// ==================================================================

template <typename seq_file_in_traits = seq_file_in_traits_default>
    requires seq_file_in_traits_concept<seq_file_in_traits>
class seq_file_in
{
public:
    using seq_file_in_traits::stream_type;               // e.g. std::ostream concept
    using seq_file_in_traits::valid_compression_formats; // = std::vector<std::pair<std::string, t>>
    using seq_file_in_traits::valid_formats;             // = std::variant<gzip_compressor, seq_file_in_format_bam>

    // constructor with arg
    seq_file_in(std::string const & _file_name);

    // copy construction and assignment are deleted
    // implicitly because we don't want multiple access to file
    seq_file_in(seq_file_in const &) = delete;
    seq_file_in & operator=(seq_file_in const &) = delete;

    // move construction and assignment are defaulted
    seq_file_in(seq_file_in &&) = default;
    seq_file_in & operator=(seq_file_in &&) = default;

    //TODO make the requirements stricter
    template <typename seq_type, typename id_type, typename qual_type>
        requires container_concept<seq_type> &&
                 container_concept<id_type> &&
                 container_concept<qual_type>
    void read(seq_type && seq, id_type && id, qual_type && qual);

    template <typename seqs_type, typename ids_type, typename quals_type, size_t max_records = 0>
        requires container_concept<typename seqs_type::value> &&
                 container_concept<typename ids_type::value> &&
                 container_concept<typename quals_type::value>
    void read(seqs_type && seqs, ids_type && ids, quals_type && quals);

    /* options */
    struct options_type
    {
        // post-processing filters that operate on buffer before assignment to out-value
        std::function<void(std::string &)>  seq_filter = [] (std::string & seq)  {};
        std::function<void(std::string &)>   id_filter = [] (std::string & id)   {};
        std::function<void(std::string &)> qual_filter = [] (std::string & qual) {};
    };
    options_type options;

protected:
    ~seq_file_in() = default;

private:
    /* file format */
    std::string file_name;
    stream_type stream;
    valid_formats format;

    /* private functions */
    void select_decompression(std::string const & compress_ext);
    template <size_t index>
    void assign_format(std::string const & ext);

    /* buffers */
};

// ------------------------------------------------------------------
// public API
// ------------------------------------------------------------------

seq_file_in::seq_file_in(std::string const & _file_name) :
        file_name(_file_name)
{
    // open stream
    stream.open(_file_name, std::ios::binary);

    // initialize format handler
    std::string ext{get_file_extension(file_name)};
    select_format<0>(format, ext);
}

template <typename seq_type, typename id_type, typename qual_type>
    requires container_concept<seq_type> &&
             container_concept<id_type> &&
             container_concept<qual_type>
inline void seq_file_in::read(seq_type && seq, id_type && id = std::string{}, qual_type && qual = std::string{})
{
    assert(!format.valueless_by_exception);
    std::visit([&] (seq_file_in_format_concept & f) { f->read(seq, id, qual, stream, options); }, format);
}

template <typename seqs_type, typename ids_type, typename quals_type, size_t num_records = 0>
    requires container_concept<typename seqs_type::value> &&
             container_concept<typename ids_type::value> &&
             container_concept<typename quals_type::value>
inline void seq_file_in::read(seqs_type  && seqs,
                              ids_type   && ids   = std::vector<std::string>{},
                              quals_type && quals = std::vector<std::string>{})
{
    assert(!format.valueless_by_exception);
    std::visit([&] (seq_file_in_format_concept & f)
    {
        f->read<num_records>(seqs, ids, quals, stream, options);
    }, format);
}

// ------------------------------------------------------------------
// private functions
// ------------------------------------------------------------------

inline void
seq_file_in::select_decompression(std::string const & compress_ext)
{
    for (auto const & pair : valid_compression_formats)
    {
        if (compress_ext == std::get<0>(pair))
        {
            std::visit([&stream] (auto & compressor) { stream.push(compressor); }, std::get<1>(pair));
            break;
        }
    }
}

template <size_t index>
inline void seq_file_in::select_format(std::string const & ext)
{
    if (index == variant_size_v<valid_formats>)
        throw std::runtime_error("No valid format found for this extension");
    else if (variant_alternative_t<index, valid_formats>::file_extensions().contains(ext))
        format = variant_alternative_t<index, valid_formats>{};
    else
        select_format<index+1>(format, ext);
}

} // namespace seqan
