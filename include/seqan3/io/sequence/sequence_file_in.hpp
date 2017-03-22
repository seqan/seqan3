#pragma once

#include <fstream>
#include <string>
#include <variant>
#include <vector>
#include <functional>
#include "for_test.hpp"
#include "sequence_file_format.hpp"
#include "sequence_file_format_fasta.hpp"

namespace seqan3
{

// ==================================================================
// sequence_file_in_traits
// ==================================================================
template <typename t>
concept bool sequence_file_in_traits_concept = requires (t v)
{
    typename t::stream_type;
    typename t::valid_formats;
    // //TODO valid_formats shall be variant and all of variants types
    // // must meet format_concept

    t::valid_compression_formats;
};

struct sequence_file_in_default_traits
{
    using stream_type = std::ifstream;
    using valid_formats = std::variant<sequence_file_format_fasta>;
                                       // sequence_file_format_fastq,
                                       // sequence_file_format_embl,
                                       // sequence_file_format_genbank,
                                       // sequence_file_format_raw>;

    using valid_compressions = std::variant<char/*specify compression formats*/>;

    static inline std::vector<std::pair<std::string, valid_compressions>> valid_compression_formats
    {
    //     {".gz", variant_alternative_t<0, valid_compressions>{}},
    //     {".bz2", variant_alternative_t<1, valid_compressions>{}}
    };
};

static_assert(sequence_file_in_traits_concept<sequence_file_in_default_traits>);

// ==================================================================
// sequence_file_in
// ==================================================================

template <typename sequence_file_in_traits = sequence_file_in_default_traits>
    requires sequence_file_in_traits_concept<sequence_file_in_traits>
class sequence_file_in
{
public:
    using stream_type = typename sequence_file_in_traits::stream_type;               // e.g. std::ostream concept
    using valid_formats = typename sequence_file_in_traits::valid_formats;             // = std::variant<gzip_compressor, sequence_file_format_bam>

    // constructor with arg
    sequence_file_in(std::string const & _file_name)
    {
        file_name = _file_name;
        // open stream
        stream.open(_file_name, std::ios::binary);

        // initialize format handler
        std::string ext{"fasta"};
        // select_format<0>(ext);
        format = sequence_file_format_fasta {};
    }
    // copy construction and assignment are deleted
    // implicitly because we don't want multiple access to file
    sequence_file_in(sequence_file_in const &) = delete;
    sequence_file_in & operator=(sequence_file_in const &) = delete;

    // move construction and assignment are defaulted
    sequence_file_in(sequence_file_in &&) = default;
    sequence_file_in & operator=(sequence_file_in &&) = default;

    //TODO make the requirements stricter
    template <typename sequence_type, typename id_type, typename qual_type>
    requires container_concept<sequence_type> &&
             container_concept<id_type> &&
             container_concept<qual_type>
    inline void read(sequence_type && seq,
                     id_type && id = std::string{},
                     qual_type && qual = std::string{});


    template <typename seqs_type,
              typename ids_type,
              typename quals_type,
              size_t max_records = 0>
    requires container_concept<typename seqs_type::value> &&
             container_concept<typename ids_type::value> &&
             container_concept<typename quals_type::value>
    inline void read(seqs_type && seqs,
                     ids_type && ids = std::vector<std::string>{},
                     quals_type && quals = std::vector<std::string>{});

    /* options */
    struct options_type
    {
        // post-processing filters that operate on buffer before assignment to out-value
        std::function<void(std::string &)>  sequence_filter = [] (std::string & seq)  {};
        std::function<void(std::string &)>   id_filter = [] (std::string & id)   {};
        std::function<void(std::string &)> qual_filter = [] (std::string & qual) {};
    };
    options_type options;

    ~sequence_file_in() = default;

protected:
    /* file format */
    std::string file_name;
    stream_type stream;
    valid_formats format;

    /* private functions */
    void select_decompression(std::string const & compress_ext);
    template <size_t index>
    void select_format(std::string const & ext);
    /* buffers */
};

// ------------------------------------------------------------------
// public API
// ------------------------------------------------------------------
template <typename sequence_file_in_traits>
template <typename sequence_type, typename id_type, typename qual_type>
requires container_concept<sequence_type> &&
         container_concept<id_type> &&
         container_concept<qual_type>
inline void sequence_file_in<sequence_file_in_traits>::read(sequence_type && seq,
                                                            id_type && id,
                                                            qual_type && qual)
{
    if (!format.valueless_by_exception())
    {
        throw "format not set!";
    }
    std::visit([&] (sequence_file_format_concept & f) {f.read(seq, id, qual, stream, options);}, format);
}

template <typename sequence_file_in_traits>
template <typename seqs_type,
          typename ids_type,
          typename quals_type,
          size_t num_records>
requires container_concept<typename seqs_type::value> &&
         container_concept<typename ids_type::value> &&
         container_concept<typename quals_type::value>
inline void sequence_file_in<sequence_file_in_traits>::read(seqs_type && seqs,
                                                            ids_type && ids,
                                                            quals_type && quals)
{
    if (!format.valueless_by_exception())
    {
        throw "format not set!";
    }
    std::visit([&] (sequence_file_format_concept & f)
    {
        f.read<num_records>(seqs, ids, quals, stream, options);
    }, format);
}

// ------------------------------------------------------------------
// private functions
// ------------------------------------------------------------------

template <typename sequence_file_in_traits>
void sequence_file_in<sequence_file_in_traits>::select_decompression(std::string const & compress_ext)
{
    for (auto const & pair : sequence_file_in_traits::valid_compression_formats)
    {
        if (compress_ext == std::get<0>(pair))
        {
            std::visit([&] (auto & compressor) { stream.push(compressor); }, std::get<1>(pair));
            break;
        }
    }
}

template <typename sequence_file_in_traits>
template <size_t index>
void sequence_file_in<sequence_file_in_traits>::select_format(std::string const & ext)
{
    if (index == std::variant_size_v<valid_formats>)
        throw std::runtime_error("No valid format found for this extension");
    else if (std::variant_alternative_t<index, valid_formats>::file_extensions().contains(ext))
        format = std::variant_alternative_t<index, valid_formats>{};
    else
        select_format<index+1>(format, ext);
}
} // namespace seqan