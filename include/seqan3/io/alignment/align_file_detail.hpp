#pragma once

#include <string>
#include <variant>
#include <vector>

namespace seqan3::detail
{

// ==================================================================
// align_file
// ==================================================================

template <typename align_file_traits>
class align_file
{
public:
    /* types */
    using align_file_traits::stream_type;               // e.g. std::ostream concept
    using align_file_traits::valid_compression_formats; // = std::vector<std::pair<std::string, t>>
    using align_file_traits::valid_formats;             // = std::variant<gzip_compressor, align_file_format_bam>

    struct options_type
    {
        std::string program_name;    // SAM/BAM/CRAM @PG: ID + PN | BLAST header first line
        std::string program_version; // SAM/BAM/CRAM @PG: VN | BLAST header first line
        std::string command_line;    // SAM/BAM/CRAM @PG: CL | BLAST header first line
        std::vector<std::string> additional_comment_lines; // SAM/BAM/CRAM @CO lines | BLAST additional header lines
        bool sorted_by_query; // SAM/BAM/CRAM "@GO query" | enforced for blast m9 and m0, optional for m8
        //TODO more stuff
    };

    struct align_record
    {
        using field_variant = std::variant<int8_t, uint8_t, int16_t, uint16_t, int32_t, uint32_t, int64_t, uint64_t,
                                           float, double, std::string, std::string_view,
                                           ranges::any_random_access_view<char&>, ranges::any_random_access_view<char>,
                                           std::vector<int8_t>, std::vector<uint8_t>, std::vector<int16_t>,
                                           std::vector<uint16_t>, std::vector<int32_t>, std::vector<uint32_t>,
                                           std::vector<int64_t>, std::vector<uint64_t>, std::vector<float>,
                                           std::vector<double>>;
        size_t qry_no;
        size_t qry_begin;
        size_t qry_end;

        size_t sbj_no;
        size_t sbj_begin;
        size_t sbj_end;

        align_file_traits::query_gaps_type qry_gaps;
        align_file_traits::subject_gaps_type sbj_gaps;

        std::vector<std::pair<align_record_fields, field_variant> additional_fields;
        std::vector<std::tuple<char[2] /*sam_bam_tag_id*/,
                               std::string /*columnlabel*/,
                               field_variant> custom_fields;
    };

    /* constructors*/

    // constructor with arg
    align_file(std::string const & _file_name);
    // TODO add contructor with stream object that needs manual
    // type selection

    // copy construction and assignment are deleted
    // implicitly because we don't want multiple access to file
    align_file(align_file const &) = delete;
    align_file & operator=(align_file const &) = delete;

    // move construction and assignment are defaulted
    align_file(align_file &&) = default;
    align_file & operator=(align_file &&) = default;

    /* member variables */
    options_type options;

    /* member functions */

    // set store
    void set_store(std::tuple<query_seqs_type const &,
                              query_ids_type const &,
                              subject_seqs_type const &,
                              subject_ids_type const &> const &);
protected:
    ~align_file() = default;

    /* member types */

    /* member variables */
    std::string file_name;
    stream_type stream;
    valid_formats format;
    store_type store;

    /* member functions */
    void select_decompression(std::string const & compress_ext);
    template <size_t index>
    void assign_format(std::string const & ext);
};

// ------------------------------------------------------------------
// public API
// ------------------------------------------------------------------

align_file::align_file(std::string const & _file_name) :
        file_name(_file_name)
{
    // open stream
    stream.open(_file_name, std::ios::binary);

    // initialize format handler
    std::string ext{get_file_extension(file_name)};
    select_format<0>(format, ext);
}

//TODO set_store

// ------------------------------------------------------------------
// private functions
// ------------------------------------------------------------------

inline void
align_file::select_decompression(std::string const & compress_ext)
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
inline void align_file::select_format(std::string const & ext)
{
    if (index == variant_size_v<valid_formats>)
        throw std::runtime_error("No valid format found for this extension");
    else if (variant_alternative_t<index, valid_formats>::file_extensions().contains(ext))
        format = variant_alternative_t<index, valid_formats>{stream, seq_filter, id_filter, qual_filter};
    else
        select_format<index+1>(format, ext);
}

} // namespace seqan::detail