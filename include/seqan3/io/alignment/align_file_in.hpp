#pragma once

#include <string>
#include <variant>
#include <vector>

namespace seqan3
{

// ==================================================================
// align_file_in_traits
// ==================================================================

template <typename t>
concept bool align_file_in_traits_concept = requires (t v)
{
    t::stream_type;
    t::valid_formats;
    //TODO valid_formats shall be variant and all of variants types
    // must meet format_concept

    t::valid_compression_formats;

    // store types
    t::query_seqs_type;
    t::query_ids_type;
    t::subject_seqs_type;
    t::subject_ids_type;

    // align_record
    t::qry_gaps_type;
    t::sbj_gaps_type;
};

struct align_file_in_default_dna_traits
{
    using stream_type = std::ofstream;
    using valid_formats = std::variant<align_file_in_format_sam,
                                       align_file_in_format_bam,
                                       align_file_in_format_blast_tabular
                                       align_file_in_format_blast_tabular_comments,
                                       align_file_in_format_blast_report>;
    static constexpr std::vector<std::pair<std::string, void>> valid_compression_formats{};

    using query_seqs_type   = std::vector<dna_vector>;
    using query_ids_type    = std::vector<std::string>;
    using subject_seqs_type = std::vector<dna_vector>;
    using subject_ids_type  = std::vector<std::string>;
};

struct align_file_in_default_aa_traits : public align_file_in_default_dna_traits
{
    using query_seqs_type   = std::vector<aa27_vector>;
    using subject_seqs_type = std::vector<aa27_vector>;
};

// ==================================================================
// align_file_in
// ==================================================================

template <typename align_file_in_traits = align_file_in_default_dna_traits>
    requires align_file_in_traits_concept<align_file_in_traits>
class align_file_in
{
public:
    /* types */
    using align_file_in_traits::stream_type;               // e.g. std::ostream concept
    using align_file_in_traits::valid_compression_formats; // = std::vector<std::pair<std::string, t>>
    using align_file_in_traits::valid_formats;             // = std::variant<gzip_compressor, align_file_in_format_bam>

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
                                           float, double, std::string, std::string_view, std::vector<int8_t>,
                                           std::vector<uint8_t>, std::vector<int16_t>, std::vector<uint16_t>,
                                           std::vector<int32_t>, std::vector<uint32_t>, std::vector<int64_t>,
                                           std::vector<uint64_t>, std::vector<float>, std::vector<double>>;
        size_t qry_no;
        size_t qry_begin;
        size_t qry_end;

        size_t sbj_no;
        size_t sbj_begin;
        size_t sbj_end;

        align_file_in_traits::gaps_type qry_gaps;
        align_file_in_traits::gaps_type sbj_gaps;

        std::vector<std::pair<our_enum, field_variant> additional_fields;
        std::vector<std::tuple<std::string /*sam_bam_tag_id*/,
                               std::string /*columnlabel*/,
                               field_variant> custom_fields;
    };

    /* constructors*/

    // constructor with arg
    align_file_in(std::string const & _file_name);

    // copy construction and assignment are deleted
    // implicitly because we don't want multiple access to file
    align_file_in(align_file_in const &) = delete;
    align_file_in & operator=(align_file_in const &) = delete;

    // move construction and assignment are defaulted
    align_file_in(align_file_in &&) = default;
    align_file_in & operator=(align_file_in &&) = default;

    /* member variables */
    options_type options;

    /* member functions */

    // high level API
    template <typename record_type>
        requires align_record_concept<record_type> ||
        (requires container_concept<record_type> &&
        requires align_record_concept<typename record_type::value_type>)
    void read_record(record_type && r);

    // low-level API
    template <typename ...types>
    void read_raw(...);

    // set store
    void set_store(std::tuple<query_seqs_type const &,
                              query_ids_type const &,
                              subject_seqs_type const &,
                              subject_ids_type const &> const &);
    void set_store(std::tuple<subject_seqs_type const &,
                              subject_ids_type const &> const &);
protected:
    ~align_file_in() = default;

private:
    /* member types */
    struct store_type
    {
        // set by set_store
        align_file_in_traits::query_seqs_type   * qry_seqs{nullptr};
        align_file_in_traits::query_ids_type    * qry_ids{nullptr};
        align_file_in_traits::subject_seqs_type * sbj_seqs{nullptr};
        align_file_in_traits::subject_ids_type  * sbj_ids{nullptr};

        // possibly filled from file
        align_file_in_traits::query_seqs_type   qry_seqs_from_file;
        align_file_in_traits::query_ids_type    qry_ids_from_file;
        align_file_in_traits::subject_seqs_type sbj_seqs_from_file;
        align_file_in_traits::subject_ids_type  sbj_ids_from_file;

        // required for correct association if reading from file
        align_file_in_traits::query_id_map_type   query_id_map;
        align_file_in_traits::subject_id_map_type subject_id_map;
    };

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

align_file_in::align_file_in(std::string const & _file_name) :
        file_name(_file_name)
{
    // open stream
    stream.open(_file_name, std::ios::binary);

    // initialize format handler
    std::string ext{get_file_extension(file_name)};
    select_format<0>(format, ext);
}

template <typename record_type>
    requires align_record_concept<record_type> ||
    (requires container_concept<record_type> && requires align_record_concept<typename record_type::value_type>)
inline void align_file_in::read_record(record_type && r)
{
    assert(!format.valueless_by_exception);
    std::visit([&] (align_file_in_format_concept & f) { f->read(r, stream, options, store); }, format);
}

// ------------------------------------------------------------------
// private functions
// ------------------------------------------------------------------

inline void
align_file_in::select_decompression(std::string const & compress_ext)
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
inline void align_file_in::select_format(std::string const & ext)
{
    if (index == variant_size_v<valid_formats>)
        throw std::runtime_error("No valid format found for this extension");
    else if (variant_alternative_t<index, valid_formats>::file_extensions().contains(ext))
        format = variant_alternative_t<index, valid_formats>{stream, seq_filter, id_filter, qual_filter};
    else
        select_format<index+1>(format, ext);
}

} // namespace seqan