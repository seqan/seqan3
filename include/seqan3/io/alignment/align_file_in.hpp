#pragma once

#include <string>
#include <variant>
#include <vector>

namespace seqan3
{
// ==================================================================
// align_file_in_traits
// ==================================================================

struct align_file_in_default_dna_traits
{
    using stream_type = std::ifstream;
    using valid_formats = std::variant<align_file_in_format_sam,
                                       align_file_in_format_bam,
                                       align_file_in_format_blast_tabular
                                       align_file_in_format_blast_tabular_comments>;
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

static_assert(align_file_traits_concept<align_file_in_default_dna_traits>);
static_assert(align_file_traits_concept<align_file_in_default_aa_traits>);

// ==================================================================
// align_file_in
// ==================================================================

template <typename align_file_in_traits = align_file_in_default_dna_traits>
    requires align_file_traits_concept<align_file_in_traits>
class align_file_in : public detail::align_file<align_file_in_traits>
{
public:
    /* types */

    /* constructors*/

    // constructor with arg
    align_file_in(std::string const & _file_name) :
        align_file{_file_name}
    {}

    // copy construction and assignment are deleted
    // implicitly because we don't want multiple access to file
    align_file_in(align_file_in const &) = delete;
    align_file_in & operator=(align_file_in const &) = delete;

    // move construction and assignment are defaulted
    align_file_in(align_file_in &&) = default;
    align_file_in & operator=(align_file_in &&) = default;

    /* member variables */

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

    // set store (two parameter version in addition to 4 parameter)
    void set_store(std::tuple<subject_seqs_type const &,
                              subject_ids_type const &> const &);

protected:
    ~align_file_in() = default;

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
    store_type store;

    /* member functions */
};

// ------------------------------------------------------------------
// public API
// ------------------------------------------------------------------

template <typename record_type>
    requires align_record_concept<record_type> ||
    (requires container_concept<record_type> && requires align_record_concept<typename record_type::value_type>)
inline void align_file_in::read_record(record_type && r)
{
    assert(!format.valueless_by_exception);
    std::visit([&] (align_file_in_format_concept & f) { f->read_record(r, stream, options, store); }, format);
}


template <typename ...arg_types>
inline void align_file_in::read_raw(arg_types && ... args)
{
    assert(!format.valueless_by_exception);
    std::visit([&] (align_file_in_format_concept & f)
    { 
        f->read_raw(stream, options, store, std::forwarard<decltype(args)>(args)...);
    }, format);
}

//TODO set_store


} // namespace seqan