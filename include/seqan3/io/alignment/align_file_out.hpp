#pragma once

#include <string>
#include <variant>
#include <vector>

#if 0
//TODO(rrahn): this is a prototype and needs more refinement, disabling for now
//!\cond

namespace seqan3
{
// ==================================================================
// align_file_out_traits
// ==================================================================

struct align_file_out_default_dna_traits
{
    using stream_type = std::ifstream;
    using valid_formats = std::variant<align_file_out_format_sam,
                                       align_file_out_format_bam,
                                       align_file_out_format_blast_tabular
                                       align_file_out_format_blast_tabular_comments,
                                       align_file_out_format_blast_report>;
    static constexpr std::vector<std::pair<std::string, void>> valid_compression_formats{};

    using query_seqs_type   = std::vector<dna_vector>;
    using query_ids_type    = std::vector<std::string>;
    using subject_seqs_type = std::vector<dna_vector>;
    using subject_ids_type  = std::vector<std::string>;
};

struct align_file_out_default_aa_traits : public align_file_out_default_dna_traits
{
    using query_seqs_type   = std::vector<aa27_vector>;
    using subject_seqs_type = std::vector<aa27_vector>;
};

static_assert(align_file_traits_concept<align_file_out_default_dna_traits>);
static_assert(align_file_traits_concept<align_file_out_default_aa_traits>);

// ==================================================================
// align_file_out
// ==================================================================

template <typename align_file_out_traits = align_file_out_default_dna_traits>
    requires align_file_traits_concept<align_file_out_traits>
class align_file_out : public detail::align_file<align_file_out_traits>
{
public:
    /* types */

    /* constructors*/

    // constructor with arg
    align_file_out(std::string const & _file_name) :
        align_file{_file_name}
    {}

    // copy construction and assignment are deleted
    // implicitly because we don't want multiple access to file
    align_file_out(align_file_out const &) = delete;
    align_file_out & operator=(align_file_out const &) = delete;

    // move construction and assignment are defaulted
    align_file_out(align_file_out &&) = default;
    align_file_out & operator=(align_file_out &&) = default;

    /* member variables */

    /* member functions */

    // high level API
    template <typename record_type>
        requires align_record_concept<record_type> ||
        (requires container_concept<record_type> &&
        requires align_record_concept<typename record_type::value_type>)
    void write_record(record_type && r);

    // low-level API
    template <typename ...types>
    void write_raw(...);

protected:
    ~align_file_out() = default;

    /* member types */
    struct store_type
    {
        // set by set_store
        align_file_out_traits::query_seqs_type   * qry_seqs{nullptr};
        align_file_out_traits::query_ids_type    * qry_ids{nullptr};
        align_file_out_traits::subject_seqs_type * sbj_seqs{nullptr};
        align_file_out_traits::subject_ids_type  * sbj_ids{nullptr};
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
inline void align_file_out::write_record(record_type && r)
{
    assert(!format.valueless_by_exception);
    std::visit([&] (align_file_out_format_concept & f) { f->write_record(r, stream, options, store); }, format);
}


template <typename ...arg_types>
inline void align_file_out::write_raw(arg_types && ... args)
{
    assert(!format.valueless_by_exception);
    std::visit([&] (align_file_out_format_concept & f)
    {
        f->write_raw(stream, options, store, std::forwarard<decltype(args)>(args)...);
    }, format);
}

//TODO set_store


} // namespace seqan

//!\endcond
#endif
