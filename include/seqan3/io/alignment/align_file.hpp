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
concept bool align_file_traits_concept = requires (t v)
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

}