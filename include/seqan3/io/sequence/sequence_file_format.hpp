#pragma once

namespace seqan3
{

template <typename t>
concept bool sequence_file_format_concept = requires (t v)
{
    t::file_extensions;

    // TODO stream_type
    // TODO constructor

    //TODO how do you specifiy a function template?
    v.read(sequence_type & seq, id_type & id, qual_type & qual);

};

namespace detail
{

template <size_t index, typename variant_type>
constexpr bool meets_concept_sequence_file_format()
{
    if constexpr (index == variant_size_v<variant_type>)
        return true;
    else if constexpr (!sequence_file_format_concept<variant_alternative_t<index, variant_type>>)
        return false;
    else
        return meets_concept_sequence_file_format<index+1, variant_type>();
}

} // namespace detail
} // namespace seqan3
