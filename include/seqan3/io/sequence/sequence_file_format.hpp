#pragma once

#include <utility>
#include <variant>

namespace seqan3
{

template <typename t>
concept bool sequence_file_format_concept = requires (t v)
{
    t::file_extensions;

    // TODO stream_type
    // TODO constructor

    //TODO how do you specifiy a function template?
    // v.read(sequence_type & seq, id_type & id, qual_type & qual);

};

namespace detail
{

template <typename variant_type, std::size_t ...I>
constexpr bool meets_sequence_file_format_concept(std::index_sequence<I...>)
{
    return (sequence_file_format_concept<std::variant_alternative_t<I, variant_type>> && ...);
}

} // namespace detail
} // namespace seqan3
