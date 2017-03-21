#pragma once

namespace seqan3
{

template <typename t>
concept bool seq_file_in_format_concept = requires (t v)
{
    t::file_extensions;

    // TODO stream_type
    // TODO constructor

    //TODO how do you specifiy a function template?
    v.read(seq_type & seq, id_type & id, qual_type & qual);

};

} // namespace seqan3