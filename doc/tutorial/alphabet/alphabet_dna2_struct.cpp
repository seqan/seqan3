//! [struct]
#include <seqan3/alphabet/concept.hpp>   // for seqan3::Alphabet concept checks
#include <seqan3/alphabet/exception.hpp> // for seqan3::invalid_char_assignment

struct dna2
{
    using rank_type = uint8_t;
    rank_type rank{};
};
//! [struct]
