//! [type_overload]
#include <seqan3/io/alignment_file/sam_tag_dictionary.hpp>
using namespace seqan3;

template <>                           // no template parameter since the tag is known
struct seqan3::sam_tag_type<"XX"_tag> // here comes your tag
{
    using type = int32_t;             // specify the type of your tag
};
//! [type_overload]

//! [tag]
#include <seqan3/io/alignment_file/sam_tag_dictionary.hpp>

using namespace seqan3;

// ...

uint16_t tag_id = "NM"_tag; // tag_id = 10061
//! [tag]


//! [tag_type_t]
#include <seqan3/io/alignment_file/sam_tag_dictionary.hpp>

using namespace seqan3;

// ...

using nm_tag_type = sam_tag_type_t<"NM"_tag>;
//! [tag_type_t]


//! [tag_type]
using nm_tag_type2 = sam_tag_type<"NM"_tag>::type;
//! [tag_type]
