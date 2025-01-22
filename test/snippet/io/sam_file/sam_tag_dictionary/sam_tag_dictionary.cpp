// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

//! [type_overload]
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>

using namespace seqan3::literals;

template <>                           // no template parameter since the tag is known
struct seqan3::sam_tag_type<"XX"_tag> // here comes your tag
{
    using type = int32_t; // specify the type of your tag
};
//! [type_overload]

//! [tag]
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>

using namespace seqan3::literals;

// ...

uint16_t tag_id = "NM"_tag; // tag_id = 10061
//! [tag]

//! [tag_type_t]
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>

using namespace seqan3::literals;

// ...

using nm_tag_type = seqan3::sam_tag_type_t<"NM"_tag>;
//! [tag_type_t]

//! [tag_type]
#include <seqan3/io/sam_file/sam_tag_dictionary.hpp>

using namespace seqan3::literals;

// ...

using nm_tag_type2 = seqan3::sam_tag_type<"NM"_tag>::type;
//! [tag_type]
