// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Contains an hacky overload for the googletest PrintTo function.
 * \author Svenja Mehringer <svenja.mehringer AT fu-berlin.de>
 */

#pragma once

#include<seqan3/io/stream/debug_stream.hpp>

namespace seqan3
{

//!\cond DEV
//!\brief Overload for the googletest PrintTo function that always delegates to our debug_stream.
template <typename t>
    requires requires (t const & v) { {debug_stream << v}; }
void PrintTo (t const & v, std::ostream * out)
{
    debug_stream_type my_stream{*out};
    my_stream << v;
}

template <typename t>
    requires requires (t && v) { {debug_stream << v}; }
void PrintTo (t && v, std::ostream * out)
{
    debug_stream_type my_stream{*out};
    my_stream << v;
}
//!\endcond

} // namespace seqan3
