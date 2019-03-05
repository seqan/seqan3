// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2019, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2019, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE
// -----------------------------------------------------------------------------------------------------

//! [hello]
#include <seqan3/io/stream/debug_stream.hpp> // for debug_stream

using namespace seqan3;

int main()
{
    debug_stream << "Hello World!\n";
}
//! [hello]
