// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2018, Knut Reinert, FU Berlin
// Copyright (c) 2016-2018, Knut Reinert & MPI Molekulare Genetik
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

#include <sstream>

#include <gtest/gtest.h>

#include <range/v3/view/zip.hpp>

#include <seqan3/core/concept/iterator.hpp>
#include <seqan3/io/sequence/sequence_file_in.hpp>
#include <seqan3/range/view/to_char.hpp>

using namespace seqan3;
using namespace seqan3::literal;

TEST(sequence_file_in_iterator, concepts)
{
    using it_t = typename sequence_file_in<>::iterator;
    using sen_t = typename sequence_file_in<>::sentinel;

    EXPECT_TRUE((seqan3::input_iterator_concept<it_t>));
    EXPECT_TRUE((seqan3::sentinel_concept<sen_t, it_t>));

    // const_iterator is not a an iterator, i.e. files are not const-iterable:
    using cit_t = typename sequence_file_in<>::const_iterator;
    EXPECT_FALSE((seqan3::input_iterator_concept<cit_t>));
}

TEST(sequence_file_in, concepts)
{
    using t = sequence_file_in<>;
    EXPECT_TRUE((seqan3::input_range_concept<t>));

    using ct = sequence_file_in<> const;
    // not const-iterable
    EXPECT_FALSE((seqan3::input_range_concept<ct>));
}

TEST(sequence_file_in, first_test)
{
    std::string input
    {
        "> TEST1\n"
        "ACGT\n"
        "> Test2\n"
        "AGGCTGN\n"
        "> Test3\n"
        "GGAGTATAATATATATATATATAT\n"
    };

    std::istringstream iss(input);

    sequence_file_in fin{std::move(iss), sequence_file_format_fasta{}};

    for (auto [ seq, id ] : fin)
    {
        std::cout << "ID:  " << id << '\n';
        std::cout << "SEQ: " << (seq | view::to_char) << '\n';

        // ensure that record elements where moved:
        auto && r = fin.back();
        EXPECT_TRUE((empty(std::get<0>(r))));
        EXPECT_TRUE((empty(std::get<1>(r))));
    }

}
