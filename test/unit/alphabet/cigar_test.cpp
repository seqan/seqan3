// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2018, Knut Reinert & Freie Universitaet Berlin
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
// ============================================================================

#include <iostream>

#include <gtest/gtest.h>

#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

using namespace seqan3;

// the cigar_op alphabet is tested in alphabet_test.cpp
// the general behaviour of the cigar alphabet is tested in cartesian_composition_test.hpp

template <sequence_container_concept container_type>
void initialize(container_type & container, std::string const && target)
{
    for (auto & val : target)
    {
        typename container_type::value_type cval;
        assign_char(cval, val);
        container.push_back(cval);
    }
}

TEST(cigar, to_string)
{
    cigar cig{cigar_op::M, 20u};

    auto cig_str = cig.to_string();

    EXPECT_EQ(cig_str, "20M");

    cigar cig2{cigar_op::EQ, 240u};

    auto cig_str2 = cig2.to_string();

    EXPECT_EQ(cig_str2, "240=");
}

TEST(cigar, stream_operator)
{
    // single cigar element
    cigar elem{cigar_op::M, 20u};

    std::string expected = "20M";

    std::ostringstream os;
    os << elem;

    EXPECT_EQ(os.str(), expected);

    // cigar vector
    cigar_vector cigv{{cigar_op::M, 20u}, {cigar_op::D, 2u}, {cigar_op::EQ, 240u}};

    std::string expected2 = "20M2D240=";

    std::ostringstream os2;
    os2 << cigv;

    EXPECT_EQ(os2.str(), expected2);
}

TEST(alignment_to_cigar, unequal_equal_length_error)
{
    std::vector<gapped<dna5>> ref;
    std::vector<gapped<dna5>> read;
    initialize(ref,  "AC");
    initialize(read, "ACGT-CGTAGTG");

    EXPECT_THROW(get_cigar_vector(std::make_pair(ref, read)), std::logic_error);
}


TEST(alignment_to_cigar, empty_sequences)
{
    std::vector<gapped<dna5>> ref{};
    std::vector<gapped<dna5>> read{};

    cigar_vector cigv = get_cigar_vector(std::make_pair(ref, read));

    EXPECT_EQ(cigv.size(), 0u);
}

TEST(alignment_to_cigar, default_parameters)
{
    // no gaps at the end
    std::vector<gapped<dna5>> ref;
    std::vector<gapped<dna5>> read;
    initialize(ref,  "ACGTGAT--CTG");
    initialize(read, "ACGT-CGTAGTG");

    cigar_vector cigv = get_cigar_vector(std::make_pair(ref, read));

    std::string expected = "4M1D2M2I3M";

    std::ostringstream os;
    os << cigv;

    EXPECT_EQ(os.str(), expected);

    // gaps at the end
    std::vector<gapped<dna5>> ref2;
    std::vector<gapped<dna5>> read2;
    initialize(ref2,  "---ACGTGAT--CTG--");
    initialize(read2, "-ACGT-CGTAGTG----");

    cigar_vector cigv2 = get_cigar_vector(std::make_pair(ref2, read2));

    std::string expected2 = "1P2I2M1D4M2I1M2D2P";

    std::ostringstream os2;
    os2 << cigv2;

    EXPECT_EQ(os2.str(), expected2);
}

TEST(alignment_to_cigar, cigar_with_soft_clipping)
{
    // no gaps at the end
    std::vector<gapped<dna5>> ref;
    std::vector<gapped<dna5>> read;
    initialize(ref,  "ACGTGAT--CTG");
    initialize(read, "ACGT-CGTAGTG");

    cigar_vector cigv = get_cigar_vector(std::make_pair(ref, read), 5, 60);

    std::string expected = "5S4M1D2M2I3M60S";

    std::ostringstream os;
    os << cigv;

    EXPECT_EQ(os.str(), expected);

    // gaps at the end
    std::vector<gapped<dna5>> ref2;
    std::vector<gapped<dna5>> read2;
    initialize(ref2,  "---ACGTGAT--CTG--");
    initialize(read2, "-ACGT-CGTAGTG----");

    cigar_vector cigv2 = get_cigar_vector(std::make_pair(ref2, read2), 3, 5);

    std::string expected2 = "3S1P2I2M1D4M2I1M2D2P5S";

    std::ostringstream os2;
    os2 << cigv2;

    EXPECT_EQ(os2.str(), expected2);
}

TEST(alignment_to_cigar, extended_cigar)
{
    // no gaps at the end
    std::vector<gapped<dna5>> ref;
    std::vector<gapped<dna5>> read;
    initialize(ref,  "ACGTGAT--CAG");
    initialize(read, "ACGT-CGTACTG");

    cigar_vector cigv1 = get_cigar_vector(std::make_pair(ref, read), 0, 0, true);
    cigar_vector cigv2 = get_cigar_vector(std::make_pair(ref, read), 5, 60, true);

    std::string expected1 =    "4=1D2X2I1=1X1=";
    std::string expected2 = "5S4=1D2X2I1=1X1=60S";

    std::ostringstream os1;
    std::ostringstream os2;

    os1 << cigv1;
    os2 << cigv2;

    EXPECT_EQ(os1.str(), expected1);
    EXPECT_EQ(os2.str(), expected2);
}
