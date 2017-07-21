// ============================================================================
//                 SeqAn - The Library for Sequence Analysis
// ============================================================================
//
// Copyright (c) 2006-2017, Knut Reinert & Freie Universitaet Berlin
// Copyright (c) 2016-2017, Knut Reinert & MPI Molekulare Genetik
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

/*!\file
 * \brief Contains test cases for sequence io (different formats).
 * \author Temesgen H. Dadi <temesgen.dadi@fu-berlin.de>
 * \ingroup test/io/sequence
 */

#include <gtest/gtest.h>
#include <fstream>
#include <experimental/filesystem>
#include <string>
#include <vector>
#include <iterator>
namespace filesystem = std::experimental::filesystem;

#include <range/v3/utility/iterator.hpp>
#include <range/v3/istream_range.hpp>

#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/io/sequence/sequence_file_in.hpp>

using namespace seqan3;

std::string inline create_temp_file(std::string content, std::string ext = ".txt")
{
    std::string filename = filesystem::temp_directory_path();
    filename += "temp_test_file";
    filename += ext;
    std::ofstream tmp_file(filename, tmp_file.binary);
    if (!tmp_file.is_open())
        throw ("[ERROR] failed to create a temporary test file!");
    tmp_file << content;
    tmp_file.close();
    return filename;
}

std::string fasta_text{">seq1\nCGATCGATAAT\n>seq2\nCCTCTCTCTCCCT\n>seq3\nCCCCCCCC\n"};
std::vector<std::string> expected_ids{"seq1", "seq2", "seq3"};
std::vector<std::string> expected_seqs{"CGATCGATAAT", "CCTCTCTCTCCCT", "CCCCCCCC"};
std::string filename = create_temp_file(fasta_text, ".fa");

TEST(fasta_format_io_test, read_single)
{
    {
        sequence_file_in fasta_file(filename);
        std::vector<std::string> ids;
        int i = 0;
        while(!fasta_file.eof())
        {
            std::string id;
            std::string seq;
            fasta_file.read(seq, id);
            EXPECT_EQ(id, expected_ids[i]);
            EXPECT_EQ(seq, expected_seqs[i]);
            i++;
        }
    }
}

TEST(fasta_format_io_test, read_batch)
{
    {
        sequence_file_in fasta_file(filename);
        std::vector<std::string> ids;
        std::vector<std::string> seqs;
        fasta_file.read(seqs, ids, 3);
        EXPECT_EQ(ids, expected_ids);
        EXPECT_EQ(seqs, expected_seqs);
    }
}
