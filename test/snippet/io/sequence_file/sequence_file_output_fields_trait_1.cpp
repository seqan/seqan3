// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/quality/phred42.hpp>
#include <seqan3/alphabet/quality/qualified.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/utility/views/elements.hpp>

int main()
{
    using namespace seqan3::literals;

    seqan3::sequence_file_output fout{std::ostringstream{},
                                      seqan3::format_fastq{},
                                      seqan3::fields<seqan3::field::id, seqan3::field::seq, seqan3::field::qual>{}};

    for (int i = 0; i < 5; i++)
    {
        std::string id{"test_id"};
        // vector of combined data structure:
        std::vector<seqan3::qualified<seqan3::dna5, seqan3::phred42>> seq_qual{{'N'_dna5, '7'_phred42},
                                                                               {'A'_dna5, '1'_phred42},
                                                                               {'C'_dna5, '3'_phred42}};

        auto view_on_seq = seqan3::views::elements<0>(seq_qual);
        auto view_on_qual = seqan3::views::elements<1>(seq_qual);

        // ...

        // Note that the order of the arguments is different from the default `seq, id, qual`,
        // because you specified that ID should be first in the fields template argument.
        fout.emplace_back(id, view_on_seq, view_on_qual);

        // or:
        fout.push_back(std::tie(id, view_on_seq, view_on_qual));
    }
}
