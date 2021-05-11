#include <seqan3/test/snippet/create_temporary_snippet_file.hpp>
// std::filesystem::current_path() / "my.qq" will be deleted after the execution
create_temporary_snippet_file my_qq{"my.qq", "\n"};

//![main]
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>

int main()
{

    seqan3::debug_stream << seqan3::format_fastq::file_extensions << '\n'; // prints [fastq,fq]
    seqan3::format_fastq::file_extensions.push_back("qq");
    seqan3::debug_stream << seqan3::format_fastq::file_extensions << '\n'; // prints [fastq,fq,qq]

    seqan3::sequence_file_input fin{std::filesystem::current_path() / "my.qq"}; // detects FASTQ format
}
//![main]
