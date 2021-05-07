#include <seqan3/test/snippet/create_temporary_snippet_file.hpp>
// std::filesystem::current_path() / "my.fastq" will be deleted after the execution
create_temporary_snippet_file my_fastq{"my.fastq", "\n"};

//![main]
#include <seqan3/io/sequence_file/all.hpp>

int main()
{
    seqan3::sequence_file_input fin{std::filesystem::current_path() / "my.fastq"};
    using record_type = typename decltype(fin)::record_type;

    record_type rec = std::move(*fin.begin()); // avoid copying
}
//![main]
