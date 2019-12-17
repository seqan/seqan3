#include <sstream>
#include <string>

#include <seqan3/io/alignment_file/output.hpp>
#include <seqan3/core/type_list/type_list.hpp>

int main()
{
    seqan3::alignment_file_output fout{std::ostringstream{}, seqan3::format_sam{}};

    seqan3::record<seqan3::type_list<uint32_t, std::string>, seqan3::fields<seqan3::field::mapq, seqan3::field::id>> r;

    // ...

    fout.push_back(r);
}
