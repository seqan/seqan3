#include <cstdlib>                              // std::rand
#include <future>                               // std::async
#include <string>                               // std::string

#include <seqan3/core/debug_stream.hpp>         // seqan3::debug_stream
#include <seqan3/io/sequence_file/input.hpp>    // seqan3::sequence_file_input
#include <seqan3/range/views/async_input_buffer.hpp>   // seqan3::views::async_input_buffer

std::string fasta_file =
R"(> seq1
ACGACTACGACGATCATCGATCGATCGATCGATCGATCGATCGATCGTACTACGATCGATCG
> seq2
ACGACTACGACGATCATCGATCGATCGATCGATCGATCGATCGATCGTACTACGATCGATCG
> seq3
ACGACTACGACGATCATCGATCGATCGATCGATCGATCGATCGATCGTACTACGATCGATCG
> seq4
ACGACTACGACGATCATCGATCGATCGATCGATCGATCGATCGATCGTACTACGATCGATCG
> seq5
ACGACTACGACGATCATCGATCGATCGATCGATCGATCGATCGATCGTACTACGATCGATCG
> seq6
ACGACTACGACGATCATCGATCGATCGATCGATCGATCGATCGATCGTACTACGATCGATCG
> seq7
ACGACTACGACGATCATCGATCGATCGATCGATCGATCGATCGATCGTACTACGATCGATCG
> seq8
ACGACTACGACGATCATCGATCGATCGATCGATCGATCGATCGATCGTACTACGATCGATCG
> seq9
ACGACTACGACGATCATCGATCGATCGATCGATCGATCGATCGATCGTACTACGATCGATCG
> seq10
ACGACTACGACGATCATCGATCGATCGATCGATCGATCGATCGATCGTACTACGATCGATCG
> seq11
ACGACTACGACGATCATCGATCGATCGATCGATCGATCGATCGATCGTACTACGATCGATCG
> seq12
ACGACTACGACGATCATCGATCGATCGATCGATCGATCGATCGATCGTACTACGATCGATCG
)";

int main()
{
    // initialise random number generator, only needed for demonstration purposes
    std::srand(std::time(nullptr));

    // create an input file from the string above
    seqan3::sequence_file_input fin{std::istringstream{fasta_file}, seqan3::format_fasta{}};

    // create the async buffer around the input file
    // spawns a background thread that tries to keep four records in the buffer
    auto v = fin | seqan3::views::async_input_buffer(4);

    // create a lambda function that iterates over the async buffer when called
    // (the buffer gets dynamically refilled as soon as possible)
    auto worker = [&v] ()
    {
        for (auto & record : v)
        {
            // pretend we are doing some work
            std::this_thread::sleep_for(std::chrono::milliseconds(std::rand() % 1000));
            // print current thread and sequence ID
            seqan3::debug_stream << "Thread: " << std::this_thread::get_id()             << '\t'
                                 << "Seq:    " << seqan3::get<seqan3::field::id>(record) << '\n';

        }
    };

    // launch two threads and pass the lambda function to both
    auto f0 = std::async(std::launch::async, worker);
    auto f1 = std::async(std::launch::async, worker);
}
