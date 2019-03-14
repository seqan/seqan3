#include <seqan3/std/filesystem>

namespace seqan3
{

//!\cond
class cleanup
{
public:
    cleanup() = delete;
    cleanup(cleanup const &) = delete;
    cleanup & operator=(cleanup const &) = delete;
    cleanup(cleanup &&) = default;
    cleanup & operator=(cleanup &&) = default;

    cleanup(char const * const str) : file(str) {};

    ~cleanup()
    {
        std::filesystem::remove(file);
    }

private:
    std::string file;
};
//!\endcond

} // namespace seqan3
