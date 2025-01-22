// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC-BY-4.0

#include <filesystem>

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
