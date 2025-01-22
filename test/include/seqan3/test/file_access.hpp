// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Detects read and write access for a path.
 * \author Enrico Seiler <enrico.seiler AT fu-berlin.de>
 */

#pragma once

#include <filesystem>
#include <fstream>

#include <seqan3/core/platform.hpp>

namespace seqan3::test
{

//!\brief Checks wheter there is read access to a path.
inline bool read_access(std::filesystem::path const & file)
{
    std::fstream stream;
    stream.open(file, std::ios::in);
    return !stream.fail();
}

//!\brief Checks wheter there is write access to a path.
inline bool write_access(std::filesystem::path const & file)
{
    if (std::filesystem::is_directory(file))
    {
        std::filesystem::path test_file{file};
        test_file /= "seqan3_test_write_access";
        std::fstream stream;
        stream.open(test_file, std::ios::out);
        bool result = !stream.fail();
        if (result)
        {
            stream.close();
            std::filesystem::remove(test_file);
        }
        return result;
    }
    else
    {
        std::fstream stream;
        stream.open(file, std::ios::out);
        return !stream.fail();
    }
}

} // namespace seqan3::test
