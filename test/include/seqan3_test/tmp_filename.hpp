// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
//
// Copyright (c) 2006-2017, Knut Reinert, FU Berlin
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
// ==========================================================================

#pragma once

/*
 * Internal test infrastructure.
 *
 * Define some helper classes and functions for the tests, that would be misplaced in the seqan3/include directory.
 * @author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#include <experimental/filesystem>

// TODO(rrahn): At support for Windows platforms, when we support it.
#if defined(__APPLE__)
#include <unistd.h>
#else  // other unix systems
#include <stdlib.h>
#endif

namespace fs = std::experimental::filesystem;

namespace seqan3
{
namespace test
{

/// \cond
/*!\brief Creates and maintains a filesystem::path to a temporary file.
 * Creates a temporary unique file directory and adds the given file name to construct a filesystem::path.
 * It automatically removes the temporary directory and all contained files and subdirectories on destruction.
 *
 * \par Example
 *
 * ```cpp
 * tmp_file_name fn{"my_file"};
 * std::cout << fn.get_path() << std::endl;
 * ```
 *
 * \par Exceptions
 *
 * Might throw a std::filesystem::filesystem_error on failure to create a temporary file directory.
 *
 * \par Thread safety
 *
 *  According to "https://www.gnu.org/software/libc/manual/html_node/Temporary-Files.html" the call to
 * \a mkdtemp is thread safe, such that creating multiple parallel instances of this class will
 * not induce a data race on the creation of temporary file path.
 */
class tmp_file_name
{
public:

    /* rule of six */
    /*!\name Constructors, destructor and assignment
     * \{
     */
    //!\brief Deleted default constructor.
    tmp_file_name() = delete;
    //!\brief Copy constructor.
    tmp_file_name(tmp_file_name const &) = default;
    //!\brief Move constructor.
    tmp_file_name(tmp_file_name &&) = default;
    //!\brief Copy assignment.
    tmp_file_name & operator=(tmp_file_name const &) = default;
    //!\brief Move assignment.
    tmp_file_name & operator=(tmp_file_name &&) = default;

    /*!\brief Constructs temp path with given file name.
     * \param f_name The name of the file.
     *
     * The generated file name is unique due to a call to \a mkdtemp.
     *
     * \par Exceptions
     * Might throw std::filesystem::filesystem_error.
     */
    explicit tmp_file_name(const char * f_name)
    {
        if (f_name == nullptr)
            throw fs::filesystem_error("Empty file name!",
                                        std::make_error_code(std::errc::invalid_argument));

        auto tmp_base_dir = fs::temp_directory_path();
        tmp_base_dir += fs::path{"seqan_test_XXXXXXXX"};
        // Sanity check, to test that the data stream of the tmp base name is in fact not const char *.
        static_assert(std::is_same_v<char*, decltype(tmp_base_dir.string().data())>);
        // We have to use mkdtemp, which is not deprecated. We place it into the dedicated tmp_dir
        // returned by temp_directory_path. Within this path we can safely create files, that would be
        // unique per test instance as the parent directory is.
        if (char * f = mkdtemp(tmp_base_dir.string().data()); f != nullptr)
        {
            file_path = f;
            file_path /= fs::path{f_name};
            return;
        }
        throw fs::filesystem_error("Could not create temporary directory with mkdtemp!",
                                   tmp_base_dir,
                                   std::make_error_code(std::errc::bad_file_descriptor));
    }

    /*!\brief Destructs the temporary file path.
     * Removes the temporary directory and all it's subdirectories and files contained.
     */
    ~tmp_file_name()
    {
        [[maybe_unused]] std::error_code ec;
        remove_all(file_path.parent_path(), ec);
    }
    //!\}

    /*!\brief Returns a reference to the path object.
     * \returns std::filesystem::path containing the path of the file.
     */
    auto & get_path() const
    {
        return file_path;
    }

private:
    //!\brief The object storing the path to the temporary file.
    fs::path file_path{};
};
/// \endcond
} // namespace test
} // namespace seqan3
