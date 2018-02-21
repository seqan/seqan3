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

// TODO(rrahn): At support for Windows platforms, when we support it.
#if defined(__APPLE__)
#include <unistd.h>
#else  // other unix systems
#include <stdlib.h>
#endif

#if __has_include(<filesystem>)
#include <filesystem>
namespace fs = std::filesystem;
#else
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#endif // __has_include(experimental/filesystem)

namespace seqan3
{
namespace test
{

/// \cond
/*!\brief Creates and maintains a filesystem::path to a temporary file.
 * Creates a temporary unique file directory and adds the given file name to construct a filesystem::path.
 * It automatically removes the temporary directory and all contained files and subdirectories on destruction.
 * The class manages the life time of the associated directory. This means, when the instance is destructed
 * the associated filesystem directory and all it's contents will be deleted automatically.
 * Hence an instance of this class cannot be copied.
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
    tmp_file_name(tmp_file_name const &) = delete;  // NOTE: We could store the path in a shared_ptr and by thus reactivating the .
    //!\brief Move constructor.
    tmp_file_name(tmp_file_name &&) = default;
    //!\brief Copy assignment.
    tmp_file_name & operator=(tmp_file_name const &) = delete;
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
        tmp_base_dir /= fs::path{"seqan_test_XXXXXXXX"};
        // We have to use mkdtemp, which is not deprecated. We place it into the dedicated tmp_dir
        // returned by temp_directory_path. Within this path we can safely create files, that would be
        // unique per test instance as the parent directory is.
        auto path_str = tmp_base_dir.string();  // Copy the underlying path to get access to the raw char *.
        if (char * f = mkdtemp(path_str.data()); f != nullptr)  // mkdtemp replaces XXXXXXXX in a safe and unique way.
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
        fs::remove_all(file_path.parent_path(), ec);
    }
    //!\}

    /*!\brief Returns a const reference to the path object.
     * \returns std::filesystem::path containing the path of the file.
     */
    fs::path const & get_path() const
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
