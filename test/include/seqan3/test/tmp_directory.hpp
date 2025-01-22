// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

/*!\file
 * \author Simon Gene Gottlieb <simon.gottlieb AT fu-berlin.de>
 * \brief Internal test infrastructure.
 *
 * Define some helper classes and functions for the tests, that would be misplaced in the seqan3/include directory.
 */

#if defined(__APPLE__)
#    include <unistd.h>
#elif defined(_WIN32)
#    include <cstring>
#    include <io.h>
#else // other unix systems
#    include <cstdlib>
#endif

#include <cassert>
#include <filesystem>
#include <iostream>
#include <optional>
#include <utility>

#include <seqan3/core/platform.hpp>
#include <seqan3/test/sandboxed_path.hpp>

namespace seqan3::test
{
#if defined(_WIN32)
namespace
{
/*!\brief Helper function that implements mkdtemp for windows.
 *
 * Caveat:
 *
 * Race condition between processes exists.
 * There exists the possibility that between generating
 * a unique file name and creating the directory another
 * process also generates exactly the same filename.
 */
char * mkdtemp(char * template_name)
{
    if (_mktemp_s(template_name, strlen(template_name) + 1))
        return nullptr;

    if (std::filesystem::create_directories(template_name))
        return template_name;

    return nullptr;
}
} // namespace
#endif

//!\cond
/*!\brief Creates and maintains a unique temporary directory.
 *
 * Creates a temporary unique directory. It automatically removes the temporary directory and all contained files and
 * subdirectories on destruction. The class manages the life time of the associated directory. This means, when the
 * instance is destructed the associated filesystem directory and all it's contents will be deleted automatically.
 * Hence an instance of this class cannot be copied.
 * The life time of the associated directory also ends if the move-operator assigns a new associated directory.
 *
 * ### Example
 *
 * \include test/snippet/test/tmp_directory.cpp
 *
 * ### Exceptions
 *
 * Might throw a std::filesystem::filesystem_error on failure to create a temporary file directory.
 *
 * ### Thread safety
 *
 * According to "https://www.gnu.org/software/libc/manual/html_node/Temporary-Files.html" the call to
 * `mkdtemp` is thread safe, such that creating multiple parallel instances of this class will
 * not induce a data race on the creation of temporary file path.
 */
class tmp_directory
{
public:
    /* rule of six */
    /*!\name Constructors, destructor and assignment
     * \{
     */
    tmp_directory(tmp_directory const &) = delete; //!< Deleted.
    tmp_directory(tmp_directory && other) : directory_path{std::exchange(other.directory_path, std::nullopt)}
    {}

    tmp_directory & operator=(tmp_directory const &) = delete; //!< Deleted.
    tmp_directory & operator=(tmp_directory && other)
    {
        /* The current hold directory is cleaned. A simple std::swap is not
         * performed to avoid prolonging the life of the temporary directory.
         */
        clean();
        directory_path = std::exchange(other.directory_path, std::nullopt);
        return *this;
    }

    /*!\brief Constructs temp path with given file name.
     *
     * The generated file name is unique due to a call to `mkdtemp`.
     *
     * ### Exceptions
     * Might throw std::filesystem::filesystem_error.
     */
    tmp_directory()
    {
        auto tmp_base_dir = std::filesystem::temp_directory_path();
        tmp_base_dir /= std::filesystem::path{"seqan_test_XXXXXXXX"};

        auto path_str = tmp_base_dir.string(); // Copy the underlying path to get access to the raw char *.

        if (char * f = mkdtemp(path_str.data()); f == nullptr) // mkdtemp replaces XXXXXXXX in a safe and unique way.
        {
            throw std::filesystem::filesystem_error("Could not create temporary directory with mkdtemp!",
                                                    tmp_base_dir,
                                                    std::make_error_code(std::errc::bad_file_descriptor));
        }
        directory_path = path_str;
    }

    /*!\brief Destructs the temporary directory path.
     *
     * Removes the temporary directory and all its subdirectories and files contained.
     */
    ~tmp_directory()
    {
        clean();
    }
    //!\}

    /*!\brief Returns a reference to the path object.
     * \return seqan3::test::sandboxed_path containing the path of the file.
     */
    sandboxed_path path() const
    {
        assert(directory_path);

        return directory_path.value();
    }

    /*!\brief Check if the temporary directory is empty.
     * \return True if the directory is empty, false otherwise.
     */
    bool empty() const
    {
        assert(directory_path);

        return exists(directory_path.value()) && is_empty(directory_path.value());
    }

private:
    /*!\brief Warns and cleans if directory is not empty
     */
    void clean()
    {
        if (!directory_path)
        {
            return;
        }

        // Deletes all directories recursivly without following symlinks.
        // Function will throw on error.
        std::filesystem::remove_all(directory_path.value());
        directory_path = std::nullopt;
    }

private:
    std::optional<sandboxed_path> directory_path;
};
//!\endcond
} // namespace seqan3::test
