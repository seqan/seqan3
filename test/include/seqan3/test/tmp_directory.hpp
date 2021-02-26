// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2021-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2021-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

#pragma once

/*
 * Internal test infrastructure.
 *
 * Define some helper classes and functions for the tests, that would be misplaced in the seqan3/include directory.
 * @author Simon Gene Gottlieb <simon.gottlieb AT fu-berlin.de>
 */

#if defined(__APPLE__)
#include <unistd.h>
#elif defined(_WIN32)
#include <cstring>
#include <io.h>
#else  // other unix systems
#include <cstdlib>
#endif

#include <seqan3/std/filesystem>
#include <iostream>

#include <seqan3/core/platform.hpp>
#include <seqan3/test/sandboxed_path.hpp>

namespace seqan3
{
namespace test
{
#if defined(_WIN32)
namespace
{
/*
 * Helper function that implements mkdtemp for windows.
 * Caveat:
 *
 * Race condition between processes exists.
 * There exists the possibility that between generating
 * an unique file name and creating the directory another
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
}
#endif

/// \cond
/*!\brief Creates and maintains a unique temporary directory.
 * Creates a temporary unique directory. It automatically removes the temporary directory and all contained files and
 * subdirectories on destruction. The class manages the life time of the associated directory. This means, when the
 * instance is destructed the associated filesystem directory and all it's contents will be deleted automatically.
 * Hence an instance of this class cannot be copied.
 *
 * ###Example
 *
 * ```cpp
 * tmp_directory fn;
 * std::cout << fn.get_path() << std::endl;
 * ```
 *
 * or in a typicall unittest
 *
 * ```cpp
 * tmp_directory fn;
 *
 *  ... some method that creates temporary files ...
 *
 * EXPECT_TRUE(fn.empty()); // or call fn.clean() to delete all files inside the temporary directory
 * ```
 *
 * ###Exceptions
 *
 * Might throw a std::filesystem::filesystem_error on failure to create a temporary file directory.
 *
 * ###Thread safety
 *
 *  According to "https://www.gnu.org/software/libc/manual/html_node/Temporary-Files.html" the call to
 * \a mkdtemp is thread safe, such that creating multiple parallel instances of this class will
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
    tmp_directory(tmp_directory &&) = default; //!< Defaulted.
    tmp_directory & operator=(tmp_directory const &) = delete; //!< Deleted.
    tmp_directory & operator=(tmp_directory &&) = default; //!< Defaulted.

    /*!\brief Constructs temp path with given file name.
     * \param f_name The name of the file.
     *
     * The generated file name is unique due to a call to \a mkdtemp.
     *
     * ###Exceptions
     * Might throw std::filesystem::filesystem_error.
     */
    explicit tmp_directory()
    {
        auto tmp_base_dir = std::filesystem::temp_directory_path();
        tmp_base_dir /= std::filesystem::path{"seqan_test_XXXXXXXX"};

        std::filesystem::create_directory(tmp_base_dir);
        // We have to use mkdtemp, which is not deprecated. We place it into the dedicated tmp_dir
        // returned by temp_directory_path. Within this path we can safely create files, that would be
        // unique per test instance as the parent directory is.
        auto path_str = tmp_base_dir.string();  // Copy the underlying path to get access to the raw char *.
        if (char * f = mkdtemp(path_str.data()); f == nullptr)  // mkdtemp replaces XXXXXXXX in a safe and unique way.
        {
            throw std::filesystem::filesystem_error("Could not create temporary directory with mkdtemp!",
                                                    tmp_base_dir,
                                                    std::make_error_code(std::errc::bad_file_descriptor));
        }
        directory_path = std::make_unique<sandboxed_path>(path_str);
    }

    /*!\brief Destructs the temporary directory path.
     *
     * \details
     *
     * Removes the temporary directory and all its subdirectories and files contained.
     */
    ~tmp_directory()
    {
        if (directory_path)
        {
            if (!exists(*directory_path))
            {
                std::cerr << "temporary directory " << *directory_path << " was deleted externally. "
                             " This is discourage program behaviour\n";
                return;
            }
            if (!empty())
            {
                std::cerr << "temporary directory " << *directory_path << " has some files that should be deleted\n";
                for (auto & p: std::filesystem::recursive_directory_iterator(*directory_path))
                {
                    std::cerr << "- " << p << "\n";
                }
            }
            std::filesystem::remove_all(*directory_path);
        }
    }
    //!\}


    /*!\brief Returns a reference to the path object.
     * \return seqan3::test::sandboxed_path containing the path of the file.
     */
    sandboxed_path const& get_path() const
    {
        assert(directory_path);

        return *directory_path;
    }

    /*!\brief Check if the temporary directory is empty.
     * \return True if the directory is empty, false otherwise.
     */
    bool empty() const
    {
        assert(directory_path);

        return exists(*directory_path) && is_empty(*directory_path);
    }

    /*!\brief Removes all files from the temporary directory
     */
    void clean()
    {
        assert(directory_path);

        for (auto & p: std::filesystem::directory_iterator(*directory_path))
        {
            std::filesystem::remove_all(p);
        }
    }

private:
    std::unique_ptr<sandboxed_path> directory_path;


};
/// \endcond
} // namespace test
} // namespace seqan3
