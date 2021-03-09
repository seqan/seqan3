// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::test::sandboxed_path and related free functions.
 * \author Simon Gene Gottlieb <simon.gottlieb AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/filesystem>


#include <seqan3/core/platform.hpp>

namespace seqan3::test
{

/*!\brief Utility class to stay inside a sandbox path.
 *
 * seqan3::sandboxed_path provides the same functionality as std::filesystem::path, but restricts
 * the access to a specified directory. This results in the following invariant, which is checked
 * for at appropiate places, and some caveats.
 *
 *  Invariant:
 *  - seqan3::sandboxed_path is always converted to an absolute path
 *  - seqan3::sandboxed_path always points to a file or directory inside a given sandbox directory
 *  - The sandbox directory is immutable during the life cycle of a sandboxed_path
 *
 * Caveat:
 *  - Relative paths are not possible
 *  - Some functions will leave the sandboxed environment.
 *    - calling relative_path() leaves the environment of sandboxed_path.
 */
class sandboxed_path : public std::filesystem::path
{
private:
    std::filesystem::path const sandbox_directory;

public:
    /*!\brief Construction of a sandboxed_path.
     * \param directory must be an absolute path.
     *
     * After construction, the sandboxed_path will point to `path`.
     */
    explicit sandboxed_path(std::filesystem::path directory)
        : std::filesystem::path{directory}
        , sandbox_directory{std::move(directory)}
    {
        normalise();
        check_invariant();
    }

    /*!\brief Construction from a given sandbox directory and a path within the sandbox directory.
     * \param sandbox_directory The absolute path to the sandbox directory.
     * \param path The relative or absolute path that must be inside the sandbox directory.
     */
    explicit sandboxed_path(std::filesystem::path sandbox_directory, std::filesystem::path path)
        : std::filesystem::path {std::move(path)}
        , sandbox_directory     {std::move(sandbox_directory)}
    {
        normalise();
        check_invariant();
    }

    sandboxed_path() = delete; //!< Deleted.
    sandboxed_path(sandboxed_path const&) = default; //!< Defaulted.
    sandboxed_path(sandboxed_path&&) noexcept = default; //!< Defaulted.
    ~sandboxed_path() = default; //!< Defaulted.

    /*!\brief Replaces the path with a new path.
     * \tparam path_t The type of the new_path.
     * \param  new_path The new path.
     *
     * This works the same way as std::filesystem::path::operator=
     * and additionally checks the invariant.
     *
     * ### Exceptions
     *
     * Basic exception guarantee.
     */
    template <typename path_t>
    sandboxed_path & operator=(path_t const & new_path)
    {
        std::filesystem::path::operator=(new_path);
        normalise();
        check_invariant();
        return *this;
    }

    /*!\brief Replaces the path with a new path.
     * \tparam path_t The type of the new_path.
     * \param  new_path The new path.
     *
     * This works the same as std::filesystem::path::operator=
     * and additionally checks the invariant.
     *
     * ### Exceptions
     *
     * Basic exception guarantee
     */
    template <typename path_t>
    sandboxed_path & operator=(path_t && new_path)
    {
        std::filesystem::path::operator=(std::forward<path_t>(new_path));
        normalise();
        check_invariant();
        return *this;
    }


private:
    /*!\brief Normalises the path.
     *
     * Normalisation means that the path is converted to an absolute path and
     * is lexically normalised as described in std::filesystem::path.
     *
     * \sa https://en.cppreference.com/w/cpp/filesystem/path
     */
    void normalise()
    {
#if SEQAN3_WORKAROUND_GCC_INCOMPLETE_FILESYSTEM
        // convert into an absolute path
        auto const path_string = [this]()
        {
            if (is_relative())
            {
                return (sandbox_directory / *this).string();
            } else
            {
                return string();
            }
        }();
        // Now we need to manually collapse any "." and ".."
        size_t current_pos = path_string.find('/', 0); // find first "/" character
        std::vector<std::string> path_parts;
        while (current_pos < path_string.size())
        {
            auto end_pos = path_string.find('/', current_pos + 1);
            auto word = path_string.substr(current_pos+1, end_pos-current_pos-1);
            if (word == "." || word == "")
            {
                if (!path_parts.empty() and path_parts.back() != "")
                {
                    path_parts.emplace_back("");
                }
            }
            else if (word == "..")
            {
                if (path_parts.empty())
                {
                    throw std::filesystem::filesystem_error("Path can not be normalised",
                                                            *this,
                                                            std::make_error_code(std::errc::invalid_argument));
                }
                path_parts.pop_back();
            }
            else
            {
                if (!path_parts.empty() and path_parts.back() == "")
                {
                    path_parts.back() = word;
                }
                else
                {
                    path_parts.emplace_back(word);
                }
            }
            current_pos = end_pos;
        }
        std::string normalised_path;
        for (auto const& p : path_parts)
        {
            normalised_path += '/' + p;
        }
        std::filesystem::path::operator=(normalised_path);
#else
        auto normalised_path = std::filesystem::weakly_canonical(sandbox_directory / *this);
        std::filesystem::path::operator=(normalised_path);
#endif
    }

    /*!\brief Checks the invariant.
     *
     * Checks that the invariant of the class sandboxed_path is kept.
     * See class description for invariant.
     * \throws std::filesystem::filesystem_error if invariant was violated.
     */
    void check_invariant() const
    {
        // Check that sandbox_directory is an absolute path
        if (!sandbox_directory.is_absolute())
        {
            throw std::filesystem::filesystem_error("sandbox path must be an absolute path",
                                                sandbox_directory,
                                                *this,
                                                std::make_error_code(std::errc::invalid_argument));
        }
        // Checking that *this is an absolute path
        if (!is_absolute())
        {
            throw std::filesystem::filesystem_error("sandbox path must be an absolute path",
                                                sandbox_directory,
                                                *this,
                                                std::make_error_code(std::errc::invalid_argument));
        }

#if SEQAN3_WORKAROUND_GCC_INCOMPLETE_FILESYSTEM
        auto current_dir = string();
        auto sandbox_dir = sandbox_directory.string();

        // add trailing '/'
        if (sandbox_dir.back() != '/')
        {
            sandbox_dir += '/';
        }
        if (current_dir.back() != '/')
        {
            current_dir += '/';
        }


        // Checks sandbox_dir is a prefix of current_dir
        bool starts_with_sandbox_dir = (current_dir.rfind(sandbox_dir, 0) == 0);
        if (!starts_with_sandbox_dir)
        {
            throw std::filesystem::filesystem_error("Leaving temporary directory is not allowed!",
                                                sandbox_directory,
                                                *this,
                                                std::make_error_code(std::errc::invalid_argument));
        }

#else
        auto rel_path = lexically_relative(sandbox_directory);

        // Leaving the temporary directory is not allowed.
        if (rel_path.string().find("..") == 0)
        {
            throw std::filesystem::filesystem_error("Leaving temporary directory is not allowed!",
                                                sandbox_directory,
                                                *this,
                                                std::make_error_code(std::errc::invalid_argument));
        }
#endif
    }

public:
    /*!\brief Replaces the path with a new path.
     * \tparam path_t The type of the new_path.
     * \param  new_path The new path.
     *
     * This works the same as std::filesystem::path::assign
     * and additionally checks the invariant.
     *
     * ### Exceptions
     *
     * Basic exception guarantee
     */
    template <typename path_t>
    sandboxed_path & assign(path_t const & new_path)
    {
        std::filesystem::path::assign(new_path);
        normalise();
        check_invariant();
        return *this;
    }

    /*!\brief Replaces the path with a new path.
     * \tparam input_iter_t The type of the input iterators.
     * \param  first The begin of a given range.
     * \param  last  The end of a given range.
     *
     * This works the same as std::filesystem::path::assign
     * and additionally checks the invariant.
     *
     * ### Exceptions
     *
     * Basic exception guarantee
     */
    template <typename input_iter_t>
    sandboxed_path & assign(input_iter_t first, input_iter_t last)
    {
        std::filesystem::path::assign(first, last);
        normalise();
        check_invariant();
        return *this;
    }

    /*!\brief Extends the path.
     * \tparam path_t The type of the new_path.
     * \param  new_path The new path.
     *
     * This works the same as std::filesystem::path::operator/=
     * and additionally checks the invariant.
     *
     * ### Exceptions
     *
     * Basic exception guarantee
     */
    template <typename path_t>
    sandboxed_path & operator/=(path_t const & new_path)
    {
        return append(new_path);
    }

    /*!\brief Extends the path.
     * \tparam path_t The type of the new_path.
     * \param  new_path The new path.
     *
     * This works the same as std::filesystem::path::append
     * and additionally checks the invariant.
     *
     * ### Exceptions
     *
     * Basic exception guarantee
     */
    template <typename path_t>
    sandboxed_path & append(path_t const & new_path)
    {
        std::filesystem::path::append(new_path);
        normalise();
        check_invariant();
        return *this;
    }

    /*!\brief Extends the path.
     * \tparam input_iter_t The type of the input iterators.
     * \param  first The begin of a given range.
     * \param  last  The end of a given range.
     *
     * This works the same as std::filesystem::path::append
     * and additionally checks the invariant.
     *
     * ### Exceptions
     *
     * Basic exception guarantee
     */
    template <typename input_iter_t>
    sandboxed_path & append(input_iter_t first, input_iter_t second)
    {
        std::filesystem::path::append(first, second);
        normalise();
        check_invariant();
        return *this;
    }

    /*!\brief Extends the path.
     * \tparam path_t The type of the new_path.
     * \param  new_path The new path.
     *
     * This works the same as std::filesystem::path::operator+=
     * and additionally checks the invariant.
     *
     * ### Exceptions
     *
     * Basic exception guarantee
     */
    template <typename path_t>
    sandboxed_path & operator+=(path_t const & new_path)
    {
        return concat(new_path);
    }

    /*!\brief Extends the path.
     * \tparam path_t The type of the new_path.
     * \param  new_path The new path.
     *
     * This works the same as std::filesystem::path::concat
     * and additionally checks the invariant.
     *
     * ### Exceptions
     *
     * Basic exception guarantee
     */
    template <typename path_t>
    sandboxed_path & concat(path_t const & new_path)
    {
        std::filesystem::path::concat(new_path);
        normalise();
        check_invariant();
        return *this;
    }

    /*!\brief Extends the path.
     * \tparam input_iter_t The type of the input iterators.
     * \param  first The begin of a given range.
     * \param  last  The end of a given range.
     *
     * This works the same as std::filesystem::path::concat
     * and additionally checks the invariant.
     *
     * ### Exceptions
     *
     * Basic exception guarantee
     */
    template <typename input_iter_t>
    sandboxed_path & concat(input_iter_t first, input_iter_t last)
    {
        std::filesystem::path::concat(first, last);
        normalise();
        check_invariant();
        return *this;
    }

    /*!\brief Removes the filename.
     *
     * This works the same as std::filesystem::path::remove_filename
     * and additionally checks the invariant.
     *
     * ### Exceptions
     *
     * Basic exception guarantee
     */
    sandboxed_path & remove_filename()
    {
        std::filesystem::path::remove_filename();
        normalise();
        check_invariant();
        return *this;
    }

    /*!\brief Replaces the file name.
     *
     * This works the same as std::filesystem::path::replace_filename
     * and additionally checks the invariant.
     *
     * ### Exceptions
     *
     * Basic exception guarantee
     */
    sandboxed_path & replace_filename(std::filesystem::path const & p)
    {
        std::filesystem::path::replace_filename(p);
        normalise();
        check_invariant();
        return *this;
    }

    /*!\brief Replaces the extension.
     *
     * This works the same as std::filesystem::path::replace_extension
     * and additionally checks the invariant.
     *
     * ### Exceptions
     *
     * Basic exception guarantee
     */
    sandboxed_path & replace_extension(std::filesystem::path const & replacement = std::filesystem::path{})
    {
        std::filesystem::path::replace_extension(replacement);
        normalise();
        check_invariant();
        return *this;
    }

    /*!\brief Returns sandboxed path to the parent path.
     *
     * This works the same as std::filesystem::path::parent_path
     * and additionally checks the invariant.
     *
     * ### Exceptions
     *
     * Basic exception guarantee
     */
    sandboxed_path parent_path() const
    {
        auto parent_path = std::filesystem::path::parent_path();
        return sandboxed_path{sandbox_directory, parent_path};
    }


    void swap(sandboxed_path & other)
    {
        std::filesystem::path::swap(other);
        check_invariant();
        other.check_invariant();
    }

    void clear() = delete; //!< Not implemented. Invariant requires the path to be an absolute path.
};

/*!\brief Append a path to a seqan3::test::sanboxed_path.
 * \tparam path_t The type of the path to append.
 * \param lhs The seqan3::test::sandboxed_path.
 * \param rhs The path to append.
 * \relates seqan3::test::sandboxed_path *
 *
 * This work the same as std::filesystem::operator/(std::filesystem::path&)
 * and additionally checks the invariant.
 */
template <typename Rhs>
sandboxed_path operator/(sandboxed_path lhs, Rhs const& rhs)
{
    lhs /= rhs;
    return lhs;
}

}