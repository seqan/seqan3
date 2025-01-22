// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

/*!\file
 * \brief Provides seqan3::test::sandboxed_path and related free functions.
 * \author Simon Gene Gottlieb <simon.gottlieb AT fu-berlin.de>
 */

#pragma once

#include <algorithm>
#include <filesystem>
#include <numeric>

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
    /*!\name Constructors, destructor and assignment
     * \{
     */

    /*!\brief Construction of a sandboxed_path.
     * \param directory must be an absolute path.
     *
     * After construction, the sandboxed_path will point to `path`.
     */
    explicit sandboxed_path(std::filesystem::path directory) :
        std::filesystem::path{directory},
        sandbox_directory{std::move(directory)}
    {
        normalise();
        check_invariant();
    }

    /*!\brief Construction from a given sandbox directory and a path within the sandbox directory.
     * \param sandbox_directory The absolute path to the sandbox directory.
     * \param path The relative or absolute path that must be inside the sandbox directory.
     */
    explicit sandboxed_path(std::filesystem::path sandbox_directory, std::filesystem::path path) :
        std::filesystem::path{std::move(path)},
        sandbox_directory{std::move(sandbox_directory)}
    {
        normalise();
        check_invariant();
    }

    sandboxed_path() = delete;                            //!< Deleted.
    sandboxed_path(sandboxed_path const &) = default;     //!< Defaulted.
    sandboxed_path(sandboxed_path &&) noexcept = default; //!< Defaulted.
    ~sandboxed_path() = default;                          //!< Defaulted.

    /*!\brief Replaces the path with a new path.
     * \param  new_path The new path.
     *
     * This works the same way as std::filesystem::path::operator=
     * and additionally checks the invariant.
     *
     * ### Exceptions
     *
     * Basic exception guarantee.
     */
    sandboxed_path & operator=(sandboxed_path const & new_path)
    {
        std::filesystem::path::operator=(new_path);
        check_invariant();
        return *this;
    }

    /*!\brief Replaces the path with a new path.
     * \param  new_path The new path.
     *
     * This works the same as std::filesystem::path::operator=
     * and additionally checks the invariant.
     *
     * ### Exceptions
     *
     * Basic exception guarantee
     */
    sandboxed_path & operator=(sandboxed_path && new_path)
    {
        std::filesystem::path::operator=(std::move(new_path));
        check_invariant();
        return *this;
    }
    //!\}

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
        auto normalised_path = (sandbox_directory / *this).lexically_normal();
        std::filesystem::path::operator=(normalised_path);
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

        auto rel_path = lexically_relative(sandbox_directory);

        // Leaving the temporary directory is not allowed.
        if (rel_path.string().find("..") == 0)
        {
            throw std::filesystem::filesystem_error("Leaving temporary directory is not allowed!",
                                                    sandbox_directory,
                                                    *this,
                                                    std::make_error_code(std::errc::invalid_argument));
        }
    }

public:
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

    /*!\brief Replaces the path with a new path.
     * \tparam path_t The type of the new_path.
     * \param  new_path The new path.
     * \return Reference to the sandboxed_path.
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
     * \return Reference to the sandboxed_path.
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
     * \return Reference to the sandboxed_path.
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
     * \return Reference to the sandboxed_path.
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
     * \return Reference to the sandboxed_path.
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
     * \return Reference to the sandboxed_path.
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
     * \return Reference to the sandboxed_path.
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
     * \return Reference to the sandboxed_path.
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
     * \return Reference to the sandboxed_path.
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
     * \param filename The replacement filename.
     * \return Reference to the sandboxed_path.
     *
     * This works the same as std::filesystem::path::replace_filename
     * and additionally checks the invariant.
     *
     * ### Exceptions
     *
     * Basic exception guarantee
     */
    sandboxed_path & replace_filename(std::filesystem::path const & filename)
    {
        std::filesystem::path::replace_filename(filename);
        normalise();
        check_invariant();
        return *this;
    }

    /*!\brief Replaces the extension.
     * \param replacement The replacement extension.
     * \return Reference to the sandboxed_path.
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
     * \return A sandboxed_path pointing to the parent path.
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

    /**!\brief Swaps the current path
     * \param other Reference to a sandboxed_path.
     */
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
sandboxed_path operator/(sandboxed_path lhs, Rhs const & rhs)
{
    lhs /= rhs;
    return lhs;
}

} // namespace seqan3::test
