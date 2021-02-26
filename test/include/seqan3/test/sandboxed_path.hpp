// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/* \file
 * \brief Provides the sandboxed_path and related free functions
 * \author Simon Gene <simon.gottlieb AT fu-berlin.de>
 */

#pragma once

#include <seqan3/std/filesystem>

#include <seqan3/core/platform.hpp>

namespace seqan3::test
{

/** \brief Utility class to stay inside a sandbox path.
 *
 * sandboxed_path inherits from std::filesystem::path and behaves mostly
 * like it. In addition, it receives a sandbox directory at construction time.
 * Functions are overloaded and a check for the invariant is added.
 *
 *  Invariant:
 *  - sandboxed_path is always converted to an absolute path
 *  - sandboxed_path always points to a file inside a given sandbox directory.
 *  - sandbox directory is immutable during the life cycle of a sandboxed_path
 *
 * Caveat:
 *  - relative paths are not possible
 *  - some functions will leave the sandboxed environment.
 *    - e.g.: calling relative_path() leaves the environment of sandboxed_path.
 */
class sandboxed_path : public std::filesystem::path
{
private:
    std::filesystem::path const sandbox_path;

public:
    /*!\brief Construction of a sandboxed_path.
     * \param path must be an absolute path.
     *
     * A sandboxed_path initialised with this constructor will
     * point to `path` and disallow extension that leave `path`.
     */
    explicit sandboxed_path(std::filesystem::path path)
        : std::filesystem::path{path}
        , sandbox_path{std::move(path)}
    {
        normalise();
        checkInvariant();
    }

    /*!\brief Construction of a sandboxed path.
     * \param sandbox_path must be an absolute path.
     * \param path must be a path that is inside of sandbox_path.
     *
     * path is allowed to be relative
     */
    explicit sandboxed_path(std::filesystem::path sandbox_path, std::filesystem::path path)
        : std::filesystem::path {path}
        , sandbox_path             {std::move(sandbox_path)}
    {
        normalise();
        checkInvariant();
    }

    sandboxed_path() = delete; //!< Deleted.
    sandboxed_path(sandboxed_path const&) = default; //!< Defaulted.
    sandboxed_path(sandboxed_path&&) noexcept = default; //!< Defaulted.
    ~sandboxed_path() = default; //!< Defaulted.

    /*!\brief Replaces the path with a new path.
     * \tparam path_t The type of the path.
     *
     * This works the same as std::filesystem::path::operator=
     *
     * This works the same way as std::filesystem::path::operator=
     * and additionally checks the invariant.
     *
     * ### Exceptions
     *
     * Basic exception guarantee.
     */
    template <typename path_t>
    sandboxed_path & operator=(path_t const & path)
    {
        std::filesystem::path::operator=(path);
        normalise();
        checkInvariant();
        return *this;
    }

    /*!\brief Replaces the path with a new path.
     *
     * This works the same as std::filesystem::path::operator=
     * and additionally checks the invariant.
     *
     * Basic exception guarantee
     */
    template <typename path_t>
    sandboxed_path & operator=(path_t && path)
    {
        std::filesystem::path::operator=(std::move(path));
        normalise();
        checkInvariant();
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
#if SEQAN3_WORKAROUND_GCC7_INCOMPLETE_FILESYSTEM
        // convert into an absolute path
        auto path = [this]() {
            if (is_relative()) {
                return (sandbox_path / *this).string();
            } else {
                return string();
            }
        }();
        // Now we need to manually collapse any "." and ".."
        size_t current_pos = path.find("/", 0); // find first "/" character
        std::vector<std::string> path_parts;
        while (current_pos < path.size()) {
            auto end_pos = path.find("/", current_pos + 1);
            auto word = path.substr(current_pos+1, end_pos-current_pos-1);
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
                    throw std::filesystem::filesystem_error("Path can not be normalised", *this,
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
        for (auto const& p : path_parts) {
            normalised_path += "/" + p;
        }
        std::filesystem::path::operator=(normalised_path.data());
#else
        auto normalised_path = std::filesystem::weakly_canonical(sandbox_path / *this);
        std::filesystem::path::operator=(normalised_path);
#endif
    }

public:
    /*!\brief Checks the invariant.
     *
     * Checks that the invariant of the class sandboxed_path is kept.
     * See class description for invariant.
     * \throws std::filesystem::filesystem_error if invariant was violated.
     */
    void checkInvariant() const
    {
        // Checking that sandbox_path is an absolute path
        if (!sandbox_path.is_absolute()) {
            throw std::filesystem::filesystem_error("sandbox path must be an absolute path",
                                                sandbox_path, *this,
                                                std::make_error_code(std::errc::invalid_argument));
        }

#if SEQAN3_WORKAROUND_GCC7_INCOMPLETE_FILESYSTEM
        auto current_dir = string();
        auto sandbox_dir = sandbox_path.string();

        // remove trailing '/' of path
        if (!sandbox_dir.empty() and sandbox_dir.back() == '/') {
            sandbox_dir.pop_back();
        }

        // Leaving the temporary directory is not allowed.
        bool starts_with_sandbox_dir = (current_dir.rfind(sandbox_dir, 0) == 0);
        if (!starts_with_sandbox_dir or
            (current_dir.size() > sandbox_dir.size() and current_dir.at(sandbox_dir.size()) != '/'))
        {
            throw std::filesystem::filesystem_error("Leaving temporary directory is not allowed!",
                                                sandbox_path, *this,
                                                std::make_error_code(std::errc::invalid_argument));
        }

#else
        auto rel_path = lexically_relative(sandbox_path);

        // Leaving the temporary directory is not allowed.
        if (rel_path.string().find("..") == 0) {
            throw std::filesystem::filesystem_error("Leaving temporary directory is not allowed!",
                                                sandbox_path, *this,
                                                std::make_error_code(std::errc::invalid_argument));
        }
#endif
    }

    /*!\brief Replaces the path with a new path.
     *
     * This works the same as std::filesystem::path::assign
     * and additionally checks the invariant.
     *
     * Basic exception guarantee
     */
    template <typename source_t>
    sandboxed_path & assign(source_t const & source)
    {
        std::filesystem::path::assign(source);
        normalise();
        checkInvariant();
        return *this;
    }

    /*!\brief Replaces the path with a new path.
     *
     * This works the same as std::filesystem::path::assign
     * and additionally checks the invariant.
     *
     * Basic exception guarantee
     */
    template <typename input_iter_t>
    sandboxed_path & assign(input_iter_t first, input_iter_t last)
    {
        std::filesystem::path::assign(first, last);
        normalise();
        checkInvariant();
        return *this;
    }

    /*!\brief Extends the path.
     *
     * This works the same as std::filesystem::path::operator/=
     * and additionally checks the invariant.
     *
     * Basic exception guarantee
     */
    template <typename source_t>
    sandboxed_path & operator/=(source_t const & source)
    {
        return append(source);
    }

    /*!\brief Extends the path.
     *
     * This works the same as std::filesystem::path::append
     * and additionally checks the invariant.
     *
     * Basic exception guarantee
     */
    template <typename source_t>
    sandboxed_path & append(source_t const & source)
    {
        std::filesystem::path::append(source);
        normalise();
        checkInvariant();
        return *this;
    }

    /*!\brief Extends the path.
     *
     * This works the same as std::filesystem::path::append
     * and additionally checks the invariant.
     *
     * Basic exception guarantee
     */
    template <typename input_iter_t>
    sandboxed_path & append(input_iter_t first, input_iter_t second)
    {
        std::filesystem::path::append(first, second);
        normalise();
        checkInvariant();
        return *this;
    }

    /*!\brief Extends the path.
     *
     * This works the same as std::filesystem::path::operator+=
     * and additionally checks the invariant.
     *
     * Basic exception guarantee
     */
    template <typename source_t>
    sandboxed_path & operator+=(source_t const & source)
    {
        return concat(source);
    }

    /*!\brief Extends the path.
     *
     * This works the same as std::filesystem::path::concat
     * and additionally checks the invariant.
     *
     * Basic exception guarantee
     */
    template <typename source_t>
    sandboxed_path & concat(source_t const & source)
    {
        std::filesystem::path::concat(source);
        normalise();
        checkInvariant();
        return *this;
    }

    /*!\brief Extends the path.
     *
     * This works the same as std::filesystem::path::concat
     * and additionally checks the invariant.
     *
     * Basic exception guarantee
     */
    template <typename input_iter_t>
    sandboxed_path & concat(input_iter_t first, input_iter_t last)
    {
        std::filesystem::path::concat(first, last);
        normalise();
        checkInvariant();
        return *this;
    }

    /*!\brief Removes the filename.
     *
     * This works the same as std::filesystem::path::remove_filename
     * and additionally checks the invariant.
     *
     * Basic exception guarantee
     */
    sandboxed_path & remove_filename()
    {
        std::filesystem::path::remove_filename();
        normalise();
        checkInvariant();
        return *this;
    }

    /*!\brief Replaces the filename.
     *
     * This works the same as std::filesystem::path::replace_filename
     * and additionally checks the invariant.
     *
     * Basic exception guarantee
     */
    sandboxed_path & replace_filename(std::filesystem::path const & p)
    {
        std::filesystem::path::replace_filename(p);
        normalise();
        checkInvariant();
        return *this;
    }

    /*!\brief Replaces the extension.
     *
     * This works the same as std::filesystem::path::replace_extension
     * and additionally checks the invariant.
     *
     * Basic exception guarantee
     */
    sandboxed_path & replace_extension(std::filesystem::path const & replacement = std::filesystem::path{})
    {
        std::filesystem::path::replace_extension(replacement);
        normalise();
        checkInvariant();
        return *this;
    }

    /*!\brief Returns sandboxed path to the parent path.
     *
     * This works the same as std::filesystem::path::parent_path
     * and additionally checks the invariant.
     *
     * Basic exception guarantee
     */
    sandboxed_path parent_path() const
    {
        auto parent_path = std::filesystem::path::parent_path();
        return sandboxed_path{sandbox_path, parent_path};
    }


    void swap(sandboxed_path & other)
    {
        std::filesystem::path::swap(other);
        checkInvariant();
        other.checkInvariant();
    }

    void clear() = delete; //!< Not implemented. Invariant requires the path to be an absolute path.
};

/** Free sandboxed_path append operator.
 *  This work the same as std::filesystem::operator/(std::filesystem::path&)
 *  and additionally checks the invariant.
 */
template <typename Rhs>
sandboxed_path operator/(sandboxed_path lhs, Rhs const& rhs) {
    lhs /= rhs;
    return lhs;
}

}
