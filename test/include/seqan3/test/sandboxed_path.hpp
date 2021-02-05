// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/* \file
 * \brief Provides the sandboxed_path and surrounding free functions
 * \author Simon Gene <simon.gottlieb AT fu-berlin.de>
 */

#pragma once

#include <seqan3/core/platform.hpp>
#include <seqan3/std/filesystem>

namespace seqan3::test
{

/** \brief sandboxed_path is a path that can not leave a certain path
 *
 * sandboxed_path inherits from std::filesystem::path and behaves mostly
 * like it. In addition it receives a sandbox directory at construction time.
 * Functions are overloaded and a check for the invariant is added.
 *
 *  Invariant:
 *  - sandboxed_path is always converted to an absolute path
 *  - sandboxed_path always points to a file inside a given sandbox directory.
 *  - sandbox directory is unmutable during the life cycle of a sandboxed_path
 *
 * Caveat:
 *  - relatives paths are not possible
 *  - sandboxed_path is not 100% tight
 *    - calling relative_path() leaves the environment of sandboxed_path
 */
class sandboxed_path : public std::filesystem::path
{
private:
    std::filesystem::path const sandbox_path;

public:
    /*!\brief Construction of a sandboxed_path
     * \param path must be an absolute path.
     *
     * A sandboxed_path initialized with this constructor will
     * point to path and disallow extension that leave path.
     */
    explicit sandboxed_path(std::filesystem::path path)
        : std::filesystem::path{path}
        , sandbox_path{std::move(path)}
    {
        normalize();
        checkInvariant();
    }

    /*!\brief Construction of a sandboxed path
     * \param sandbox_path must be an absolute path.
     * \param path must be a path that is inside of sandbox_path.
     *
     * path is allowed to be relative
     */
    explicit sandboxed_path(std::filesystem::path sandbox_path, std::filesystem::path path)
        : std::filesystem::path {path}
        , sandbox_path             {std::move(sandbox_path)}
    {
        normalize();
        checkInvariant();
    }

    sandboxed_path() = delete; //!< Default construction not possible
    sandboxed_path(sandboxed_path const&) = default;    //!< Default copy constructor
    sandboxed_path(sandboxed_path&&) noexcept = default; //!< Default move constructor


private:
    void normalize()
    {
        auto concated_path   = sandbox_path / std::filesystem::path{*this};
        auto normalized_path = concated_path.lexically_normal();
        std::filesystem::path::operator=(normalized_path);
    }

    void checkInvariant() const
    {
        // Checking that sandbox_path is an absolute path
        if (!sandbox_path.is_absolute()) {
            throw std::filesystem::filesystem_error("sandbox path must be an absolute path",
                                                sandbox_path, relative(sandbox_path),
                                                std::make_error_code(std::errc::invalid_argument));
        }

        auto rel_path = lexically_relative(sandbox_path);

        // Leaving the temporary directory is not allowed.
        if (rel_path.string().find("..") == 0) {
            throw std::filesystem::filesystem_error("Leaving temporary directory is not allowed!",
                                                sandbox_path, relative(sandbox_path),
                                                std::make_error_code(std::errc::invalid_argument));
        }
    }

public:
    /*!\brief Replaces the path with a new path
     *
     * This works the same as std::filesystem::path::operator=
     * and additionally checks the invariant.
     *
     * Basic exception guarantee
     */
    template <class Path>
    sandboxed_path & operator=(Path const & path)
    {
        std::filesystem::path::operator=(path);
        normalize();
        checkInvariant();
        return *this;
    }

    /*!\brief Replaces the path with a new path
     *
     * This works the same as std::filesystem::path::operator=
     * and additionally checks the invariant.
     *
     * Basic exception guarantee
     */
    template <class Path>
    sandboxed_path & operator=(Path && path)
    {
        std::filesystem::path::operator=(std::move(path));
        normalize();
        checkInvariant();
        return *this;
    }


    /*!\brief Replaces the path with a new path
     *
     * This works the same as std::filesystem::path::assign
     * and additionally checks the invariant.
     *
     * Basic exception guarantee
     */
    template <typename Source>
    sandboxed_path & assign(Source const & source)
    {
        std::filesystem::path::assign(source);
        normalize();
        checkInvariant();
        return *this;
    }

    /*!\brief Replaces the path with a new path
     *
     * This works the same as std::filesystem::path::assign
     * and additionally checks the invariant.
     *
     * Basic exception guarantee
     */
    template <typename InputIt>
    sandboxed_path & assign(InputIt first, InputIt last)
    {
        std::filesystem::path::assign(first, last);
        normalize();
        checkInvariant();
        return *this;
    }

    /*!\brief Extends the path
     *
     * This works the same as std::filesystem::path::operator/=
     * and additionally checks the invariant.
     *
     * Basic exception guarantee
     */
    template <typename Source>
    sandboxed_path & operator/=(Source const & source)
    {
        return append(source);
    }

    /*!\brief Extends the path
     *
     * This works the same as std::filesystem::path::append
     * and additionally checks the invariant.
     *
     * Basic exception guarantee
     */
    template <typename Source>
    sandboxed_path & append(Source const & source)
    {
        std::filesystem::path::append(source);
        normalize();
        checkInvariant();
        return *this;
    }

    /*!\brief Extends the path
     *
     * This works the same as std::filesystem::path::append
     * and additionally checks the invariant.
     *
     * Basic exception guarantee
     */
    template <typename InputIter>
    sandboxed_path & append(InputIter first, InputIter second)
    {
        std::filesystem::path::append(first, second);
        normalize();
        checkInvariant();
        return *this;
    }

    /*!\brief Extends the path
     *
     * This works the same as std::filesystem::path::operator+=
     * and additionally checks the invariant.
     *
     * Basic exception guarantee
     */
    template <typename Source>
    sandboxed_path & operator+=(Source const & source)
    {
        return concat(source);
    }

    /*!\brief Extends the path
     *
     * This works the same as std::filesystem::path::concat
     * and additionally checks the invariant.
     *
     * Basic exception guarantee
     */
    template <typename Source>
    sandboxed_path & concat(Source const & source)
    {
        std::filesystem::path::concat(source);
        normalize();
        checkInvariant();
        return *this;
    }

    /*!\brief Extends the path
     *
     * This works the same as std::filesystem::path::concat
     * and additionally checks the invariant.
     *
     * Basic exception guarantee
     */
    template <typename InputIt>
    sandboxed_path & concat(InputIt first, InputIt last)
    {
        std::filesystem::path::concat(first, last);
        normalize();
        checkInvariant();
        return *this;
    }

    /*!\brief Removes the filename
     *
     * This works the same as std::filesystem::path::remove_filename
     * and additionally checks the invariant.
     *
     * Basic exception guarantee
     */
    sandboxed_path & remove_filename()
    {
        std::filesystem::path::remove_filename();
        normalize();
        checkInvariant();
        return *this;
    }

    /*!\brief Replaces the filename
     *
     * This works the same as std::filesystem::path::replace_filename
     * and additionally checks the invariant.
     *
     * Basic exception guarantee
     */
    sandboxed_path & replace_filename(std::filesystem::path const & p)
    {
        std::filesystem::path::replace_filename(p);
        normalize();
        checkInvariant();
        return *this;
    }

    /*!\brief Replaces the extension
     *
     * This works the same as std::filesystem::path::replace_extension
     * and additionally checks the invariant.
     *
     * Basic exception guarantee
     */
    sandboxed_path & replace_extension(std::filesystem::path const & replacement = std::filesystem::path{})
    {
        std::filesystem::path::replace_extension(replacement);
        normalize();
        checkInvariant();
        return *this;
    }

    /*!\brief Returns sandboxed path to the parent path
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


    void swap(sandboxed_path &) noexcept = delete; //!< Swap operator is not possible

private:
    void clear() = delete; //!< Unlikly operator
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
