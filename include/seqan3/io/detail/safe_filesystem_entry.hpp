// -----------------------------------------------------------------------------------------------------
// Copyright (c) 2006-2020, Knut Reinert & Freie Universität Berlin
// Copyright (c) 2016-2020, Knut Reinert & MPI für molekulare Genetik
// This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
// shipped with this file and also available at: https://github.com/seqan/seqan3/blob/master/LICENSE.md
// -----------------------------------------------------------------------------------------------------

/*!\file
 * \brief Provides seqan3::detail::safe_filesystem_entry.
 * \author Rene Rahn <rene.rahn AT fu-berlin.de>
 */

#pragma once

#include <cassert>
#include <system_error>

#include <seqan3/core/platform.hpp>
#include <seqan3/std/filesystem>

namespace seqan3::detail
{

/*!\brief A safe guard to manage a filesystem entry, e.g. a file or a directory.
 * \ingroup io
 *
 * \details
 *
 * This raii-wrapper class allows for a safe removal of a created filesystem entry such as a directory or file.
 * This wrapper class assumes owning semantics. It is not copy-constructible or copy-assignable. In order to
 * prevent misuse also the default constructor is deleted.
 *
 * The following example demonstrates the use case.
 *
 * \include test/snippet/io/detail/safe_filesystem_entry_snippet.cpp
 */
class safe_filesystem_entry
{
public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    safe_filesystem_entry() = delete;                                             //!< Deleted.
    safe_filesystem_entry(safe_filesystem_entry const &) = delete;                //!< Deleted.
    safe_filesystem_entry(safe_filesystem_entry &&) = default;                    //!< Defaulted.
    safe_filesystem_entry & operator=(safe_filesystem_entry const &) = delete;    //!< Deleted.
    safe_filesystem_entry & operator=(safe_filesystem_entry &&) = default;        //!< Defaulted.

    /*!\brief Constructs the safe guard from a std::filesystem::path.
     * \param p The path pointing to a filesystem entry.
     */
    safe_filesystem_entry(std::filesystem::path p) : entry(std::move(p))
    {}

    //!\brief Calls std::filesystem::remove_all on the wrapped entry.
    ~safe_filesystem_entry()
    {
        std::error_code ec;
        std::filesystem::remove_all(entry, ec);

        assert(!static_cast<bool>(ec));
    }
    //!\}

    /*!\brief Removes a file or empty directory.
     * \returns `true` if the file was deleted, `false` if it did not exist.
     * \throws std::filesystem::filesystem_error on underlying OS API errors.
     *
     * \details
     *
     * Internally calls std::filesystem::remove on the stored std::filesystem::path.
     */
    bool remove()
    {
        return std::filesystem::remove(entry);
    }

    //!\copydoc remove
    bool remove_no_throw() const noexcept
    {
        std::error_code ec;
        return std::filesystem::remove(entry, ec);
    }

    /*!\brief Removes a file or directory and all its contents, recursively.
    * \returns Returns the number of files and directories that were deleted (which may be zero if p did not exist to
    *          begin with).
    * \throws std::filesystem::filesystem_error on underlying OS API errors.
    *
    * \details
    *
    * Internally calls std::filesystem::remove_all on the stored std::filesystem::path.
    */
    std::uintmax_t remove_all()
    {
        return std::filesystem::remove_all(entry);
    }

private:

    //!\brief The managed resource.
    std::filesystem::path entry;
};

} // namespace seqan3::detail
