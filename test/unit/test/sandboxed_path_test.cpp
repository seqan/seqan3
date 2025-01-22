// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/test/sandboxed_path.hpp>

namespace fs = std::filesystem;
using sandboxed_path = seqan3::test::sandboxed_path;

// initialization
TEST(sandboxed_path_init, init)
{
    // Valid paths
    EXPECT_NO_THROW(sandboxed_path("/dir"));
    EXPECT_NO_THROW(sandboxed_path("/dir", "/dir"));
    EXPECT_NO_THROW(sandboxed_path("/dir", "/dir/anotherdir"));
    EXPECT_NO_THROW(sandboxed_path("/dir", "."));
    EXPECT_NO_THROW(sandboxed_path("/dir", "./anotherdir"));
    EXPECT_NO_THROW(sandboxed_path("/dir", "../dir"));
    EXPECT_NO_THROW(sandboxed_path("/dir", "../dir/anotherdir"));
    EXPECT_NO_THROW(sandboxed_path("/dir", "anotherdir/../../dir/someotherdir"));

    EXPECT_EQ(sandboxed_path("/dir"), "/dir");
    EXPECT_EQ(sandboxed_path("/dir", "/dir"), "/dir");
    EXPECT_EQ(sandboxed_path("/dir", "/dir/anotherdir"), "/dir/anotherdir");
    EXPECT_EQ(sandboxed_path("/dir", "."), "/dir/");
    EXPECT_EQ(sandboxed_path("/dir", "./anotherdir"), "/dir/anotherdir");
    EXPECT_EQ(sandboxed_path("/dir", "../dir"), "/dir");
    EXPECT_EQ(sandboxed_path("/dir", "../dir/anotherdir"), "/dir/anotherdir");
    EXPECT_EQ(sandboxed_path("/dir", "anotherdir/../../dir/someotherdir"), "/dir/someotherdir");

    // Valid paths with suffix "/"
    EXPECT_NO_THROW(sandboxed_path("/dir/"));
    EXPECT_NO_THROW(sandboxed_path("/dir/", "/dir"));
    EXPECT_NO_THROW(sandboxed_path("/dir", "/dir/"));
    EXPECT_NO_THROW(sandboxed_path("/dir/", "/dir/"));
    EXPECT_EQ(sandboxed_path("/dir/"), "/dir/");
    EXPECT_EQ(sandboxed_path("/dir/", "/dir"), "/dir");
    EXPECT_EQ(sandboxed_path("/dir", "/dir/"), "/dir/");
    EXPECT_EQ(sandboxed_path("/dir/", "/dir/"), "/dir/");

    // Leaving the sandbox path is not allowed
    EXPECT_THROW(sandboxed_path("/dir", "/"), fs::filesystem_error);
    EXPECT_THROW(sandboxed_path("/dir", ".."), fs::filesystem_error);
    EXPECT_THROW(sandboxed_path("/dir", "/dir/.."), fs::filesystem_error);
    EXPECT_THROW(sandboxed_path("/dir", "somedir/../.."), fs::filesystem_error);
    EXPECT_THROW(sandboxed_path(".", ""), fs::filesystem_error);

    // Leaving the root directory
    EXPECT_THROW(sandboxed_path("/dir", "../.."), fs::filesystem_error);

    // Misc
    EXPECT_THROW(sandboxed_path("", ""), fs::filesystem_error);
}

// constructors
TEST(sandboxed_path_copy_constructor, operator_copy_constructor)
{
    sandboxed_path sandbox("/dir", "/dir/anotherdir");

    EXPECT_NO_THROW(sandboxed_path{sandbox});
    EXPECT_NO_THROW(sandboxed_path{std::move(sandbox)});
}
//

// assign
TEST(sandboxed_path_operator_assign, operator_assign)
{
    sandboxed_path path{"/dir"};

    EXPECT_NO_THROW(path = "/dir/charptr");
    EXPECT_EQ(path, "/dir/charptr");

    EXPECT_NO_THROW(path = std::string("/dir/string"));
    EXPECT_EQ(path, "/dir/string");

    EXPECT_NO_THROW(path = fs::path("/dir/path"));
    EXPECT_EQ(path, "/dir/path");

    EXPECT_NO_THROW(path = sandboxed_path("/dir/sandboxed_path"));
    EXPECT_EQ(path, "/dir/sandboxed_path");

    // check invalid assignment
    EXPECT_THROW(path = "/dir2", fs::filesystem_error);
}

// assign
TEST(sandboxed_path_assign, assign1)
{
    sandboxed_path path{"/dir"};

    // check assignment
    EXPECT_NO_THROW(path.assign("/dir/dir2"));
    EXPECT_EQ(path, "/dir/dir2");

    // check invalid assignment
    EXPECT_THROW(path.assign("/invalidDir"), fs::filesystem_error);
}

TEST(sandboxed_path_assign, assign2)
{
    sandboxed_path path{"/dir"};

    // check assignment
    std::string s1{"/dir/dir2"};
    EXPECT_NO_THROW(path.assign(begin(s1), end(s1)));
    EXPECT_EQ(path, "/dir/dir2");

    // check invalid assignment
    std::string s2{"/invalidDir"};
    EXPECT_THROW(path.assign(begin(s2), end(s2)), fs::filesystem_error);
}

// operator=/, append
TEST(sandboxed_path_operator_append, operator_append)
{
    sandboxed_path path{"/dir"};

    // check append
    EXPECT_NO_THROW(path /= "dir2");
    EXPECT_EQ(path, "/dir/dir2");

    // check invalid append
    EXPECT_THROW(path /= "../..", fs::filesystem_error);
}

TEST(sandboxed_path_append, append1)
{
    sandboxed_path path{"/dir"};

    // check append
    EXPECT_NO_THROW(path.append("dir2"));
    EXPECT_EQ(path, "/dir/dir2");

    // check invalid append
    EXPECT_THROW(path.append("../.."), fs::filesystem_error);
}

TEST(sandboxed_path_append, operator_append2)
{
    sandboxed_path path{"/dir"};

    // check append
    std::string s1{"dir2"};
    EXPECT_NO_THROW(path.append(begin(s1), end(s1)));
    EXPECT_EQ(path, "/dir/dir2");

    // check invalid append
    std::string s2{"../.."};
    EXPECT_THROW(path.append(begin(s2), end(s2)), fs::filesystem_error);
}

// concat
TEST(sandboxed_path_operator_concat, operator_concat)
{
    sandboxed_path path{"/dir"};

    // check concat
    EXPECT_NO_THROW(path += "/dir2");
    EXPECT_EQ(path, "/dir/dir2");

    // check invalid concat
    EXPECT_THROW(path += "/../..", fs::filesystem_error);
}

TEST(sandboxed_path_concat, concat1)
{
    sandboxed_path path{"/dir"};

    // check concat
    EXPECT_NO_THROW(path.concat("/dir2"));
    EXPECT_EQ(path, "/dir/dir2");

    // check invalid concat
    EXPECT_THROW(path.concat("/../.."), fs::filesystem_error);
}

TEST(sandboxed_path_concat, concat2)
{
    sandboxed_path path{"/dir"};

    // check concat
    std::string s1{"/dir2"};
    EXPECT_NO_THROW(path.concat(begin(s1), end(s1)));
    EXPECT_EQ(path, "/dir/dir2");

    // check invalid concat
    std::string s2{"/../.."};
    EXPECT_THROW(path.concat(begin(s2), end(s2)), fs::filesystem_error);
}

// remove_filename
TEST(sandboxed_path_remove_filename, remove_filename)
{
    sandboxed_path path{"/dir", "/dir/dir2/dir3"};

    EXPECT_EQ(path, "/dir/dir2/dir3");
    EXPECT_NO_THROW(path.remove_filename());
    EXPECT_EQ(path, "/dir/dir2/");

    path = "/dir";
    EXPECT_THROW(path.remove_filename(), fs::filesystem_error);
}

// replace_filename
TEST(sandboxed_path_replace_filename, replace_filename)
{

    sandboxed_path path{"/dir", "/dir/dir2"};

    EXPECT_EQ(path, "/dir/dir2");
    EXPECT_NO_THROW(path.replace_filename("dir3"));
    EXPECT_EQ(path, "/dir/dir3");

    path = "/dir";
    EXPECT_THROW(path.replace_filename("invalidDir"), fs::filesystem_error);
}

// replace_extension
TEST(sandboxed_path_replace_extension, replace_extension)
{
    sandboxed_path path{"/dir", "/dir/file.txt"};
    EXPECT_EQ(path, "/dir/file.txt");
    EXPECT_NO_THROW(path.replace_extension("doc"));
    EXPECT_EQ(path, "/dir/file.doc");
}

TEST(sandboxed_path_replace_extension, replace_extension2)
{
    sandboxed_path path{"/dir.txt", "/dir.txt"};
    EXPECT_EQ(path, "/dir.txt");
    EXPECT_THROW(path.replace_extension("doc"), fs::filesystem_error);
}

// parent_path
TEST(sandboxed_path_parent_path, parent_path)
{
    sandboxed_path path{"/dir", "/dir/dir2/dir3"};
    EXPECT_EQ(path, "/dir/dir2/dir3");

    EXPECT_NO_THROW(path = path.parent_path());
    EXPECT_EQ(path, "/dir/dir2");

    EXPECT_NO_THROW(path = path.parent_path());
    EXPECT_EQ(path, "/dir");

    EXPECT_THROW(path.parent_path(), fs::filesystem_error);
}

// swap
TEST(sandboxed_path_swap, swap)
{
    sandboxed_path path1{"/dir", "/dir/dir2/dir3"};
    sandboxed_path path2{"/dir", "/dir/dir_abc"};

    // check swap works
    EXPECT_NO_THROW(path1.swap(path2));
    EXPECT_EQ(path1, "/dir/dir_abc");
    EXPECT_EQ(path2, "/dir/dir2/dir3");

    // checks swap triggers invariant violation
    sandboxed_path path3{"/dir/dir2", "/dir/dir2/hallo"};
    EXPECT_THROW(path1.swap(path3), fs::filesystem_error);
}

// operator/
TEST(sandboxed_path_free_operator_append, free_operator_append)
{
    sandboxed_path path{"/dir"};
    EXPECT_NO_THROW(path = path / "dir2" / "dir3");
    EXPECT_EQ(path, "/dir/dir2/dir3");

    EXPECT_THROW(path / "../../../", fs::filesystem_error);
}

// Test special case when symbolic links are involved
TEST(sandboxed_path_symbolic_link, symbolic_link)
{
    // We create a symbolic link from /tmp/seqan3_sandboxed_path_symbolic_link_test -> /tmp
    auto tmp_base_dir = std::filesystem::temp_directory_path();
    auto tmp_dir = tmp_base_dir / "seqan3_sandboxed_path_symbolic_link_test";

    // if link already exists, remove it
    if (std::filesystem::exists(tmp_dir))
    {
        std::filesystem::remove_all(tmp_dir);
    }

    // create symlink and a sandboxed_path
    std::filesystem::create_directory_symlink(tmp_base_dir, tmp_dir);
    sandboxed_path path{tmp_dir};

    // check the sandboxed path did not resolve the symlink
    EXPECT_EQ(tmp_dir, path);

    // Cleanup, remove link
    std::filesystem::remove_all(tmp_dir);
}
