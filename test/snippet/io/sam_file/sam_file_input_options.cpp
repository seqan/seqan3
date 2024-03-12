// SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: CC0-1.0

#include <sstream>

#include <seqan3/io/sam_file/input.hpp>

// A helper struct to create a temporary file and remove it when it goes out of scope.
struct temporary_file
{
    std::filesystem::path const path{std::filesystem::temp_directory_path() / "warnings.txt"};

    temporary_file()
    {
        std::ofstream file{path}; // Create file
    }
    temporary_file(temporary_file const &) = delete;
    temporary_file & operator=(temporary_file const &) = delete;
    temporary_file(temporary_file &&) = delete;
    temporary_file & operator=(temporary_file &&) = delete;
    ~temporary_file()
    {
        std::filesystem::remove(path);
    }

    std::string read_content() const
    {
        std::ifstream file{path};
        return std::string{std::istreambuf_iterator<char>{file}, std::istreambuf_iterator<char>{}};
    }
};

static constexpr auto sam_file_raw = R"(@HD	VN:1.6	pb:5.0.0	ot:ter
@SQ	SN:ref	LN:34
)";

static auto get_sam_file_input()
{
    return seqan3::sam_file_input{std::istringstream{sam_file_raw}, seqan3::format_sam{}};
}

void defaults_to_cerr()
{
    auto fin = get_sam_file_input();
    auto it = fin.begin();
}

void redirect_to_cout()
{
    auto fin = get_sam_file_input();
    fin.options.stream_warnings_to = std::addressof(std::cout); // Equivalent to `= &std::cout;`
    auto it = fin.begin();
}

void redirect_to_file()
{
    temporary_file tmp_file{};
    auto fin = get_sam_file_input();

    { // Inner scope to close file before reading
        std::ofstream warning_file{tmp_file.path};
        fin.options.stream_warnings_to = std::addressof(warning_file); // Equivalent to `= &warning_file;`
        auto it = fin.begin();
    }

    std::cout << "File content:\n" << tmp_file.read_content();
}

void silence_warnings()
{
    auto fin = get_sam_file_input();
    fin.options.stream_warnings_to = nullptr;
    auto it = fin.begin();
}

void filter()
{
    auto fin = get_sam_file_input();
    std::stringstream stream{};
    fin.options.stream_warnings_to = std::addressof(stream); // Equivalent to `= &stream;`
    auto it = fin.begin();

    for (std::string line{}; std::getline(stream, line);)
    {
        // If "pb" is not found in the warning, print it to cerr.
        if (line.find("pb") == std::string::npos) // C++23: `!line.contains("pb")`
            std::cerr << line << '\n';
    }
}

void print_section(std::string_view const section)
{
    std::cout << "### " << section << " ###\n";
    std::cerr << "### " << section << " ###\n";
}

int main()
{
    print_section("defaults_to_cerr");
    defaults_to_cerr();

    print_section("redirect_to_cout");
    redirect_to_cout();

    print_section("redirect_to_file");
    redirect_to_file();

    print_section("silence_warnings");
    silence_warnings();

    print_section("filter");
    filter();
}
