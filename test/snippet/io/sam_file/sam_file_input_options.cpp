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

static constexpr auto sam_file_raw = R"(@HD	VN:1.6	pb:5.0.0
@SQ	SN:ref	LN:34
)";

static auto get_sam_file_input()
{
    return seqan3::sam_file_input{std::istringstream{sam_file_raw}, seqan3::format_sam{}};
}

void defaults_to_cerr()
{
    auto fin = get_sam_file_input();
    std::cerr << "Written to cerr: ";
    auto it = fin.begin(); // Prints to cerr: "Unsupported SAM header tag in @HD: pb"
}

void redirect_to_cout()
{
    auto fin = get_sam_file_input();
    fin.options.stream_warnings_to = std::addressof(std::cout); // Equivalent to `= &std::cout;`
    std::cout << "Written to cout: ";
    auto it = fin.begin(); // Prints to cout: "Unsupported SAM header tag in @HD: pb"
}

void redirect_to_file()
{
    temporary_file tmp_file{};
    auto fin = get_sam_file_input();

    { // Inner scope to close file before reading
        std::ofstream warning_file{tmp_file.path};
        fin.options.stream_warnings_to = std::addressof(warning_file); // Equivalent to `= &warning_file;`
        auto it = fin.begin(); // Prints to file: "Unsupported SAM header tag in @HD: pb"
    }

    std::cout << "Written to file: " << tmp_file.read_content();
}

void silence_warnings()
{
    auto fin = get_sam_file_input();
    fin.options.stream_warnings_to = nullptr;
    auto it = fin.begin(); // No warning emitted
}

int main()
{
    defaults_to_cerr();
    redirect_to_cout();
    redirect_to_file();
    silence_warnings();
}
