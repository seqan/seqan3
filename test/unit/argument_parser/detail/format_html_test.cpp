// SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
// SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/argument_parser/detail/format_html.hpp>
#include <seqan3/utility/char_operations/predicate.hpp>

TEST(html_format, empty_information)
{
    std::string my_stdout;
    std::string expected;

    // Empty html help page.
    char const * argv0[] = {"./help_add_test --version-check false", "--export-help", "html"};
    seqan3::argument_parser parser0{"empty_options", 3, argv0};
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser0.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    my_stdout = testing::internal::GetCapturedStdout();
    expected =
        std::string("<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01//EN\" http://www.w3.org/TR/html4/strict.dtd\">\n"
                    "<html lang=\"en\">\n"
                    "<head>\n"
                    "<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\">\n"
                    "<title>empty_options &mdash; </title>\n"
                    "</head>\n"
                    "<body>\n"
                    "<h1>empty_options</h1>\n"
                    "<div></div>\n"
                    "<h2>Options</h2>\n"
                    "<h3>Basic options:</h3>\n"
                    "<dl>\n"
                    "<dt><strong>-h</strong>, <strong>--help</strong></dt>\n"
                    "<dd>Prints the help page.</dd>\n"
                    "<dt><strong>-hh</strong>, <strong>--advanced-help</strong></dt>\n"
                    "<dd>Prints the help page including advanced options.</dd>\n"
                    "<dt><strong>--version</strong></dt>\n"
                    "<dd>Prints the version information.</dd>\n"
                    "<dt><strong>--copyright</strong></dt>\n"
                    "<dd>Prints the copyright/license information.</dd>\n"
                    "<dt><strong>--export-help</strong> (std::string)</dt>\n"
                    "<dd>Export the help page information. Value must be one of [html, man].</dd>\n"
                    "<dt><strong>--version-check</strong> (bool)</dt>\n"
                    "<dd>Whether to check for the newest app version. Default: true.</dd>\n"
                    "</dl>\n"
                    "<h2>Version</h2>\n"
                    "<p>\n"
                    "<strong>Last update: </strong>\n"
                    "<br>\n"
                    "<strong>empty_options version: </strong>\n"
                    "<br>\n"
                    "<strong>SeqAn version: </strong>"
                    + std::string{seqan3::seqan3_version_cstring}
                    + "\n"
                      "<br>\n"
                      "</p>\n"
                      "</body></html>");
    EXPECT_EQ(my_stdout, expected);

    char const * argv1[] = {"./help_add_test --version-check false", "--export-help=html"};
    seqan3::argument_parser parser1{"empty_options", 2, argv1};
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser1.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");
    my_stdout = testing::internal::GetCapturedStdout();
    EXPECT_EQ(my_stdout, expected);
}

TEST(html_format, full_information_information)
{
    std::string my_stdout;
    std::string expected;
    int option_value{5};
    bool flag_value{false};
    int8_t non_list_pos_opt_value{1};
    std::vector<std::string> list_pos_opt_value{};

    // Full html help page.
    char const * argv0[] = {"./help_add_test --version-check false", "--export-help", "html"};
    seqan3::argument_parser parser1{"program_full_options", 3, argv0};
    parser1.info.synopsis.push_back("./some_binary_name synopsis");
    parser1.info.synopsis.push_back("./some_binary_name synopsis2");
    parser1.info.description.push_back("description");
    parser1.info.description.push_back("description2");
    parser1.info.short_description = "short description";
    parser1.info.url = "https://seqan.de";
    parser1.info.short_copyright = "short copyright";
    parser1.info.long_copyright = "long_copyright";
    parser1.info.citation = "citation";
    parser1.info.author = "author";
    parser1.info.email = "email";
    parser1.add_option(option_value, 'i', "int", "this is a int option.");
    parser1.add_option(option_value, 'j', "jint", "this is a required int option.", seqan3::option_spec::required);
    parser1.add_flag(flag_value, 'f', "flag", "this is a flag.");
    parser1.add_flag(flag_value, 'k', "kflag", "this is a flag.");
    parser1.add_positional_option(non_list_pos_opt_value, "this is a positional option.");
    parser1.add_positional_option(list_pos_opt_value, "this is a positional option.");
    parser1.info.examples.push_back("example");
    parser1.info.examples.push_back("example2");
    testing::internal::CaptureStdout();
    EXPECT_EXIT(parser1.parse(), ::testing::ExitedWithCode(EXIT_SUCCESS), "");

    my_stdout = testing::internal::GetCapturedStdout();
    expected = std::string(
        "<!DOCTYPE html PUBLIC \"-//W3C//DTD HTML 4.01//EN\" http://www.w3.org/TR/html4/strict.dtd\">\n"
        "<html lang=\"en\">\n"
        "<head>\n"
        "<meta http-equiv=\"content-type\" content=\"text/html; charset=utf-8\">\n"
        "<title>program_full_options &mdash; short description</title>\n"
        "</head>\n"
        "<body>\n"
        "<h1>program_full_options</h1>\n"
        "<div>short description</div>\n"
        "<h2>Synopsis</h2>\n"
        "<p>\n"
        "<strong>./some_binary_name</strong> synopsis\n"
        "<br>\n"
        "<strong>./some_binary_name</strong> synopsis2\n"
        "<br>\n"
        "</p>\n"
        "<h2>Description</h2>\n"
        "<p>\n"
        "description\n"
        "</p>\n"
        "<p>\n"
        "description2\n"
        "</p>\n"
        "<h2>Positional Arguments</h2>\n"
        "<dl>\n"
        "<dt><strong>ARGUMENT-1</strong> (<em>signed 8 bit integer</em>)</dt>\n"
        "<dd>this is a positional option. </dd>\n"
        "<dt><strong>ARGUMENT-2</strong> (<em>List</em> of <em>std::string</em>)</dt>\n"
        "<dd>this is a positional option. Default: []. </dd>\n"
        "</dl>\n"
        "<h2>Options</h2>\n"
        "<h3>Basic options:</h3>\n"
        "<dl>\n"
        "<dt><strong>-h</strong>, <strong>--help</strong></dt>\n"
        "<dd>Prints the help page.</dd>\n"
        "<dt><strong>-hh</strong>, <strong>--advanced-help</strong></dt>\n"
        "<dd>Prints the help page including advanced options.</dd>\n"
        "<dt><strong>--version</strong></dt>\n"
        "<dd>Prints the version information.</dd>\n"
        "<dt><strong>--copyright</strong></dt>\n"
        "<dd>Prints the copyright/license information.</dd>\n"
        "<dt><strong>--export-help</strong> (std::string)</dt>\n"
        "<dd>Export the help page information. Value must be one of [html, man].</dd>\n"
        "<dt><strong>--version-check</strong> (bool)</dt>\n"
        "<dd>Whether to check for the newest app version. Default: true.</dd>\n"
        "<dt><strong>-i</strong>, <strong>--int</strong> (<em>signed 32 bit integer</em>)</dt>\n"
        "<dd>this is a int option. Default: 5. </dd>\n"
        "<dt><strong>-j</strong>, <strong>--jint</strong> (<em>signed 32 bit integer</em>)</dt>\n"
        "<dd>this is a required int option. </dd>\n"
        "<dt><strong>-f</strong>, <strong>--flag</strong></dt>\n"
        "<dd>this is a flag.</dd>\n"
        "<dt><strong>-k</strong>, <strong>--kflag</strong></dt>\n"
        "<dd>this is a flag.</dd>\n"
        "</dl>\n"
        "<h2>Examples</h2>\n"
        "<p>\n"
        "example\n"
        "</p>\n"
        "<p>\n"
        "example2\n"
        "</p>\n"
        "<h2>Version</h2>\n"
        "<p>\n"
        "<strong>Last update: </strong>\n"
        "<br>\n"
        "<strong>program_full_options version: </strong>\n"
        "<br>\n"
        "<strong>SeqAn version: </strong>"
        + std::string{seqan3::seqan3_version_cstring}
        + "\n"
          "<br>\n"
          "</p>\n"
          "<h2>Url</h2>\n"
          "<p>\n"
          "https://seqan.de\n"
          "<br>\n"
          "</p>\n"
          "<h2>Legal</h2>\n"
          "<p>\n"
          "<strong>program_full_options Copyright: </strong>short copyright\n"
          "<br>\n"
          "<strong>Author: </strong>author\n"
          "<br>\n"
          "<strong>Contact: </strong>email\n"
          "<br>\n"
          "<strong>SeqAn Copyright: </strong>2006-2025 Knut Reinert, FU-Berlin; released under the 3-clause BSDL.\n"
          "<br>\n"
          "<strong>In your academic works please cite: </strong>citation\n"
          "<br>\n"
          "For full copyright and/or warranty information see <strong>--copyright</strong>.\n"
          "<br>\n"
          "</p>\n"
          "</body></html>");
    EXPECT_EQ(my_stdout, expected);
}

TEST(export_help, parse_error)
{
    char const * argv[] = {"./help_add_test --version-check false", "--export-help"};
    char const * argv2[] = {"./help_add_test --version-check false", "--export-help=atml"};
    char const * argv3[] = {"./help_add_test --version-check false", "--export-help", "atml"};

    // no value after --export-help
    EXPECT_THROW((seqan3::argument_parser{"test_parser", 2, argv}), seqan3::argument_parser_error);

    // wrong value after --export-help
    EXPECT_THROW((seqan3::argument_parser{"test_parser", 2, argv2}), seqan3::validation_error);

    // wrong value after --export-help
    EXPECT_THROW((seqan3::argument_parser{"test_parser", 3, argv3}), seqan3::validation_error);
}
