<!--
    SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
    SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
    SPDX-License-Identifier: CC-BY-4.0
-->

# Continuous on Pull Request

## ci_cmake.yml

Runs CMake specific tests, e.g., using SeqAn3 as external project.

## ci_coverage.yml

Runs coverage test and uploads the report to Codecov.

## ci_documentation.yml

Builds the documentation and checks for failures.

## ci_lint.yml

Runs clang-format and cmake-format. Resulting changes will be pushed as a new commit.

Other CI only starts after linting does not produce changes.

Only runs on PRs.

## ci_linux.yml

Runs unit tests on Ubuntu.

## ci_macos.yml

Runs unit tests on macOS.

## ci_misc.yml

Runs snippet, performance, and header tests on Ubuntu.

# Continuous on Push

These run in addition to *Continuous on Pull Request*.

## deploy_documentation.yaml

Builds documentation and uploads it to https://docs.seqan.de/.

Uses repository secrets for the target server and authentication.

## update_cookbook.yml

Checks if any new snippets are to be added to the Cookbook. Will push a new commit if new snippets are found.

# Cron

## cron_api.yml

Runs the [API-Stability](https://github.com/seqan/seqan3/blob/main/test/api_stability/README.md) test.

In case of failure, creates an issue containing error logs.

## cron_avx2.yml

Will run all test suites on all compilers with AVX2 enabled.

In case of failure, creates an issue containing error logs.

## cron_latest_libraries.yml

Will update all submodules, gbenchmark and googletest to respective main versions. Then runs all test suites on all
compilers.

In case of failure, creates an issue containing error logs.

# On Demand

## hide_comments.yml

Can be used to hide bot comments in an issue. Intended for use with the CRON issues.

## ram_usage.yml

Tracks the RAM-Usage when compiling and generates a CSV artifact.

This workflow has two optional inputs:
  * The GCC version (default: `13`)
  * The CXX_FLAGS (default: None)
