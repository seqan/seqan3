#!/usr/bin/env python3

# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: CC0-1.0
#
# Usage process_compiler_error_log.py <log_file>
#
# Processes a build log and writes a GitHub-Markdown-formatted list of build errors to stdcout.
import html
import re
import sys

# https://regex101.com/r/aE0qX3/2
tokenise_regex = re.compile(r"(\[\s*\d+%\](?:.(?!\[\s*\d+%\]))+(?:: error:).*?(?=\[\s*\d+%\]|$))", re.DOTALL)
# Match line containg 'error:', but stop at linebreak, semicolon or parenthesis.
error_regex = re.compile(r"(error:[^(;\n]*)")

with open(sys.argv[1]) as file:
    content = file.read()

counter = 1
log = ''
# One `group` (a string) is all messages in the log associated with one failing compilation
for group in [match.groups()[0] for match in tokenise_regex.finditer(content)]:
    error_message = html.escape(re.search(error_regex, group).group().rstrip()[:110])
    log += '<details><summary>Error {}: <code>{}</code></summary>\n\n```text\n'.format(counter, error_message)
    log += '\n'.join(group.split('\n')[:30]).rstrip() # Only take first 30 lines
    log += '\n```\n</details>\n'
    counter += 1

# Truncate to 65300 to not exceed the limit of GitHub's API
print(log[:65300])
