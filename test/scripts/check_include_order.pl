#!/usr/bin/env perl

# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

# Usage check_include_order.pl <file1> [<file2>] [<file3>] ....
# Will output the names of the files that have incorrect include order.
#
# Usage example: `check_include_order.pl $(git diff --name-only HEAD origin/main)`
# Checks all files that have been touch since last main for any incorrect include ordering
#
# Usage example: `check_include_order.pl $(find include -type f)
# Checks all files of seqan3 for incorrect include ordering
#
# This is a script checking for blocks of include statements:
#
#  #include <cassert>
#  #include <concepts>
#  #include <utility>
#
#  It will extract the include paths to
#     cassert
#     concepts
#     utility
#  and check if they appear in the correct order.
use strict;
use warnings;

# Loops through all input files
foreach my $file (@ARGV) {
    # if file doesn't exists, jump to next file
    next unless (open my $info, $file);

    # variable to track current include block
    my @includeList=();

    # loop line by line through file
    while (my $line = <$info>) {
        # add included file to includeList (removing seqan3/std/)
        if ($line =~ /\#include <(?:seqan3\/std\/)?(.*)>$/) {
            push(@includeList, $1);

        # Check at end of include block if all included files are in order
        } elsif (scalar @includeList > 0) {
            # check if include block is sorted
            my $isSorted = 1;
            foreach my $i (1 ...$#includeList) {
                if ($includeList[$i-1] gt $includeList[$i]) {
                    $isSorted = 0;
                }
            }

            # print file name if unsorted
            if ($isSorted == 0) {
                print "$file\n";

                # print correct include order
                my @orderedList = sort @includeList;
                foreach my $i (0 ...$#includeList) {
                    print "*" unless ($includeList[$i] eq $orderedList[$i]);
                    print "\t$includeList[$i]";
                    print "\n";
                }

            }
            # clear include block
            @includeList=();
        }
    }

    # close file
    close $info;
}
